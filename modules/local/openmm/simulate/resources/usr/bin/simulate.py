from pathlib import Path
import sys
import argparse

from rdkit import Chem
import mdtraj as md
import openmm as mm
import openmm.app as app
from openmm import unit
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator


class MolecularDynamics:
    """分子动力学模拟的主类"""

    def __init__(
        self,
        complex_topology,
        complex_positions,
        rdkit_ligand=None,
        output_dir=None,
        temperature=300 * unit.kelvin,
        protein_ff="amber14-all.xml",
        solvent_ff="amber14/tip3pfb.xml",
    ):
        self.complex_topology = complex_topology
        self.complex_positions = complex_positions
        self.rdkit_ligand = rdkit_ligand
        self.output_dir = Path(output_dir) if output_dir else Path("./output")
        self.temperature = temperature
        self.protein_ff = protein_ff
        self.solvent_ff = solvent_ff
        self.simulation = None

    def setup_forcefield(self):
        """设置力场"""
        forcefield = app.ForceField(self.protein_ff, self.solvent_ff)

        if self.rdkit_ligand is not None:
            gaff = GAFFTemplateGenerator(
                molecules=Molecule.from_rdkit(
                    self.rdkit_ligand, allow_undefined_stereo=True
                )
            )
            forcefield.registerTemplateGenerator(gaff.generator)

        return forcefield

    def prepare_system(
        self, padding=1.0 * unit.nanometers, ionic_strength=0.15 * unit.molar
    ):
        """准备模拟系统"""
        forcefield = self.setup_forcefield()
        modeller = app.Modeller(self.complex_topology, self.complex_positions)
        modeller.addSolvent(
            forcefield,
            padding=padding,
            ionicStrength=ionic_strength,
        )

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)

        return system, modeller

    def setup_simulation(
        self,
        friction=1.0 / unit.picoseconds,
        stepsize=2.0 * unit.femtoseconds,
        padding=1.0 * unit.nanometers,
        ionic_strength=0.15 * unit.molar,
    ):
        """设置模拟器"""
        system, modeller = self.prepare_system(
            padding=padding, ionic_strength=ionic_strength
        )
        integrator = mm.LangevinIntegrator(self.temperature, friction, stepsize)
        self.simulation = app.Simulation(modeller.topology, system, integrator)
        self.simulation.context.setPositions(modeller.positions)
        return self.simulation

    def minimize_energy(self):
        """能量最小化"""
        if self.simulation is None:
            self.setup_simulation()
        self.simulation.minimizeEnergy()

    def save_topology(self):
        """保存拓扑结构"""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        topology_file = self.output_dir / "topology.pdb"

        with open(topology_file, "w", encoding="utf-8") as pdb_file:
            app.PDBFile.writeFile(
                self.simulation.topology,
                self.simulation.context.getState(
                    getPositions=True, enforcePeriodicBox=True
                ).getPositions(),
                file=pdb_file,
                keepIds=True,
            )

    def add_reporters(self, steps, write_interval=10, log_interval=250):
        """添加报告器"""
        self.simulation.reporters.append(
            md.reporters.XTCReporter(
                file=str(self.output_dir / "trajectory.xtc"),
                reportInterval=write_interval,
            )
        )
        self.simulation.reporters.append(
            app.StateDataReporter(
                sys.stdout,
                log_interval,
                step=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=steps,
                separator="\t",
            )
        )

    def run_simulation(self, steps=5000, write_interval=10, log_interval=250):
        """运行模拟"""
        if self.simulation is None:
            self.setup_simulation()

        self.minimize_energy()
        self.save_topology()
        self.add_reporters(steps, write_interval, log_interval)

        self.simulation.context.setVelocitiesToTemperature(self.temperature)
        self.simulation.step(steps)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="分子动力学模拟程序",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 必需参数
    parser.add_argument("--pdb", required=True, help="蛋白质复合物PDB文件路径")

    # 可选参数
    parser.add_argument("--ligand", help="配体分子SDF文件路径")
    parser.add_argument("--output", default="./output", help="输出目录路径")
    parser.add_argument("--steps", type=int, default=5000, help="模拟步数")
    parser.add_argument("--temperature", type=float, default=300.0, help="模拟温度(K)")
    parser.add_argument("--write-interval", type=int, default=10, help="轨迹写入间隔")
    parser.add_argument("--log-interval", type=int, default=250, help="日志输出间隔")
    parser.add_argument("--protein-ff", default="amber14-all.xml", help="蛋白质力场")
    parser.add_argument("--solvent-ff", default="amber14/tip3pfb.xml", help="溶剂力场")
    parser.add_argument(
        "--padding", type=float, default=1.0, help="溶剂盒子填充距离(nm)"
    )
    parser.add_argument(
        "--ionic-strength", type=float, default=0.15, help="离子强度(M)"
    )
    parser.add_argument(
        "--friction", type=float, default=1.0, help="朗之万摩擦系数(1/ps)"
    )
    parser.add_argument("--stepsize", type=float, default=2.0, help="时间步长(fs)")

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 加载PDB文件
    pdb = app.PDBFile(args.pdb)

    # 加载配体分子(如果提供)
    rdkit_ligand = None
    if args.ligand:
        rdkit_ligand = Chem.SDMolSupplier(args.ligand)[0]
        if rdkit_ligand is None:
            print(f"错误: 无法加载配体文件 {args.ligand}", file=sys.stderr)
            sys.exit(1)

    # 设置单位
    temperature = args.temperature * unit.kelvin
    padding = args.padding * unit.nanometers
    ionic_strength = args.ionic_strength * unit.molar
    friction = args.friction / unit.picoseconds
    stepsize = args.stepsize * unit.femtoseconds

    # 创建分子动力学实例
    md_sim = MolecularDynamics(
        complex_topology=pdb.topology,
        complex_positions=pdb.positions,
        rdkit_ligand=rdkit_ligand,
        output_dir=args.output,
        temperature=temperature,
        protein_ff=args.protein_ff,
        solvent_ff=args.solvent_ff,
    )

    # 设置模拟参数并运行
    md_sim.setup_simulation(
        friction=friction,
        stepsize=stepsize,
        padding=padding,  # 添加这一行
        ionic_strength=ionic_strength,  # 添加这一行
    )
    md_sim.run_simulation(
        steps=args.steps,
        write_interval=args.write_interval,
        log_interval=args.log_interval,
    )

    print(f"模拟完成! 结果保存在 {args.output} 目录")


if __name__ == "__main__":
    main()

#  python modules/local/openmm/simulate/resources/usr/bin/simulate.py --pdb data/openmm/6w70_complex.pdb --ligand data/openmm/6w70_prepared.sdf --output sim_results --steps 10000 --temperature 310
