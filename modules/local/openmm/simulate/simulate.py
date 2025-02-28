from pathlib import Path
import sys

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
        self, friction=1.0 / unit.picoseconds, stepsize=2.0 * unit.femtoseconds
    ):
        """设置模拟器"""
        system, modeller = self.prepare_system()
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


# # 使用示例


if __name__ == "__main__":
    pdb = app.PDBFile("/root/learn/nf-mol/xx.pdb")

    rdkit_ligand = Chem.SDMolSupplier(
        "/root/learn/nf-mol/data/openmm/6w70_prepared.sdf"
    )[0]

    md_sim = MolecularDynamics(
        complex_topology=pdb.topology,
        complex_positions=pdb.positions,
        rdkit_ligand=rdkit_ligand,
        output_dir="simulation_results",
    )

    # 运行完整模拟
    md_sim.run_simulation(steps=5000)

    # # 或者分步运行
    # md_sim.setup_simulation()
    # md_sim.minimize_energy()
    # md_sim.save_topology()
    # md_sim.run_simulation()
