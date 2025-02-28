#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openmm import app


def rdkit_to_openmm(rdkit_mol, name="LIG"):
    """
    将RDKit分子转换为OpenMM分子。

    参数
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        需要转换的RDKit分子
    name: str
        分子名称

    返回值
    -------
    omm_molecule: openmm.app.Modeller
        包含目标分子的OpenMM建模器对象
    """
    # 第一步：RDKit转OpenFF
    # 直觉理解：就像将文件从一种格式转换为中间格式
    # 不这么做会怎样：无法利用OpenFF提供的功能和工具
    off_mol = Molecule.from_rdkit(rdkit_mol)

    # 第二步：设置分子名称
    # 直觉理解：给分子贴上标签，方便后续识别
    # 不这么做会怎样：在复杂体系中难以区分不同分子
    off_mol.name = name

    # 第三步：为原子添加名称
    # 直觉理解：给分子中的每个原子编号，像给团队成员分配工号
    # 不这么做会怎样：无法在模拟中追踪和区分同类型的原子
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        # 计数器更新或初始化
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        # 命名格式：元素符号+编号（如C1, C2, O1等）
        off_atom.name = element + str(element_counter_dict[element])

    # 第四步：OpenFF转OpenMM
    # 直觉理解：将中间格式转换为最终需要的格式
    # 不这么做会怎样：无法在OpenMM中使用这个分子进行模拟
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0]

    # 第五步：单位转换（埃转纳米）
    # 直觉理解：将度量单位统一到OpenMM使用的标准
    # 不这么做会怎样：模拟中的距离会差1000倍，导致完全错误的结果
    mol_positions = mol_positions.to("nanometers")

    # TODO 补丁做法，防止名字丢失
    for chain in mol_topology.chains():
        for residue in chain.residues():
            residue.name = name  # 将所有残基的名称都设置为指定的名称
    # 第六步：合并拓扑和位置信息
    # 直觉理解：将分子的"连接关系"和"空间位置"组合在一起
    # 不这么做会怎样：OpenMM无法同时知道原子如何连接以及在哪里
    omm_mol = app.Modeller(mol_topology, mol_positions)

    return omm_mol


def main():
    parser = argparse.ArgumentParser(description="将RDKit分子转换为OpenMM格式")
    parser.add_argument("input", type=Path, help="输入SDF文件路径")
    parser.add_argument("--name", "-n", default="LIG", help="分子名称（默认：LIG）")
    parser.add_argument(
        "--output", "-o", type=Path, required=True, help="输出PDB文件路径"
    )
    args = parser.parse_args()

    try:
        # 读取输入分子
        rdkit_mol = Chem.SDMolSupplier(str(args.input))[0]

        # 转换为OpenMM格式
        omm_mol = rdkit_to_openmm(rdkit_mol, args.name)

        # 保存为PDB文件
        app.PDBFile.writeFile(
            omm_mol.topology, omm_mol.positions, open(str(args.output), "w")
        )

    except Exception as e:
        print(f"错误：{str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
