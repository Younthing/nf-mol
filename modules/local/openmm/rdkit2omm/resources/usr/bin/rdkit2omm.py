#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openmm import app
from collections import defaultdict


def rdkit_to_openmm(rdkit_mol, name="GG2", add_hydrogens=False, conformer_id=0):
    """
    将RDKit分子转换为OpenMM分子。

    参数
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        需要转换的RDKit分子
    name: str
        分子名称
    add_hydrogens: bool
        是否添加氢原子
    conformer_id: int
        使用的构象ID（如果分子有多个构象）

    返回值
    -------
    omm_molecule: openmm.app.Modeller
        包含目标分子的OpenMM建模器对象
    """
    if rdkit_mol is None:
        raise ValueError("提供的RDKit分子为None")

    # 可选添加氢原子
    if add_hydrogens:
        from rdkit import Chem

        rdkit_mol = Chem.AddHs(rdkit_mol)

    # 第一步：RDKit转OpenFF
    off_mol = Molecule.from_rdkit(rdkit_mol)

    # 第二步：设置分子名称
    off_mol.name = name

    # 第三步：为原子添加名称 - 避免添加后缀
    element_counter = defaultdict(int)
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        element_counter[element] += 1
        # 确保名称不超过4个字符（PDB规范）
        atom_number = element_counter[element]
        if atom_number < 100:
            off_atom.name = f"{element}{atom_number}"
        else:
            off_atom.name = f"{element[0]}{atom_number % 100}"

    # 第四步：创建OpenMM拓扑结构
    off_mol_topology = off_mol.to_topology()

    # 直接创建OpenMM拓扑，不尝试修改off_mol_topology中的残基
    mol_topology = off_mol_topology.to_openmm()

    # 检查构象是否存在
    if len(off_mol.conformers) == 0:
        raise ValueError("分子没有构象数据")

    if conformer_id >= len(off_mol.conformers):
        raise ValueError(
            f"构象ID {conformer_id} 超出范围，分子只有 {len(off_mol.conformers)} 个构象"
        )

    mol_positions = off_mol.conformers[conformer_id]

    # 第五步：单位转换（埃转纳米）
    mol_positions = mol_positions.to("nanometers")

    # 第六步：合并拓扑和位置信息
    omm_mol = app.Modeller(mol_topology, mol_positions)

    # 重要：直接修改OpenMM拓扑中的原子名称
    atom_index = 0
    element_counter = defaultdict(int)
    for atom in omm_mol.topology.atoms():
        element = atom.element.symbol
        element_counter[element] += 1
        # 不添加"x"后缀，仅使用元素符号+数字
        atom.name = f"{element}{element_counter[element]}"
        atom_index += 1

    # 修改残基名称
    for chain in omm_mol.topology.chains():
        for residue in chain.residues():
            residue.name = name

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
