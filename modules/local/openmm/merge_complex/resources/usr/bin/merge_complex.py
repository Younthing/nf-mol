#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from openmm.app import PDBFile
import mdtraj as md
from openmm import unit
import numpy as np


def merge_protein_and_ligand(protein_pdb, ligand_pdb):
    """
    合并蛋白质和配体的PDB结构。

    参数说明
    --------
    protein_pdb: openmm.app.PDBFile
        蛋白质PDB文件对象
    ligand_pdb: openmm.app.PDBFile
        配体PDB文件对象

    返回值
    -------
    tuple:
        - complex_topology: openmm.app.topology.Topology
            合并后的拓扑结构
        - complex_positions: openmm.unit.quantity.Quantity
            合并后的原子坐标（单位：纳米）
    """
    # 使用mdtraj转换拓扑结构
    md_protein_topology = md.Topology.from_openmm(protein_pdb.topology)
    md_ligand_topology = md.Topology.from_openmm(ligand_pdb.topology)

    # 验证输入结构
    if md_protein_topology.n_atoms == 0:
        raise ValueError("蛋白质结构不包含任何原子")
    if md_ligand_topology.n_atoms == 0:
        raise ValueError("配体结构不包含任何原子")

    # 合并拓扑结构
    md_complex_topology = md_protein_topology.join(md_ligand_topology)
    complex_topology = md_complex_topology.to_openmm()

    # 合并原子坐标
    total_atoms = len(protein_pdb.positions) + len(ligand_pdb.positions)
    complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
    complex_positions[: len(protein_pdb.positions)] = protein_pdb.positions
    complex_positions[len(protein_pdb.positions) :] = ligand_pdb.positions

    return complex_topology, complex_positions


def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="合并蛋白质和配体结构文件")
    parser.add_argument("protein", type=Path, help="输入蛋白质PDB文件路径")
    parser.add_argument("ligand", type=Path, help="输入配体PDB文件路径")
    parser.add_argument(
        "--output", "-o", type=Path, required=True, help="输出复合物PDB文件路径"
    )
    args = parser.parse_args()

    try:
        # 验证输入文件是否存在
        if not args.protein.exists():
            raise FileNotFoundError(f"找不到蛋白质文件：{args.protein}")
        if not args.ligand.exists():
            raise FileNotFoundError(f"找不到配体文件：{args.ligand}")

        # 直接使用PDBFile读取结构
        protein_pdb = PDBFile(str(args.protein))
        ligand_pdb = PDBFile(str(args.ligand))

        # 合并结构
        complex_topology, complex_positions = merge_protein_and_ligand(
            protein_pdb, ligand_pdb
        )

        # 保存合并后的复合物结构
        with open(args.output, "w") as f:
            PDBFile.writeFile(complex_topology, complex_positions, f)
            print(f"成功将复合物结构保存到：{args.output}")

    except Exception as e:
        print(f"错误：{str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
