#!/usr/bin/env python3
import argparse
import sys
import copy
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


def prepare_ligand(sdf_file, smiles, save=None):
    """
    从SDF文件准备配体结构，包括添加氢原子和分配键级。
    可以将配体的准备前后状态保存为图片文件以便检查结果。

    参数说明
    ----------
    sdf_file: pathlib.PosixPath
        包含目标配体的SDF文件
    smiles: str
        配体的SMILES字符串，用于提供正确的质子化状态和键级信息
    save: str or pathlib.Path, 可选
        如果提供，将比对图保存到指定路径，例如 'ligand_comparison.png'

    返回值
    -------
    prepared_ligand: rdkit.Chem.rdchem.Mol
        准备好的配体分子对象
    """
    ligand = Chem.SDMolSupplier(sdf_file)[0]  # 读取sdf中第一个分子

    ligand = Chem.RemoveHs(ligand)

    # 第四步：使用模板分配键级
    # 直觉理解：用标准图纸（SMILES）来修正模型中的连接关系
    # 不这么做会怎样：分子中的化学键性质将不准确，影响后续所有依赖于电子结构的计算
    reference_mol = Chem.MolFromSmiles(smiles)
    prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference_mol, ligand)
    # 保持原有的3D结构
    # 直觉理解：保留分子的"姿势"，只改变内部连接方式
    prepared_ligand.AddConformer(ligand.GetConformer(0))

    # 第五步：添加氢原子
    # 直觉理解：在正确的位置添加上所有缺失的氢原子
    # 不这么做会怎样：分子的价键将不完整，影响分子动力学模拟的准确性
    prepared_ligand = Chem.rdmolops.AddHs(prepared_ligand, addCoords=True)
    prepared_ligand = Chem.MolFromMolBlock(Chem.MolToMolBlock(prepared_ligand))

    # 第六步：保存比对图（可选）
    # 直觉理解：将处理前后的对比图保存到文件中，方便后续查看和分析
    # 不这么做会怎样：无法保留处理结果的可视化记录
    if save:
        save_path = Path(save)
        ligand_2d = copy.deepcopy(ligand)
        prepared_ligand_2d = copy.deepcopy(prepared_ligand)
        AllChem.Compute2DCoords(ligand_2d)
        AllChem.Compute2DCoords(prepared_ligand_2d)
        img = Draw.MolsToGridImage(
            [ligand_2d, prepared_ligand_2d],
            molsPerRow=2,
            legends=["before", "after"],
            returnPNG=False,
        )
        img.save(save_path)

    return prepared_ligand

def main():
    parser = argparse.ArgumentParser(description='准备配体SDF文件，包括添加氢原子和分配键级')
    parser.add_argument('input', type=Path, help='输入SDF文件路径')
    parser.add_argument('--smiles', '-s', required=True, help='配体的SMILES字符串')
    parser.add_argument('--output', '-o', type=Path, required=True, help='输出SDF文件路径')
    parser.add_argument('--image', '-i', type=Path, help='可选：保存处理前后的对比图')
    args = parser.parse_args()

    try:
        # 准备配体
        prepared_mol = prepare_ligand(args.input, args.smiles, args.image)
        
        # 保存结果
        writer = Chem.SDWriter(str(args.output))
        writer.write(prepared_mol)
        writer.close()
        
    except Exception as e:
        print(f"错误：{str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
