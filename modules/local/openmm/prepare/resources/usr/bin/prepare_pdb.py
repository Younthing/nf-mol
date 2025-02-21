#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
蛋白质结构处理CLI工具
适用于生产环境的蛋白质PDB文件预处理工具
"""

import pathlib
import sys
import logging
from typing import Union

import click
import pdbfixer


# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("protein-cli")


def prepare_protein(
    pdb_file: Union[pathlib.Path, str],
    ignore_missing_residues: bool = True,
    ignore_terminal_missing_residues: bool = True,
    ph: float = 7.0,
):
    """
    使用pdbfixer从PDB文件中准备蛋白质。去除诸如配体之类的杂原子，并替换非标准残基。
    为现有残基添加缺失的原子。默认情况下忽略缺失的残基，但可以包含。

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        包含要模拟的系统的PDB文件。
    ignore_missing_residues: bool, optional
        是否应该忽略或构建缺失的残基。
    ignore_terminal_missing_residues: bool, optional
        是否应该忽略或构建链头和链尾缺失的残基。
    ph: float, optional
        用于确定残基质子化状态的pH值。

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        准备好的蛋白质系统。
    """
    logger.info(f"开始处理PDB文件: {pdb_file}")

    fixer = pdbfixer.PDBFixer(str(pdb_file))
    logger.debug("已创建PDBFixer对象")

    fixer.removeHeterogens()
    logger.debug("已移除杂原子")

    fixer.findMissingResidues()
    logger.debug(f"已发现缺失残基: {len(fixer.missingResidues)}")

    # 如果要忽略末端的缺失残基，则将其从字典中删除
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        removed_count = 0
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]
                removed_count += 1
        logger.debug(f"已忽略 {removed_count} 个末端缺失残基")

    # 如果要忽略所有缺失的残基，则清空字典
    if ignore_missing_residues:
        missing_count = len(fixer.missingResidues)
        fixer.missingResidues = {}
        logger.debug(f"已忽略所有 {missing_count} 个缺失残基")

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    logger.debug("已替换非标准残基")

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    logger.debug("已添加缺失原子")

    fixer.addMissingHydrogens(ph)
    logger.debug(f"已在pH {ph}下添加缺失的氢原子")

    logger.info("蛋白质准备完成")
    return fixer


def save_prepared_protein(fixer, output_path: str, file_format: str = "pdb"):
    """保存处理过的蛋白质结构到指定文件"""
    from openmm.app import PDBFile, PDBxFile

    output_path = pathlib.Path(output_path)
    logger.info(f"正在保存处理后的结构到: {output_path}")

    if file_format.lower() == "pdb":
        with open(output_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
    elif file_format.lower() == "cif" or file_format.lower() == "pdbx":
        with open(output_path, "w") as f:
            PDBxFile.writeFile(fixer.topology, fixer.positions, f)
    else:
        raise ValueError(f"不支持的文件格式: {file_format}，支持的格式：pdb, cif, pdbx")

    logger.info(f"文件已保存: {output_path}")


@click.group()
@click.option("--verbose", "-v", is_flag=True, help="启用详细日志输出")
def cli(verbose):
    """蛋白质结构预处理工具 - 用于净化和准备PDB文件进行分子模拟"""
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("已启用详细日志输出")


@cli.command("prepare")
@click.argument("pdb_file", type=click.Path(exists=True))
@click.option("--output", "-o", required=True, help="输出文件路径")
@click.option(
    "--format",
    "-f",
    type=click.Choice(["pdb", "cif", "pdbx"]),
    default="pdb",
    help="输出文件格式 (默认: pdb)",
)
@click.option(
    "--ignore-missing/--include-missing",
    default=True,
    help="是否忽略缺失残基 (默认: 忽略)",
)
@click.option(
    "--ignore-terminal-missing/--include-terminal-missing",
    default=True,
    help="是否忽略链两端的缺失残基 (默认: 忽略)",
)
@click.option(
    "--ph", type=float, default=7.0, help="用于确定残基质子化状态的pH值 (默认: 7.0)"
)
def prepare_command(
    pdb_file, output, format, ignore_missing, ignore_terminal_missing, ph
):
    """准备蛋白质结构文件，清理并修复常见问题"""
    try:
        # 将字符串路径转换为Path对象
        pdb_path = pathlib.Path(pdb_file)

        # 打印处理参数摘要
        logger.info("处理参数:")
        logger.info(f"  输入文件: {pdb_path}")
        logger.info(f"  输出文件: {output} (格式: {format})")
        logger.info(f"  忽略缺失残基: {ignore_missing}")
        logger.info(f"  忽略末端缺失残基: {ignore_terminal_missing}")
        logger.info(f"  pH值: {ph}")

        # 执行蛋白质准备
        fixer = prepare_protein(
            pdb_file=pdb_path,
            ignore_missing_residues=ignore_missing,
            ignore_terminal_missing_residues=ignore_terminal_missing,
            ph=ph,
        )

        # 保存结果
        save_prepared_protein(fixer, output, format)

        logger.info("处理完成!")

    except Exception as e:
        logger.error(f"处理过程中发生错误: {str(e)}")
        import traceback

        logger.debug(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    cli()
