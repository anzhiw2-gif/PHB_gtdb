#!/usr/bin/env python3
"""
PHB_gtdb — Step 4: 系统发育分析

基于 PHB 基因比对结果构建系统发育树。
- 基因树：PHB 降解/合成酶的进化关系
- 物种树：已有 GTDB bac120 树 (R226)

支持方法:
- IQ-TREE: 最大似然法 (ML)，自动模型选择
- FastTree: 快速近似 ML (备用)

用法:
    python 04_phylogeny.py [--gene phaZ] [--method iqtree]
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional

import numpy as np
import pandas as pd
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import setup_logging, run_cmd, read_fasta, count_sequences

PHAZ_SUBTYPES = [
    "phaz_bacillus_type",
    "phaz_extracellular",
    "phaz_extracellular_lemoignei",
    "phaz_intracellular",
    "phaz_ralstonia",
]


def run_iqtree(alignment: Path, output_prefix: str, threads: int = 30,
               bootstrap: int = 1000, model: str = "MFP",
               logger: logging.Logger = None) -> Path:
    """IQ-TREE 最大似然系统发育树。

    Key features:
    - ModelFinder Plus: 自动选择最优氨基酸替换模型
    - UFBoot: 超快 Bootstrap
    - SH-aLRT: 分支检验
    """
    cmd = (
        f"iqtree -s {alignment} "
        f"-pre {output_prefix} "
        f"-m {model} "
        f"-B {bootstrap} "
        f"-T {threads} "
        f"-alrt 1000 "
        f"--keep-ident "
        f"--quiet"
    )

    logger.info(f"IQ-TREE (bootstrap={bootstrap}, model={model})...")
    logger.info("  此步骤可能耗时较长，取决于序列数量和长度...")

    run_cmd(cmd, f"IQ-TREE: {alignment.name}", timeout=172800, logger=logger)

    tree_file = Path(f"{output_prefix}.treefile")
    if not tree_file.exists():
        raise RuntimeError(f"IQ-TREE 输出不存在: {tree_file}")

    # 读取并报告结果
    tree = Phylo.read(tree_file, "newick")
    n_leaves = len(list(tree.get_terminals()))
    logger.info(f"IQ-TREE 完成: {n_leaves} 个叶节点")

    # 报告最优模型
    iqtree_log = Path(f"{output_prefix}.log")
    if iqtree_log.exists():
        with open(iqtree_log) as f:
            for line in f:
                if "Best-fit model" in line:
                    logger.info(f"最优模型: {line.strip()}")
                    break

    return tree_file


def run_fasttree(alignment: Path, output_tree: Path, threads: int = 30,
                 logger: logging.Logger = None) -> Path:
    """FastTree 快速近似 ML 树 (备用方案，适合超大序列集)。"""
    cmd = (
        f"fasttree -lg -gamma -spr 4 -mlacc 2 "
        f"-thread {threads} "
        f"{alignment} > {output_tree}"
    )

    logger.info(f"FastTree (LG+gamma)...")
    run_cmd(cmd, f"FastTree: {alignment.name}", timeout=86400, logger=logger)

    return output_tree


def compute_bootstrap_stats(tree_file: Path) -> Dict:
    """统计 Bootstrap 值分布。"""
    tree = Phylo.read(tree_file, "newick")
    bootstrap_values = []
    for clade in tree.find_clades():
        if clade.confidence is not None:
            bootstrap_values.append(clade.confidence)

    if not bootstrap_values:
        return {"n_branches": 0}

    return {
        "n_branches": len(bootstrap_values),
        "mean": np.mean(bootstrap_values),
        "median": np.median(bootstrap_values),
        "min": np.min(bootstrap_values),
        "q25": np.percentile(bootstrap_values, 25),
        "q75": np.percentile(bootstrap_values, 75),
        "above_70": sum(1 for v in bootstrap_values if v >= 70),
        "above_90": sum(1 for v in bootstrap_values if v >= 90),
    }


def build_genome_phb_matrix(tree_file: Path, phb_counts: Path) -> Path:
    """将 PHB 基因计数映射到 GTDB 物种树的叶节点。

    用于后续在系统发育树上可视化 PHB 基因分布。
    """
    # 读取 PHB 计数
    df = pd.read_csv(phb_counts, sep="\t")
    gene_cols = [c for c in df.columns if c != "genome_id"]

    # 读取物种树
    tree = Phylo.read(tree_file, "newick")
    leaf_ids = {leaf.name for leaf in tree.get_terminals()}
    # GTDB 树使用 RS_GCF_... 或 GB_GCA_... 格式
    # 需要映射

    # 构建基因组 → PHB 计数映射
    genome_to_phb = {}
    for _, row in df.iterrows():
        gid = row["genome_id"]
        counts = {col: int(row[col]) for col in gene_cols}
        genome_to_phb[gid] = counts
        # 也尝试 RS_ 和 GB_ 前缀
        genome_to_phb[f"RS_{gid}"] = counts
        genome_to_phb[f"GB_{gid}"] = counts

    # 统计树中有多少基因组有 PHB 基因
    n_in_tree = 0
    n_with_phb = 0
    for leaf_id in leaf_ids:
        n_in_tree += 1
        if leaf_id in genome_to_phb:
            n_with_phb += 1

    return genome_to_phb, n_in_tree, n_with_phb


def main():
    parser = argparse.ArgumentParser(description="系统发育分析")
    parser.add_argument("--gene", type=str, default="all",
                        help="指定基因 (phaZ/phaC/all)")
    parser.add_argument("--method", type=str, default="iqtree",
                        choices=["iqtree", "fasttree"],
                        help="树构建方法 (默认: iqtree)")
    parser.add_argument("--threads", type=int, default=PHYLO_THREADS,
                        help=f"线程数 (默认: {PHYLO_THREADS})")
    parser.add_argument("--bootstrap", type=int, default=PHYLO_BOOTSTRAP,
                        help=f"Bootstrap 数 (默认: {PHYLO_BOOTSTRAP})")
    parser.add_argument("--input", type=str, default=None,
                        help="自定义比对文件")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "04_phylogeny")

    # --- 基因树 ---
    if args.input:
        alignment_files = {"custom": Path(args.input)}
    elif args.gene == "all":
        alignment_files = {}
        for subtype in PHAZ_SUBTYPES:
            f = PROCESSED_DIR / f"{subtype}_trim.fasta"
            if not f.exists():
                f = PROCESSED_DIR / f"{subtype}_msa.fasta"
            if f.exists() and count_sequences(f) >= 3:
                alignment_files[subtype] = f
        if not alignment_files:
            for gene_name in PHB_GENES:
                f = PROCESSED_DIR / f"{gene_name}_trim.fasta"
                if not f.exists():
                    f = PROCESSED_DIR / f"{gene_name}_trimmed.fasta"
                if not f.exists():
                    f = PROCESSED_DIR / f"{gene_name}_msa.fasta"
                if not f.exists():
                    f = PROCESSED_DIR / f"{gene_name}_aligned.fasta"
                if f.exists() and count_sequences(f) >= 3:
                    alignment_files[gene_name] = f
    else:
        candidates = [
            PROCESSED_DIR / f"{args.gene}_trim.fasta",
            PROCESSED_DIR / f"{args.gene}_trimmed.fasta",
            PROCESSED_DIR / f"{args.gene}_msa.fasta",
            PROCESSED_DIR / f"{args.gene}_aligned.fasta",
            PROCESSED_DIR / f"phaz_{args.gene}_trim.fasta",
            PROCESSED_DIR / f"phaz_{args.gene}_msa.fasta",
        ]
        f = next((p for p in candidates if p.exists()), None)
        if f is None:
            logger.error(f"比对文件不存在: {args.gene}")
            sys.exit(1)
        alignment_files = {args.gene: f}

    if not alignment_files:
        logger.error("未找到比对文件，请先运行 03_msa.py")
        sys.exit(1)

    logger.info(f"待建树: {list(alignment_files.keys())}")

    for gene_name, aln_file in alignment_files.items():
        logger.info("=" * 60)
        logger.info(f"构建 {gene_name} 基因树")
        logger.info("=" * 60)

        n_seqs = count_sequences(aln_file)
        logger.info(f"序列数: {n_seqs}")

        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        tree_prefix = str(RESULTS_DIR / f"{gene_name}_tree")

        if args.method == "iqtree":
            # 检查 IQ-TREE 是否可用
            iqtree_available = os.system("which iqtree > /dev/null 2>&1") == 0
            if not iqtree_available:
                # 尝试 conda
                iqtree_path = "/home/data/haoyu/miniconda3/bin/iqtree"
                if os.path.exists(iqtree_path):
                    os.environ["PATH"] += f":{os.path.dirname(iqtree_path)}"
                    iqtree_available = True

            if iqtree_available:
                tree_file = run_iqtree(
                    aln_file, tree_prefix, args.threads,
                    args.bootstrap, logger=logger
                )
            else:
                logger.warning("IQ-TREE 未安装，使用 FastTree 替代")
                tree_file = Path(f"{tree_prefix}.treefile")
                run_fasttree(aln_file, tree_file, args.threads, logger)
        else:
            tree_file = Path(f"{tree_prefix}.treefile")
            run_fasttree(aln_file, tree_file, args.threads, logger)

        # Bootstrap 统计
        if tree_file.exists():
            stats = compute_bootstrap_stats(tree_file)
            if stats.get("n_branches", 0) > 0:
                logger.info(
                    f"Bootstrap: mean={stats['mean']:.1f}, "
                    f"median={stats['median']:.1f}, "
                    f">=70: {stats['above_70']}/{stats['n_branches']} "
                    f"({100*stats['above_70']/stats['n_branches']:.1f}%)"
                )

    # --- 物种树映射 ---
    logger.info("=" * 60)
    logger.info("基因-物种树映射")
    logger.info("=" * 60)

    if GTDB_TREE.exists():
        search_result = PROCESSED_DIR / "phb_search_results.tsv"
        if search_result.exists():
            mapping, n_in_tree, n_with_phb = build_genome_phb_matrix(
                GTDB_TREE, search_result
            )
            logger.info(
                f"GTDB 物种树: {n_in_tree} 个分类单元, "
                f"其中 {n_with_phb} 个含 PHB 基因"
            )
    else:
        logger.warning(f"GTDB 物种树不存在: {GTDB_TREE}")

    logger.info("=" * 60)
    logger.info("Step 4 完成!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
