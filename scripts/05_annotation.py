#!/usr/bin/env python3
"""
PHB_gtdb — Step 5: 功能注释

对鉴定到的 PHB 基因进行功能注释:
1. eggNOG-mapper: KEGG/COG/GO 注释
2. Pfam 结构域注释
3. 信号肽预测 (SignalP)
4. PHB 通路完整性评估

用法:
    python 05_annotation.py
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import (setup_logging, run_cmd, read_fasta, write_fasta,
                    load_gtdb_taxonomy, parse_gtdb_taxonomy)


def run_eggnog_mapper(input_fasta: Path, output_dir: Path,
                      threads: int = 30,
                      logger: logging.Logger = None) -> pd.DataFrame:
    """eggNOG-mapper v2 功能注释。

    输出包含: KEGG_ko, KEGG_Pathway, COG_category, GO_terms, PFAMs, EC
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / "eggnog"
    output_file = Path(str(output_prefix) + ".emapper.annotations")

    if output_file.exists():
        logger.info("eggNOG 结果已存在，跳过")
        return pd.read_csv(output_file, sep="\t", comment="#")

    # 检查 eggnog-mapper 是否可用
    emapper_path = None
    for p in [
        "/home/data/haoyu/miniconda3/bin/emapper.py",
        "/usr/local/bin/emapper.py",
        "emapper.py",
    ]:
        if os.path.exists(p) or os.system(f"which {p} > /dev/null 2>&1") == 0:
            emapper_path = p
            break

    if emapper_path is None:
        logger.warning("eggNOG-mapper 未安装，跳过此步骤")
        logger.info("安装方法: conda install -c bioconda eggnog-mapper")
        return pd.DataFrame()

    cmd = (
        f"{emapper_path} "
        f"-i {input_fasta} "
        f"--output {output_prefix} "
        f"--cpu {threads} "
        f"--database bact_arch "
        f"--override "
        f"--report_orthologs "
        f"--report_unannotated"
    )

    logger.info("运行 eggNOG-mapper (预计 1-6 小时)...")
    try:
        run_cmd(cmd, "eggNOG-mapper", timeout=86400, logger=logger)
    except Exception as e:
        logger.error(f"eggNOG-mapper 失败: {e}")
        return pd.DataFrame()

    if output_file.exists():
        df = pd.read_csv(output_file, sep="\t", comment="#")
        logger.info(f"eggNOG 完成: {len(df)} 条注释")
        return df
    return pd.DataFrame()


def parse_kegg_from_eggnog(eggnog_df: pd.DataFrame) -> Dict[str, List[str]]:
    """从 eggNOG 结果中提取 KEGG KO 信息。

    Returns:
        {gene_id: [KO1, KO2, ...]} 字典
    """
    gene_kos = defaultdict(list)

    ko_col = None
    for col in ["KEGG_ko", "KEGG_KOs"]:
        if col in eggnog_df.columns:
            ko_col = col
            break

    if ko_col:
        for _, row in eggnog_df.iterrows():
            if pd.notna(row[ko_col]):
                kos = str(row[ko_col]).split(",")
                gene_id = str(row.get("#query", row.get("query", "")))
                gene_kos[gene_id] = [k.strip() for k in kos if k.strip() and k.strip() != "-"]

    return dict(gene_kos)


def parse_go_from_eggnog(eggnog_df: pd.DataFrame) -> Dict[str, List[str]]:
    """从 eggNOG 结果提取 GO terms。"""
    gene_gos = defaultdict(list)

    for col in ["GO_terms", "GOs"]:
        if col in eggnog_df.columns:
            for _, row in eggnog_df.iterrows():
                if pd.notna(row[col]):
                    gos = str(row[col]).split(",")
                    gene_id = str(row.get("#query", row.get("query", "")))
                    gene_gos[gene_id] = [g.strip() for g in gos if g.strip() and g.strip() != "-"]
            break

    return dict(gene_gos)


def annotate_phb_pathway(eggnog_df: pd.DataFrame,
                         logger: logging.Logger = None) -> pd.DataFrame:
    """评估 PHB 代谢通路的完整性。

    完整的 PHB 降解通路:
        phaZ: PHB → 3-hydroxybutyrate (K01069)
        phaC: acetyl-CoA → PHB (K03821)
        phaA: acetyl-CoA + acetyl-CoA → acetoacetyl-CoA (K00626)
        phaB: acetoacetyl-CoA → 3-hydroxybutyryl-CoA (K00023)

    对每个基因组评估其 PHB 通路完整性。
    """
    # PHB 通路关键 KO
    phb_pathway_kos = {
        "phaZ": ["K01069"],  # PHB depolymerase
        "phaC": ["K03821"],  # PHA synthase
        "phaA": ["K00626"],  # beta-ketothiolase
        "phaB": ["K00023"],  # acetoacetyl-CoA reductase
    }

    # 提取 KO → gene 映射
    gene_kos = parse_kegg_from_eggnog(eggnog_df)

    # 按基因组分组
    genome_ko_counts = defaultdict(lambda: defaultdict(int))
    for gene_id, kos in gene_kos.items():
        genome_id = gene_id.split("|")[0] if "|" in gene_id else gene_id
        for ko in kos:
            for phb_gene, ko_list in phb_pathway_kos.items():
                if ko in ko_list:
                    genome_ko_counts[genome_id][phb_gene] += 1

    if not genome_ko_counts:
        if logger:
            logger.warning("未找到 PHB 通路 KO 注释")
        return pd.DataFrame()

    # 构建结果表
    rows = []
    for genome_id, counts in genome_ko_counts.items():
        row = {"genome_id": genome_id}
        for gene in phb_pathway_kos:
            row[gene] = counts.get(gene, 0)
        row["pathway_completeness"] = sum(1 for g in phb_pathway_kos if row[g] > 0)
        row["has_phaZ"] = row["phaZ"] > 0
        row["has_phaC"] = row["phaC"] > 0
        row["has_degradation"] = row["has_phaZ"]
        row["has_synthesis"] = all(row[g] > 0 for g in ["phaC", "phaA", "phaB"])
        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values("pathway_completeness", ascending=False)

    if logger:
        logger.info(f"PHB 通路分析: {len(df)} 个基因组")
        logger.info(f"  含降解通路 (phaZ): {df['has_degradation'].sum()}")
        logger.info(f"  含合成通路 (phaC+phaA+phaB): {df['has_synthesis'].sum()}")

    return df


def phylum_level_summary(phb_counts: Path, taxonomy_df: pd.DataFrame,
                         output_path: Path) -> pd.DataFrame:
    """按门水平汇总 PHB 基因分布。"""
    df = pd.read_csv(phb_counts, sep="\t")
    gene_cols = [c for c in df.columns if c.endswith("_count")]

    # 添加分类学
    taxonomy_df["genome_id"] = taxonomy_df["accession"].apply(
        lambda x: str(x).replace("RS_", "").replace("GB_", "")
    )

    df = df.merge(taxonomy_df[["genome_id", "taxonomy"]], on="genome_id", how="left")

    # 解析门
    def get_phylum(tax_str):
        if pd.isna(tax_str):
            return "Unknown"
        tax = parse_gtdb_taxonomy(str(tax_str))
        return tax.get("phylum", "Unknown")

    df["phylum"] = df["taxonomy"].apply(get_phylum)

    # 按门汇总
    summary = df.groupby("phylum").agg(
        n_genomes=("genome_id", "nunique"),
        **{f"total_{col}": (col, "sum") for col in gene_cols},
        **{f"genomes_with_{col}": (col, lambda x: (x > 0).sum()) for col in gene_cols},
    ).reset_index()

    summary = summary.sort_values("n_genomes", ascending=False)
    summary.to_csv(output_path, sep="\t", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser(description="功能注释")
    parser.add_argument("--eggnog", action="store_true", default=True,
                        help="运行 eggNOG-mapper 注释")
    parser.add_argument("--no-eggnog", action="store_true",
                        help="跳过 eggNOG-mapper")
    parser.add_argument("--threads", type=int, default=MAX_THREADS,
                        help=f"线程数 (默认: {MAX_THREADS})")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "05_annotation")

    # 输入: PHB 蛋白序列 (去冗余后)
    input_candidates = [
        PROCESSED_DIR / "phaz_proteins_validated.fasta",
        PROCESSED_DIR / "phaz_proteins_filtered.fasta",
        PROCESSED_DIR / "phaz_proteins_c95.fasta",
        PROCESSED_DIR / "phaz_proteins_all.fasta",
        PROCESSED_DIR / "phb_proteins_dedup.fasta",
        PROCESSED_DIR / "phb_proteins_annotated.fasta",
    ]
    input_fasta = next((p for p in input_candidates if p.exists() and p.stat().st_size > 0), None)
    if input_fasta is None:
        logger.error("PHB 蛋白序列不存在，请先运行 02_extract_sequences.py")
        sys.exit(1)

    n_seqs = sum(1 for _ in open(input_fasta))
    logger.info(f"待注释序列: {n_seqs}")

    eggnog_df = pd.DataFrame()

    # --- 1. eggNOG-mapper ---
    if not args.no_eggnog:
        logger.info("=" * 60)
        logger.info("Step 1: eggNOG-mapper 功能注释")
        logger.info("=" * 60)

        eggnog_dir = PROCESSED_DIR / "annotations"
        eggnog_df = run_eggnog_mapper(input_fasta, eggnog_dir, args.threads, logger)

    # --- 2. PHB 通路分析 ---
    logger.info("=" * 60)
    logger.info("Step 2: PHB 代谢通路完整性分析")
    logger.info("=" * 60)

    if not eggnog_df.empty:
        pathway_df = annotate_phb_pathway(eggnog_df, logger)

        if not pathway_df.empty:
            pathway_out = TABLES_DIR / "phb_pathway_completeness.tsv"
            pathway_df.to_csv(pathway_out, sep="\t", index=False)
            logger.info(f"PHB 通路分析保存至: {pathway_out}")

    # --- 3. 门水平分布 ---
    logger.info("=" * 60)
    logger.info("Step 3: 门水平 PHB 基因分布")
    logger.info("=" * 60)

    search_result = PROCESSED_DIR / "phb_search_results.tsv"
    if search_result.exists():
        taxonomy = load_gtdb_taxonomy()
        phylum_summary = phylum_level_summary(
            search_result, taxonomy,
            TABLES_DIR / "phb_phylum_summary.tsv"
        )

        if not phylum_summary.empty:
            logger.info(f"门水平分布 (前10):")
            for _, row in phylum_summary.head(10).iterrows():
                logger.info(f"  {row['phylum']}: {int(row['n_genomes'])} 个基因组")

    logger.info("=" * 60)
    logger.info("Step 5 完成!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
