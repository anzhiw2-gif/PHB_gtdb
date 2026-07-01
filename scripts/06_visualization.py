#!/usr/bin/env python3
"""
PHB_gtdb — Step 6: 可视化

生成论文级别的图表:
1. 系统发育树 + PHB 基因分布热图
2. 门水平 PHB 基因分布 (柱状图)
3. PHB 代谢通路完整性热图
4. PHB 基因丰度与基因组大小相关性
5. PhaZ 序列保守性 Logo (需要 WebLogo)

用法:
    python 06_visualization.py
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互后端
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import setup_logging, read_fasta, load_gtdb_taxonomy, parse_gtdb_taxonomy


# 设置中文字体 (如果可用)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "figure.dpi": FIGURE_DPI,
})


def plot_phylum_distribution(phylum_summary: Path, output_path: Path,
                             top_n: int = 15,
                             logger: logging.Logger = None):
    """门水平 PHB 基因分布柱状图。"""
    df = pd.read_csv(phylum_summary, sep="\t")

    # 取 Top N
    df_plot = df.head(top_n).copy()
    df_plot = df_plot.sort_values("n_genomes", ascending=True)

    fig, ax = plt.subplots(figsize=(12, 8))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(df_plot)))
    bars = ax.barh(range(len(df_plot)), df_plot["n_genomes"], color=colors, edgecolor="white")

    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot["phylum"])
    ax.set_xlabel("Number of Genomes with PHB Genes")
    ax.set_title(f"Distribution of PHB-Related Genes Across Phyla (Top {top_n})")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # 添加数值标签
    for i, (bar, val) in enumerate(zip(bars, df_plot["n_genomes"])):
        ax.text(bar.get_width() + bar.get_width() * 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{int(val):,}", va="center", fontsize=10)

    plt.tight_layout()
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    if logger:
        logger.info(f"门分布图: {output_path}")


def plot_phb_gene_counts_heatmap(phb_counts: Path, taxonomy_df,
                                  output_path: Path,
                                  top_n: int = 50,
                                  logger: logging.Logger = None):
    """PHB 基因分布热图 (Top 基因组)。"""
    df = pd.read_csv(phb_counts, sep="\t")
    gene_cols = [c for c in df.columns if c.endswith("_count")]

    if not gene_cols:
        return

    # 计算总 PHB 基因数
    df["total"] = df[gene_cols].sum(axis=1)
    df = df.nlargest(top_n, "total")

    # 准备热图数据
    heatmap_data = df.set_index("genome_id")[gene_cols]

    # 缩短基因组 ID
    heatmap_data.index = [x[:20] for x in heatmap_data.index]

    fig, ax = plt.subplots(figsize=(10, max(8, top_n * 0.25)))

    sns.heatmap(heatmap_data, cmap="YlOrRd", annot=False,
                cbar_kws={"label": "Gene Count"}, ax=ax,
                linewidths=0.5, linecolor="gray")

    ax.set_title(f"PHB Gene Distribution (Top {top_n} Genomes)")
    ax.set_xlabel("PHB-Related Genes")
    ax.set_ylabel("Genome")

    plt.tight_layout()
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    if logger:
        logger.info(f"基因热图: {output_path}")


def plot_pathway_completeness(pathway_file: Path, output_path: Path,
                               logger: logging.Logger = None):
    """PHB 代谢通路完整性饼图。"""
    if not pathway_file.exists():
        if logger:
            logger.warning(f"通路文件不存在: {pathway_file}")
        return

    df = pd.read_csv(pathway_file, sep="\t")

    # 分类
    categories = {
        "Both (D+S)": (df["has_degradation"] & df["has_synthesis"]).sum(),
        "Degradation only (phaZ)": (df["has_degradation"] & ~df["has_synthesis"]).sum(),
        "Synthesis only (phaCAB)": (~df["has_degradation"] & df["has_synthesis"]).sum(),
        "Partial (其他)": (~df["has_degradation"] & ~df["has_synthesis"]).sum(),
    }

    labels = list(categories.keys())
    sizes = list(categories.values())
    colors = ["#2ecc71", "#3498db", "#e74c3c", "#95a5a6"]
    explode = (0.05, 0.05, 0.05, 0.05)

    fig, ax = plt.subplots(figsize=(8, 8))

    wedges, texts, autotexts = ax.pie(
        sizes, explode=explode, labels=labels, colors=colors,
        autopct="%1.1f%%", startangle=140,
        textprops={"fontsize": 11}
    )

    ax.set_title("PHB Metabolic Pathway Completeness\nin GTDB Genomes", fontsize=14)

    # 图例
    ax.legend(wedges, [f"{l}: {s:,} genomes" for l, s in zip(labels, sizes)],
              title="Pathway Type", loc="lower right", fontsize=9)

    plt.tight_layout()
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    if logger:
        logger.info(f"通路饼图: {output_path}")


def plot_phb_summary_dashboard(phb_counts: Path, phylum_summary: Path,
                                pathway_file: Path, output_path: Path,
                                logger: logging.Logger = None):
    """综合仪表盘图 (4 个子图)。"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # --- 子图 1: 门分布 ---
    if phylum_summary.exists():
        df_phylum = pd.read_csv(phylum_summary, sep="\t").head(10)
        df_phylum = df_phylum.sort_values("n_genomes", ascending=True)
        colors1 = plt.cm.plasma(np.linspace(0.1, 0.9, len(df_phylum)))
        axes[0, 0].barh(df_phylum["phylum"], df_phylum["n_genomes"], color=colors1)
        axes[0, 0].set_title("A. PHB Genes by Phylum (Top 10)")
        axes[0, 0].set_xlabel("Number of Genomes")

    # --- 子图 2: 基因类别分布 ---
    if phb_counts.exists():
        df = pd.read_csv(phb_counts, sep="\t")
        gene_cols = [c for c in df.columns if c.endswith("_count")]
        if gene_cols:
            totals = {col: int(df[col].sum()) for col in gene_cols}
            names = list(totals.keys())
            values = list(totals.values())
            colors2 = plt.cm.viridis(np.linspace(0.1, 0.9, len(names)))
            axes[0, 1].bar(names, values, color=colors2)
            axes[0, 1].set_title("B. Total PHB Gene Counts by Category")
            axes[0, 1].set_ylabel("Count")
            axes[0, 1].tick_params(axis="x", rotation=45)

    # --- 子图 3: 基因数分布直方图 ---
    if phb_counts.exists():
        df = pd.read_csv(phb_counts, sep="\t")
        gene_cols = [c for c in df.columns if c.endswith("_count")]
        if gene_cols:
            df["total"] = df[gene_cols].sum(axis=1)
            axes[1, 0].hist(df["total"], bins=50, color="steelblue",
                           edgecolor="white", alpha=0.8)
            axes[1, 0].set_title("C. Distribution of PHB Gene Count per Genome")
            axes[1, 0].set_xlabel("Number of PHB Genes")
            axes[1, 0].set_ylabel("Number of Genomes")
            axes[1, 0].axvline(df["total"].median(), color="red", linestyle="--",
                              label=f"Median: {df['total'].median():.0f}")
            axes[1, 0].legend()

    # --- 子图 4: 通路完整性 ---
    if pathway_file.exists():
        pw = pd.read_csv(pathway_file, sep="\t")
        categories_count = pw["pathway_completeness"].value_counts().sort_index()
        axes[1, 1].bar(categories_count.index, categories_count.values,
                      color=plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(categories_count))))
        axes[1, 1].set_title("D. PHB Pathway Completeness")
        axes[1, 1].set_xlabel("Number of PHB Genes Present (0-4)")
        axes[1, 1].set_ylabel("Number of Genomes")

    plt.suptitle("PHB Gene Analysis Summary — GTDB Representative Genomes",
                 fontsize=16, fontweight="bold")
    plt.tight_layout()
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    if logger:
        logger.info(f"综合仪表盘: {output_path}")


def generate_report(phb_counts: Path, phylum_summary: Path,
                    output_path: Path, logger: logging.Logger = None):
    """生成文本分析报告。"""
    lines = []
    lines.append("=" * 70)
    lines.append("PHB Gene Analysis Report — GTDB Representative Genomes (R226)")
    lines.append("=" * 70)
    lines.append("")

    # 基本统计
    if phb_counts.exists():
        df = pd.read_csv(phb_counts, sep="\t")
        gene_cols = [c for c in df.columns if c.endswith("_count")]
        lines.append(f"Total genomes analyzed: {len(df):,}")
        lines.append(f"Genomes with PHB genes: {len(df):,}")
        lines.append("")
        lines.append("Gene counts by category:")
        for col in gene_cols:
            total = int(df[col].sum())
            genomes = int((df[col] > 0).sum())
            lines.append(f"  {col}: {total:,} genes in {genomes:,} genomes")
        lines.append("")

    # 门分布
    if phylum_summary.exists():
        ps = pd.read_csv(phylum_summary, sep="\t")
        lines.append("Top phyla by PHB-containing genomes:")
        for _, row in ps.head(10).iterrows():
            lines.append(f"  {row['phylum']}: {int(row['n_genomes']):,} genomes")
        lines.append("")

    lines.append("=" * 70)

    report_text = "\n".join(lines)

    with open(output_path, "w") as f:
        f.write(report_text)

    if logger:
        logger.info(f"分析报告: {output_path}")
        for line in lines:
            logger.info(line)

    return report_text


def main():
    parser = argparse.ArgumentParser(description="PHB 分析可视化")
    parser.add_argument("--all", action="store_true", default=True,
                        help="生成所有图表")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "06_visualization")

    # 输入文件
    phb_counts = PROCESSED_DIR / "phb_search_results.tsv"
    phylum_summary = TABLES_DIR / "phb_phylum_summary.tsv"
    pathway_file = TABLES_DIR / "phb_pathway_completeness.tsv"

    # 检查数据是否存在
    if not phb_counts.exists():
        logger.error("PHB 搜索结果不存在，请先运行 01_phb_search.py")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("生成可视化图表")
    logger.info("=" * 60)

    # 加载分类学
    taxonomy = load_gtdb_taxonomy()

    # --- 1. 门分布图 ---
    logger.info("--- 1. 门水平分布 ---")
    if phylum_summary.exists():
        plot_phylum_distribution(
            phylum_summary,
            FIGURES_DIR / "phylum_distribution.pdf",
            logger=logger
        )
    else:
        logger.warning(f"门汇总文件不存在: {phylum_summary}")

    # --- 2. 基因热图 ---
    logger.info("--- 2. PHB 基因热图 ---")
    plot_phb_gene_counts_heatmap(
        phb_counts, taxonomy,
        FIGURES_DIR / "phb_gene_heatmap.pdf",
        logger=logger
    )

    # --- 3. 通路完整性 ---
    logger.info("--- 3. 通路完整性饼图 ---")
    plot_pathway_completeness(
        pathway_file,
        FIGURES_DIR / "pathway_completeness.pdf",
        logger=logger
    )

    # --- 4. 综合仪表盘 ---
    logger.info("--- 4. 综合仪表盘 ---")
    plot_phb_summary_dashboard(
        phb_counts, phylum_summary, pathway_file,
        FIGURES_DIR / "phb_summary_dashboard.pdf",
        logger=logger
    )

    # --- 5. 分析报告 ---
    logger.info("--- 5. 文本报告 ---")
    generate_report(
        phb_counts, phylum_summary,
        TABLES_DIR / "analysis_report.txt",
        logger=logger
    )

    logger.info("=" * 60)
    logger.info("Step 6 完成! 所有图表保存在 results/figures/")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
