#!/usr/bin/env python3
"""
PHB_gtdb — Step 3: 多序列比对 (MSA)

对 PHB 蛋白序列进行多序列比对，为系统发育分析做准备。
支持 MAFFT (默认) 和 Muscle。

用法:
    python 03_msa.py [--gene phaZ] [--tool mafft]
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import List, Optional

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import setup_logging, run_cmd, read_fasta, write_fasta, count_sequences


def run_mafft(input_fasta: Path, output_fasta: Path, threads: int = 30,
              logger: logging.Logger = None) -> Path:
    """MAFFT 多序列比对。

    自动选择最优算法 (--auto): L-INS-i (>200条序列用FFT-NS-2)
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd = f"mafft --auto --thread {threads} --reorder {input_fasta} > {output_fasta}"

    if logger:
        n_seqs = count_sequences(input_fasta)
        logger.info(f"MAFFT 比对中 ({n_seqs} 条序列)...")

    run_cmd(cmd, f"MAFFT: {input_fasta.name}", timeout=86400, logger=logger)

    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        raise RuntimeError(f"MAFFT 输出为空: {output_fasta}")

    return output_fasta


def run_trimal(input_fasta: Path, output_fasta: Path,
               method: str = "automated1",
               logger: logging.Logger = None) -> Path:
    """trimAl 修剪比对中的低质量区域。

    方法选项:
    - automated1: 自动选择最优方法 (推荐)
    - gappyout: 移除 gap 比例高的列
    - strict: 严格修剪
    - nogaps: 移除所有含 gap 的列
    """
    cmd = f"trimal -in {input_fasta} -out {output_fasta} -{method}"

    # trimAl 可能未安装，尝试 conda 路径
    if logger:
        logger.info(f"trimAl ({method}) 修剪中...")

    result = run_cmd(cmd, f"trimAl: {input_fasta.name}",
                     check=False, logger=logger)

    # 如果 trimAl 失败，使用备选方案：删除全 gap 列
    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        logger.warning("trimAl 不可用，使用 Python 进行简单修剪...")
        alignment = AlignIO.read(input_fasta, "fasta")
        # 删除 >90% gap 的列
        aln_len = alignment.get_alignment_length()
        keep_cols = []
        for i in range(aln_len):
            col = alignment[:, i]
            gap_ratio = col.count("-") / len(col)
            if gap_ratio < 0.9:
                keep_cols.append(i)

        trimmed = alignment[:, keep_cols[0]:keep_cols[0]+1]
        for col in keep_cols[1:]:
            trimmed += alignment[:, col:col+1]

        AlignIO.write(trimmed, output_fasta, "fasta")
        logger.info(f"Python 修剪: {aln_len} → {len(keep_cols)} 列")

    return output_fasta


def calculate_conservation(alignment_file: Path) -> dict:
    """计算比对保守性统计。"""
    alignment = AlignIO.read(alignment_file, "fasta")
    n_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()

    # 计算 gap 比例
    total_gaps = sum(str(r.seq).count("-") for r in alignment)
    gap_pct = total_gaps / (n_seqs * aln_len) * 100

    # 计算完全保守列 (所有序列相同的列)
    conserved = 0
    for i in range(aln_len):
        col_chars = set(alignment[:, i].replace("-", ""))
        if len(col_chars) <= 1:
            conserved += 1
    conserved_pct = conserved / aln_len * 100

    return {
        "n_sequences": n_seqs,
        "alignment_length": aln_len,
        "gap_percentage": round(gap_pct, 2),
        "fully_conserved_cols": conserved,
        "conserved_percentage": round(conserved_pct, 2),
    }


def main():
    parser = argparse.ArgumentParser(description="多序列比对")
    parser.add_argument("--gene", type=str, default="all",
                        help="指定基因 (phaZ/phaC/all), 默认 all")
    parser.add_argument("--input", type=str, default=None,
                        help="自定义输入 FASTA 文件")
    parser.add_argument("--tool", type=str, default=ALIGNMENT_TOOL,
                        help=f"比对工具 (默认: {ALIGNMENT_TOOL})")
    parser.add_argument("--threads", type=int, default=ALIGNMENT_THREADS,
                        help=f"线程数 (默认: {ALIGNMENT_THREADS})")
    parser.add_argument("--no-trim", action="store_true",
                        help="跳过 trimAl 修剪")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "03_msa")

    # 确定输入文件
    if args.input:
        input_files = {"input": Path(args.input)}
    elif args.gene == "all":
        input_files = {}
        for gene_name in PHB_GENES:
            f = PROCESSED_DIR / f"phb_{gene_name}.fasta"
            if f.exists():
                input_files[gene_name] = f
        # 也处理合并文件
        combined = PROCESSED_DIR / "phb_proteins_annotated.fasta"
        if combined.exists():
            input_files["all_phb"] = combined
    else:
        f = PROCESSED_DIR / f"phb_{args.gene}.fasta"
        if not f.exists():
            logger.error(f"文件不存在: {f}")
            logger.info("请先运行 02_extract_sequences.py")
            sys.exit(1)
        input_files = {args.gene: f}

    if not input_files:
        logger.error("未找到任何输入文件")
        sys.exit(1)

    logger.info(f"待比对: {list(input_files.keys())}")

    # 对每个基因进行比对
    for gene_name, input_fasta in input_files.items():
        logger.info("=" * 60)
        logger.info(f"比对基因: {gene_name}")
        logger.info("=" * 60)

        n_seqs = count_sequences(input_fasta)
        logger.info(f"序列数: {n_seqs}")

        if n_seqs < 3:
            logger.warning(f"{gene_name} 序列数不足 (<3)，跳过比对")
            continue

        # Step 1: MAFFT 比对
        aligned_fasta = PROCESSED_DIR / f"{gene_name}_aligned.fasta"
        run_mafft(input_fasta, aligned_fasta, args.threads, logger)

        # 比对质量评估
        stats = calculate_conservation(aligned_fasta)
        logger.info(f"比对统计: {stats}")

        # Step 2: trimAl 修剪
        if not args.no_trim:
            trimmed_fasta = PROCESSED_DIR / f"{gene_name}_trimmed.fasta"
            run_trimal(aligned_fasta, trimmed_fasta, logger=logger)

            trim_stats = calculate_conservation(trimmed_fasta)
            logger.info(f"修剪后统计: {trim_stats}")
        else:
            logger.info("跳过 trimAl 修剪")

    logger.info("=" * 60)
    logger.info("Step 3 完成!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
