#!/usr/bin/env python3
"""
PHB_gtdb — Step 2: 提取 PHB 基因序列

从 Step 1 的搜索结果中提取 PHB 基因的蛋白序列，用于后续比对和系统发育分析。

功能:
1. 从 HMM 命中中提取完整蛋白序列
2. 按基因类别 (phaZ/phaC/...) 分组
3. 去除冗余序列 (CD-HIT)
4. 添加分类学注释到序列 ID

用法:
    python 02_extract_sequences.py
"""

import os
import sys
import gzip
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from typing import List, Dict, Set

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import (setup_logging, run_cmd, read_fasta, write_fasta,
                    sanitize_id, filter_by_length, load_gtdb_taxonomy,
                    parse_gtdb_taxonomy)


def extract_hit_sequences(
    genome_id: str,
    gene_names: List[str],
    min_len: int = 100,
) -> Dict[str, List[SeqRecord]]:
    """从基因组蛋白文件中提取指定基因的序列。

    通过在基因组上重新运行 Prodigal 和 hmmscan，提取命中蛋白。

    Returns:
        {gene_name: [SeqRecord, ...]}
    """
    results = defaultdict(list)

    # 定位基因组文件
    genome_path = None
    for db in ["GCA", "GCF"]:
        # 从 genome_id 提取路径前缀
        parts = genome_id.split("_")
        if len(parts) >= 3:
            # e.g., GCF_000006945.2 → GCF/000/006/945/
            prefix = parts[1]  # 000006945.2
            subdir = f"{prefix[0:3]}/{prefix[3:6]}/{prefix[6:9]}"
            candidate = GTDB_GENOMES / db / subdir / f"{genome_id}_genomic.fna.gz"
            if candidate.exists():
                genome_path = candidate
                break

    if genome_path is None:
        return results

    try:
        tmp_dir = PROCESSED_DIR / "tmp" / genome_id
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Prodigal (需要可seek文件，不能使用管道)
        protein_file = tmp_dir / f"{genome_id}_proteins.faa"
        dna_tmp = tmp_dir / f"{genome_id}_genomic.fna"
        if not protein_file.exists() or protein_file.stat().st_size == 0:
            import gzip as gz
            with gz.open(genome_path, "rb") as fin, open(dna_tmp, "wb") as fout:
                fout.write(fin.read())
            run_cmd(
                f"prodigal -i {dna_tmp} -a {protein_file} -o /dev/null -q",
                f"Prodigal: {genome_id}",
                timeout=300, check=False
            )
            if dna_tmp.exists():
                dna_tmp.unlink()

        if not protein_file.exists():
            return results

        # hmmscan
        hmm_db = EXTERNAL_DIR / "phb_hmm" / "phb_all.hmm"
        hmm_result = tmp_dir / f"{genome_id}_hmmscan.txt"

        run_cmd(
            f"hmmscan --cpu 1 --noali -E 1e-5 --domE 1e-3 "
            f"--tblout {hmm_result} {hmm_db} {protein_file} > /dev/null 2>&1",
            f"hmmscan: {genome_id}",
            timeout=600, check=False
        )

        if not hmm_result.exists():
            return results

        # 解析命中蛋白 ID
        hit_protein_ids = set()
        pfam_to_gene = {}
        for gene_id, gene_info in PHB_GENES.items():
            if "pfam" in gene_info:
                pfam_to_gene[gene_info["pfam"]] = gene_id

        with open(hmm_result) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    target = parts[0]
                    protein_id = parts[2]
                    if target in pfam_to_gene:
                        hit_protein_ids.add(protein_id)

        # 提取匹配的蛋白序列
        proteins = {r.id: r for r in SeqIO.parse(protein_file, "fasta")}
        for pid in hit_protein_ids:
            if pid in proteins:
                record = proteins[pid]
                if len(record.seq) >= min_len:
                    # 美化序列 ID
                    record.id = f"{genome_id}|{record.id}"
                    record.description = ""
                    results["all_hits"].append(record)

    except Exception as e:
        pass
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return results


def run_cdhit(
    input_fasta: Path,
    output_fasta: Path,
    identity: float = 0.95,
    logger: logging.Logger = None,
) -> Path:
    """使用 CD-HIT 去除冗余序列。

    Args:
        input_fasta: 输入 FASTA 文件
        output_fasta: 输出 FASTA 文件
        identity: 序列相似度阈值
    """
    tmp_file = output_fasta.with_suffix(".cdhit.tmp")

    run_cmd(
        f"cd-hit -i {input_fasta} -o {tmp_file} -c {identity} "
        f"-n 5 -M 16000 -d 0 -T 8 > /dev/null 2>&1",
        f"CD-HIT ({identity})",
        check=False, logger=logger
    )

    # CD-HIT 输出文件是 tmp_file (不带扩展名), 实际文件名有 .fasta 后缀
    # 取决于 cd-hit 版本
    cdhit_output = Path(str(tmp_file))
    if cdhit_output.exists():
        import shutil
        shutil.move(str(cdhit_output), str(output_fasta))
    elif Path(str(tmp_file) + ".fasta").exists():
        import shutil
        shutil.move(str(tmp_file) + ".fasta", str(output_fasta))

    return output_fasta


def add_taxonomy_to_fasta(
    input_fasta: Path,
    output_fasta: Path,
    taxonomy_df: pd.DataFrame,
) -> None:
    """将分类学信息添加到序列 ID 中，便于系统发育树着色。

    新 ID 格式: genomeID|phylum|originalID
    """
    # genome_id → phylum 映射
    id_to_phylum = {}
    for _, row in taxonomy_df.iterrows():
        genome_id = str(row.get("accession", "")).replace("RS_", "").replace("GB_", "")
        tax_str = str(row.get("taxonomy", ""))
        if tax_str:
            tax = parse_gtdb_taxonomy(tax_str)
            id_to_phylum[genome_id] = tax.get("phylum", "Unknown")

    records = read_fasta(input_fasta)
    for r in records:
        genome_id = r.id.split("|")[0]
        phylum = id_to_phylum.get(genome_id, "Unknown")
        phylum_clean = phylum.replace(" ", "_").replace("/", "_")
        r.id = f"{genome_id}|{phylum_clean}|{r.id.split('|')[-1] if '|' in r.id else r.id}"
        r.description = ""

    write_fasta(records, output_fasta)


def main():
    parser = argparse.ArgumentParser(description="提取 PHB 基因序列")
    parser.add_argument("--cdhit", type=float, default=0.95,
                        help="CD-HIT 相似度阈值 (默认: 0.95)")
    parser.add_argument("--no-cdhit", action="store_true",
                        help="跳过 CD-HIT 去冗余")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "02_extract_sequences")

    # 读取 Step 1 结果
    search_result = PROCESSED_DIR / "phb_search_results.tsv"
    if not search_result.exists():
        logger.error(f"Step 1 结果不存在: {search_result}")
        logger.info("请先运行 01_phb_search.py")
        sys.exit(1)

    df = pd.read_csv(search_result, sep="\t")
    logger.info(f"读取搜索结果: {len(df)} 个基因组含 PHB 基因")

    # 获取所有有 PHB 基因的基因组
    gene_cols = [c for c in df.columns if c != "genome_id"]
    genomes_with_phb = []
    for _, row in df.iterrows():
        for col in gene_cols:
            if row[col] > 0:
                genomes_with_phb.append(row["genome_id"])
                break

    logger.info(f"需要提取序列的基因组数: {len(genomes_with_phb)}")

    # 提取序列
    logger.info("提取 PHB 蛋白序列...")
    all_records = []
    for genome_id in tqdm(genomes_with_phb, desc="Extracting"):
        hits = extract_hit_sequences(genome_id, gene_cols)
        if "all_hits" in hits:
            all_records.extend(hits["all_hits"])

    logger.info(f"提取了 {len(all_records)} 条蛋白序列")

    # 保存原始序列
    raw_fasta = PROCESSED_DIR / "phb_proteins_raw.fasta"
    write_fasta(all_records, raw_fasta)
    logger.info(f"原始序列: {raw_fasta}")

    # CD-HIT 去冗余
    if not args.no_cdhit:
        logger.info(f"CD-HIT 去冗余 (identity={args.cdhit})...")
        dedup_fasta = PROCESSED_DIR / "phb_proteins_dedup.fasta"
        run_cdhit(raw_fasta, dedup_fasta, args.cdhit, logger)
        final_fasta = dedup_fasta
    else:
        final_fasta = raw_fasta

    n_seqs = sum(1 for _ in open(final_fasta))
    logger.info(f"去冗余后序列数: {n_seqs}")

    # 添加分类学注释
    logger.info("添加 GTDB 分类学注释...")
    taxonomy = load_gtdb_taxonomy()
    annotated_fasta = PROCESSED_DIR / "phb_proteins_annotated.fasta"
    add_taxonomy_to_fasta(final_fasta, annotated_fasta, taxonomy)
    logger.info(f"注释后序列: {annotated_fasta}")

    # 按基因类别分组
    logger.info("按基因类别分组...")
    gene_groups = defaultdict(list)
    for record in read_fasta(annotated_fasta):
        for gene_name in PHB_GENES:
            if gene_name in record.id.lower():
                gene_groups[gene_name].append(record)
                break
        else:
            gene_groups["other"].append(record)

    for gene_name, records in gene_groups.items():
        if records:
            out_path = PROCESSED_DIR / f"phb_{gene_name}.fasta"
            write_fasta(records, out_path)
            logger.info(f"  {gene_name}: {len(records)} 条序列 → {out_path.name}")

    logger.info("=" * 60)
    logger.info("Step 2 完成!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
