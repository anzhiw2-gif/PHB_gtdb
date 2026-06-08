#!/usr/bin/env python3
"""
PHB_gtdb — Step 2: 提取 PhaZ 蛋白序列

从 Step 1 的 DIAMOND 搜索结果中提取 PhaZ 蛋白序列。

对每个含 PhaZ 的基因组:
1. 运行 Pyrodigal 预测所有蛋白
2. 运行 DIAMOND 找 PhaZ 匹配
3. 提取匹配的蛋白序列
4. 添加 GTDB 分类学注释

用法:
    python 02_extract_sequences.py --threads 30
    python 02_extract_sequences.py --threads 30 --limit 50  # 测试
"""

import os, sys, gzip, shutil, argparse, logging
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from typing import List, Dict, Tuple

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import (setup_logging, run_cmd, write_fasta, read_fasta,
                    load_gtdb_taxonomy, parse_gtdb_taxonomy)


def find_genome_path(genome_id: str) -> Path:
    """根据 GTDB accession 定位基因组文件。"""
    parts = genome_id.split("_")
    if len(parts) < 2:
        return None
    prefix = parts[1]  # e.g. 000952205.1
    subdir = f"{prefix[0:3]}/{prefix[3:6]}/{prefix[6:9]}"
    for db in ["GCA", "GCF"]:
        candidate = GTDB_GENOMES / db / subdir / f"{genome_id}_genomic.fna.gz"
        if candidate.exists():
            return candidate
    return None


def process_one_genome(args: Tuple[str, str, Path]) -> Dict:
    """单个基因组: Pyrodigal → DIAMOND → 提取匹配蛋白。"""
    genome_id, diamond_db, tmp_base = args
    result = {"genome_id": genome_id, "records": []}

    try:
        genome_path = find_genome_path(genome_id)
        if genome_path is None:
            return result

        tmp_dir = tmp_base / genome_id
        tmp_dir.mkdir(parents=True, exist_ok=True)

        protein_file = tmp_dir / f"{genome_id}.faa"

        # Pyrodigal 预测
        if not protein_file.exists() or protein_file.stat().st_size < 100:
            with gzip.open(genome_path, "rb") as f:
                dna_seq = f.read()

            gf = pyrodigal.GeneFinder(meta=True)
            genes = gf.find_genes(dna_seq)

            with open(protein_file, "w") as f:
                for i, gene in enumerate(genes):
                    seq = gene.translate()  # 蛋白翻译
                    if seq and len(seq) >= 30:
                        f.write(f">{genome_id}_{i+1}\n{seq}\n")
            del dna_seq

        if not protein_file.exists() or protein_file.stat().st_size < 100:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return result

        # DIAMOND
        result_file = tmp_dir / f"{genome_id}.hits.tsv"
        run_cmd(
            f"diamond blastp -q {protein_file} -d {diamond_db} "
            f"-o {result_file} --outfmt 6 qseqid sseqid pident evalue "
            f"-e 1e-10 --id 30 --query-cover 50 "
            f"--threads 1 --quiet --max-target-seqs 5 --ignore-warnings",
            f"DIAMOND:{genome_id}", timeout=300, check=False
        )

        # 提取匹配的蛋白
        if result_file.exists() and result_file.stat().st_size > 0:
            # 获取命中的蛋白ID
            hit_ids = set()
            best_each = {}
            with open(result_file) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        qid = parts[0]
                        pident = float(parts[2])
                        hit_ids.add(qid)
                        if qid not in best_each or pident > best_each[qid][0]:
                            best_each[qid] = (pident, parts[3], parts[1])

            # 提取序列
            proteins = {}
            for rec in SeqIO.parse(protein_file, "fasta"):
                proteins[rec.id] = str(rec.seq)

            for pid in hit_ids:
                if pid in proteins:
                    best = best_each.get(pid, (0, "N/A", "N/A"))
                    rec = SeqRecord(
                        id=f"{genome_id}|{pid}",
                        seq=proteins[pid],
                        description=f"pident={best[0]:.1f}% eval={best[1]} ref={best[2]}"
                    )
                    result["records"].append(rec)

        shutil.rmtree(tmp_dir, ignore_errors=True)

    except Exception as e:
        sys.stderr.write(f"[ERROR] {genome_id}: {e}\n")
        sys.stderr.flush()
        try:
            shutil.rmtree(tmp_base / genome_id, ignore_errors=True)
        except:
            pass

    return result


def main():
    parser = argparse.ArgumentParser(description="提取 PhaZ 蛋白序列")
    parser.add_argument("--threads", type=int, default=MAX_THREADS)
    parser.add_argument("--limit", type=int, default=0,
                        help="限制处理基因组数 (测试用)")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "02_extract_sequences")

    # 读取 Step 1 结果
    result_files = [
        PROCESSED_DIR / "phb_search_results.tsv",
        PROCESSED_DIR / "archaea_phb_search_results.tsv",
    ]

    all_genomes = []
    for rf in result_files:
        if rf.exists():
            df = pd.read_csv(rf, sep="\t")
            if "phaZ_count" in df.columns:
                hits = df[df["phaZ_count"] > 0]
                logger.info(f"{rf.name}: {len(hits)} 个基因组含 PhaZ")
                all_genomes.extend(hits["genome_id"].tolist())

    logger.info(f"总计: {len(all_genomes)} 个基因组需提取序列")

    if args.limit > 0:
        all_genomes = all_genomes[:args.limit]
        logger.info(f"限制模式: 仅处理前 {args.limit} 个")

    # DIAMOND DB
    diamond_db = str(EXTERNAL_DIR / "phb_references" / "phaz_db.dmnd")

    # 并行提取
    tmp_base = PROCESSED_DIR / "tmp_extract"
    tmp_base.mkdir(parents=True, exist_ok=True)

    all_records = []

    with Pool(args.threads) as pool:
        tasks = [(gid, diamond_db, tmp_base) for gid in all_genomes]
        it = pool.imap_unordered(process_one_genome, tasks, chunksize=1)
        for result in tqdm(it, total=len(all_genomes),
                           desc="Extracting PhaZ", ncols=100):
            if result and result["records"]:
                all_records.extend(result["records"])

    logger.info(f"提取了 {len(all_records)} 条 PhaZ 蛋白序列")

    # 加载分类学
    logger.info("添加 GTDB 分类学注释...")
    taxonomy_bac = load_gtdb_taxonomy()
    taxonomy_arc = load_gtdb_taxonomy(GTDB_TAXONOMY_AR53)
    taxonomy_all = pd.concat([taxonomy_bac, taxonomy_arc], ignore_index=True)

    gid_to_tax = {}
    for _, row in taxonomy_all.iterrows():
        gid = str(row["accession"]).replace("RS_", "").replace("GB_", "")
        tax = parse_gtdb_taxonomy(str(row.get("taxonomy", "")))
        gid_to_tax[gid] = tax

    # 添加分类学到序列ID: genomeID|phylum|genus|species|originalID
    annotated = []
    for rec in all_records:
        genome_id = rec.id.split("|")[0]
        tax = gid_to_tax.get(genome_id, {})
        phylum = tax.get("phylum", "Unknown").replace(" ", "_")
        genus = tax.get("genus", "Unknown").replace(" ", "_")
        species = tax.get("species", "sp.").replace(" ", "_")[:50]

        new_id = f"{genome_id}|{phylum}|{genus}|{species}|{rec.id.split('|')[-1]}"
        rec.id = new_id
        rec.description = f"[domain={tax.get('domain','?')}] " + rec.description
        annotated.append(rec)

    # 保存
    output_fasta = PROCESSED_DIR / "phaz_proteins_all.fasta"
    write_fasta(annotated, output_fasta)
    logger.info(f"最终序列文件: {output_fasta} ({len(annotated)} 条)")

    # 统计
    phyla = defaultdict(int)
    for rec in annotated:
        phylum = rec.id.split("|")[1] if "|" in rec.id else "Unknown"
        phyla[phylum] += 1

    logger.info("门级分布:")
    for phylum, count in sorted(phyla.items(), key=lambda x: -x[1]):
        logger.info(f"  {phylum}: {count}")

    logger.info("=" * 60)
    logger.info("Step 2 完成!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
