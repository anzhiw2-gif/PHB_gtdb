#!/usr/bin/env python3
"""
PHB_gtdb — Step 1: PhaZ 降解基因搜索 (优化版)

优化点:
- Pyrodigal (进程内基因预测, 无 subprocess 开销)
- DIAMOND (比 BLASTP 快 100-1000x)
- chunksize=1 负载均衡
- 增量保存 + 断点续传

用法:
    conda run -n phb_gtdb python 01_phb_search.py --threads 30
    conda run -n phb_gtdb python 01_phb_search.py --threads 30 --limit 100  # 测试
"""

import os, sys, gzip, time, shutil, argparse, logging
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from typing import List, Dict, Tuple

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import *
from utils import setup_logging, run_cmd, load_gtdb_taxonomy, parse_gtdb_taxonomy

# ==============================================================================
# PhaZ 参考序列 — 已验证的 PHB/PHA depolymerase (NCBI 实时验证, 2026-06-05)
# 详见 docs/PHAZ_REFERENCES.md
# ==============================================================================
PHAZ_REFERENCES = {
    # === 胞外型 PhaZ (含 lipase box + 信号肽) ===
    # Stutzerimonas stutzeri (Pseudomonas stutzeri) — 576aa, S163-D240-H289
    # 文献: Ohura et al., 1999, Appl Environ Microbiol. PMID: PMC91002
    "phaZ_Pst":     "BAA32541.1",
    # Comamonas testosteroni 31A — PHA depolymerase precursor
    # 文献: Jendrossek et al., 1995. GenBank: U16275.1
    "phaZ_Cte":     "AAA87070.1",
    # Acidovorax sp. TP4 — PHB depolymerase, Type II 催化域
    # GenBank: AB015309.1
    "phaZ_Asp":     "BAA35137.1",
    # Delftia acidovorans (Comamonas acidovorans) — PHB depolymerase
    # GenBank: AB003186.1
    "phaZ_Dac":     "BAA19791.1",
    # Streptomyces exfoliatus — 胞外型, 革兰氏阳性菌
    # 文献: García-Hidalgo et al., 2012. GenBank: U58990.1
    "phaZ_Sex":     "AAB02914.1",
    # Paucimonas lemoignei — PHA depolymerase C (PhaZ3), Swiss-Prot 审阅
    # UniProt: P52090
    "phaZ_Ple3":    "P52090.1",
    # Paucimonas lemoignei — 胞外型 Type I 催化域
    "phaZ_Ple1":    "WP_243656647.1",
    "phaZ_Ple2":    "WP_207907290.1",

    # === 胞内型 PhaZ (无 lipase box) ===
    # Cupriavidus necator H16 PhaZ1 — 首个克隆的胞内 PhaZ, 419aa
    # 文献: Saegusa et al., 2001, J Bacteriol. PMID: PMC94854
    "phaZ_Cne1":    "BAA33394.1",
    # C. necator H16 PhaZ2 — 胞内型
    # 文献: Pohlmann et al., 2006, Nat Biotechnol.
    "phaZ_Cne2":    "CAJ93939.1",
    # C. necator H16 PhaZ5/PhaZd — 胞内型
    "phaZ_Cne5":    "CAJ95805.1",

    # === 其他 ===
    # Ralstonia pickettii — PHA depolymerase
    "phaZ_Rpi":     "WKZ88401.1",
    "phaZ_Rpi2":    "UCA14981.1",
    # Bacillus sp. CDB3 — 胞内型, 含 lipase box (新型 PhaZ), 革兰氏阳性菌
    # 文献: Tseng et al., 2006, J Bacteriol. PMID: PMC1636284
    "phaZ_Bsp":     "WP_128854079.1",
}

# GTDB 细菌门
BACTERIAL_PHYLA = {
    "Pseudomonadota", "Actinomycetota", "Bacillota", "Bacteroidota",
    "Myxococcota", "Cyanobacteriota", "Chloroflexota", "Bdellovibrionota",
    "Desulfobacterota", "Planctomycetota", "Verrucomicrobiota", "Acidobacteriota",
    "Spirochaetota", "Thermodesulfobacteriota", "Nitrospirota", "Gemmatimonadota",
    "Deinococcota", "Armatimonadota", "Synergistota", "Fusobacteriota",
    "Campylobacterota", "Aquificota", "Thermotogota", "Deferribacterota",
    "Dictyoglomerota", "Elusimicrobiota", "Fibrobacterota",
}


def download_phaz_references(output_dir: Path, logger: logging.Logger) -> tuple:
    """下载 PhaZ 参考序列并构建 BLAST + DIAMOND DB。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "phaz_references.fasta"
    blast_db = output_dir / "phaz_db"
    diamond_db = output_dir / "phaz_db.dmnd"

    if diamond_db.with_suffix(".dmnd").exists():
        logger.info("BLAST + DIAMOND 数据库已存在")
        return blast_db, str(diamond_db)

    logger.info(f"下载 {len(PHAZ_REFERENCES)} 条 PhaZ 参考序列...")
    if not fasta_path.exists() or fasta_path.stat().st_size < 500:
        for name, acc in tqdm(PHAZ_REFERENCES.items(), desc="NCBI download"):
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                f"db=protein&id={acc}&rettype=fasta&retmode=text"
            )
            run_cmd(f"wget -q -O - '{url}' >> {fasta_path}",
                    f"download {acc}", check=False, logger=None)
            time.sleep(0.35)

    n = sum(1 for l in open(fasta_path) if l.startswith(">"))
    logger.info(f"下载完成: {n} 条 PhaZ 参考序列")

    # BLAST DB
    run_cmd(f"makeblastdb -in {fasta_path} -dbtype prot -out {blast_db}",
            "makeblastdb", logger=logger)

    # DIAMOND DB
    run_cmd(f"diamond makedb --in {fasta_path} -d {diamond_db} --quiet",
            "diamond makedb", logger=logger)

    return blast_db, str(diamond_db)


def get_bacterial_genomes(logger: logging.Logger) -> List[Path]:
    """获取细菌基因组列表，排除古菌。"""
    logger.info("加载 GTDB 分类学...")
    taxonomy = load_gtdb_taxonomy()

    phylum_map = {}
    domain_map = {}
    for _, row in tqdm(taxonomy.iterrows(), total=len(taxonomy), desc="解析分类"):
        gid = str(row.get("accession", "")).replace("RS_", "").replace("GB_", "")
        tax_str = str(row.get("taxonomy", ""))
        if gid and tax_str:
            tax = parse_gtdb_taxonomy(tax_str)
            phylum_map[gid] = tax.get("phylum", "Unknown")
            domain_map[gid] = tax.get("domain", "Unknown")

    all_genomes = []
    for db in ["GCA", "GCF"]:
        db_dir = GTDB_GENOMES / db
        if db_dir.exists():
            all_genomes.extend(db_dir.rglob("*_genomic.fna.gz"))

    logger.info(f"基因组总数: {len(all_genomes)}")

    bacterial = []
    n_archaea = 0
    n_unknown = 0
    n_other_phylum = 0

    for g in all_genomes:
        gid = g.name.replace("_genomic.fna.gz", "")
        domain = domain_map.get(gid, "Unknown")
        phylum = phylum_map.get(gid, "")

        if domain == "Archaea":
            n_archaea += 1; continue
        elif domain != "Bacteria":
            n_unknown += 1; continue
        elif phylum and phylum not in BACTERIAL_PHYLA:
            n_other_phylum += 1; continue
        else:
            bacterial.append(g)

    logger.info(f"过滤: 排除古菌{n_archaea}, 非细菌{n_unknown}, 非目标门{n_other_phylum}")
    logger.info(f"保留细菌基因组: {len(bacterial)}")
    return sorted(bacterial)


def process_genome_pyrodigal(args: Tuple[Path, str, Path, bool]) -> Dict:
    """Pyrodigal (进程内) → DIAMOND blastp → 返回命中。"""
    genome_path, diamond_db, tmp_base, use_diamond = args
    genome_id = genome_path.name.replace("_genomic.fna.gz", "")

    try:
        tmp_dir = tmp_base / genome_id
        tmp_dir.mkdir(parents=True, exist_ok=True)

        protein_file = tmp_dir / f"{genome_id}.faa"
        dna_tmp = tmp_dir / f"{genome_id}.fna"

        # Pyrodigal (进程内基因预测, 无 subprocess 开销)
        if not protein_file.exists() or protein_file.stat().st_size < 100:
            # 解压基因组
            with gzip.open(genome_path, "rb") as fin:
                dna_seq = fin.read()

            # Pyrodigal 基因预测
            import pyrodigal
            orf_finder = pyrodigal.GeneFinder(meta=True)
            genes = orf_finder.find_genes(dna_seq)

            # 写入蛋白文件 (迭代每个基因)
            with open(protein_file, "w") as f:
                for i, gene in enumerate(genes):
                    seq = gene.translate()
                    if seq and len(seq) >= 30:  # 最少 30aa
                        f.write(f">{genome_id}_{i+1}\n{seq}\n")

            del dna_seq  # 释放内存

        if not protein_file.exists() or protein_file.stat().st_size < 100:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return {}

        # 搜索 vs PhaZ DB
        result_file = tmp_dir / f"{genome_id}.hits.tsv"

        if use_diamond:
            run_cmd(
                f"diamond blastp -q {protein_file} -d {diamond_db} "
                f"-o {result_file} --outfmt 6 qseqid sseqid pident evalue "
                f"-e 1e-10 --id 30 --query-cover 50 "
                f"--threads 1 --quiet --max-target-seqs 5",
                f"DIAMOND:{genome_id}", timeout=300, check=False
            )
        else:
            blast_db = diamond_db  # fallback: use BLAST DB path
            run_cmd(
                f"blastp -query {protein_file} -db {blast_db} "
                f"-out {result_file} -outfmt '6 qseqid sseqid pident evalue' "
                f"-evalue 1e-10 -num_threads 1 -max_target_seqs 5",
                f"BLASTP:{genome_id}", timeout=300, check=False
            )

        # 解析结果
        n_hits = 0
        best_evalue = 999
        best_pident = 0
        hit_refs = set()

        if result_file.exists() and result_file.stat().st_size > 0:
            with open(result_file) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        try:
                            pident = float(parts[2])
                            evalue = float(parts[3])
                            if evalue <= 1e-10 and pident >= 30:
                                n_hits += 1
                                hit_refs.add(parts[1])
                                if evalue < best_evalue:
                                    best_evalue = evalue
                                if pident > best_pident:
                                    best_pident = pident
                        except (ValueError, IndexError):
                            continue

        shutil.rmtree(tmp_dir, ignore_errors=True)

        if n_hits > 0:
            return {
                "genome_id": genome_id,
                "phaZ_count": n_hits,
                "best_evalue": best_evalue,
                "best_pident": best_pident,
                "hit_refs": ",".join(sorted(hit_refs)),
            }
        return {}

    except Exception as e:
        try:
            shutil.rmtree(tmp_base / genome_id, ignore_errors=True)
        except:
            pass
        return {}


def main():
    parser = argparse.ArgumentParser(description="PhaZ 降解基因搜索 (Pyrodigal+DIAMOND)")
    parser.add_argument("--threads", type=int, default=MAX_THREADS)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--no-diamond", action="store_true",
                        help="使用 BLASTP 替代 DIAMOND")
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "01_phb_search")

    use_diamond = not args.no_diamond
    logger.info("=" * 60)
    logger.info("PhaZ 降解基因搜索 (Pyrodigal + %s)",
                 "DIAMOND" if use_diamond else "BLASTP")
    logger.info("=" * 60)

    # Step 1: 参考序列 DB
    blast_db, diamond_db = download_phaz_references(
        EXTERNAL_DIR / "phb_references", logger)
    search_db = diamond_db if use_diamond else str(blast_db)

    # Step 2: 基因组列表
    logger.info("=" * 60)
    logger.info("获取细菌基因组 (排除古菌)")
    logger.info("=" * 60)
    genomes = get_bacterial_genomes(logger)
    if args.limit > 0:
        genomes = genomes[:args.limit]

    # Step 3: 断点续传
    result_file = PROCESSED_DIR / "phb_search_results.tsv"
    completed = set()
    if result_file.exists() and result_file.stat().st_size > 0:
        try:
            existing = pd.read_csv(result_file, sep="\t")
            if "genome_id" in existing.columns:
                completed = set(existing["genome_id"])
        except Exception:
            pass

    to_process = [g for g in genomes
                  if g.name.replace("_genomic.fna.gz", "") not in completed]
    logger.info(f"已完成: {len(completed)}, 待处理: {len(to_process)}")

    if not to_process:
        logger.info("全部完成！")
        # 统计
        if result_file.exists():
            df = pd.read_csv(result_file, sep="\t")
            n_genomes = len(df)
            n_phaz = int(df["phaZ_count"].sum())
            logger.info(f"结果: {n_genomes} 个基因组含 PhaZ, 共 {n_phaz} 个 PhaZ 基因")
        return

    # Step 4: 并行搜索
    logger.info("=" * 60)
    logger.info(f"并行搜索 ({args.threads} 线程, chunksize=1)")
    logger.info("=" * 60)

    tmp_base = PROCESSED_DIR / "tmp"
    tmp_base.mkdir(parents=True, exist_ok=True)
    batch = []
    batch_size = 500

    with Pool(args.threads) as pool:
        tasks = [(g, search_db, tmp_base, use_diamond) for g in to_process]
        it = pool.imap_unordered(process_genome_pyrodigal, tasks, chunksize=1)
        processed = 0
        for result in tqdm(it, total=len(to_process),
                           desc="PhaZ search", ncols=100):
            processed += 1
            if result:
                batch.append(result)

            if len(batch) >= batch_size:
                new_df = pd.DataFrame(batch)
                new_df.to_csv(result_file, sep="\t", mode="a",
                             header=not result_file.exists(), index=False)
                logger.info(f"保存 {len(batch)} 条结果 "
                          f"(已处理 {processed}/{len(to_process)} 基因组)")
                batch = []

            # 每 5000 基因组输出进度
            if processed % 5000 == 0:
                logger.info(f"进度: {processed}/{len(to_process)} "
                          f"({100*processed/len(to_process):.1f}%), "
                          f"累计命中: {len(batch) + (int(pd.read_csv(result_file, sep='\t')['phaZ_count'].sum()) if result_file.exists() else 0)} PhaZ")

    # 保存剩余
    if batch:
        pd.DataFrame(batch).to_csv(result_file, sep="\t", mode="a",
                                   header=not result_file.exists(), index=False)

    # 即使零命中也创建结果文件
    if not result_file.exists():
        pd.DataFrame(columns=["genome_id", "phaZ_count", "best_evalue",
                             "best_pident", "hit_refs"]).to_csv(
            result_file, sep="\t", index=False)

    # 去重汇总
    if result_file.exists():
        df = pd.read_csv(result_file, sep="\t").drop_duplicates(subset="genome_id")
        df.to_csv(result_file, sep="\t", index=False)
        n_genomes = len(df)
        n_phaz = int(df["phaZ_count"].sum())
        logger.info(f"最终结果: {n_genomes} 个基因组含 PhaZ, 共 {n_phaz} 个 PhaZ 基因")

        # Top 10
        if n_genomes > 0:
            logger.info("Top 10 PhaZ 命中:")
            for _, row in df.nlargest(10, "phaZ_count").iterrows():
                extras = ""
                if "best_pident" in row:
                    extras = f" (best: {row['best_pident']:.1f}%)"
                logger.info(f"  {row['genome_id']}: {int(row['phaZ_count'])} PhaZ{extras}")

    logger.info("=" * 60)
    logger.info("完成！")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
