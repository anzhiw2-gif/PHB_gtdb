#!/usr/bin/env python3
"""
PHB_gtdb — Step 1b: 古菌 PhaZ 降解基因搜索
验证古菌中是否存在 phaZ 基因 (预期: 无命中)

与 01_phb_search.py 相同逻辑，但只搜索古菌基因组。

用法:
    conda run -n phb_gtdb python 01b_archaea_search.py --threads 30
    conda run -n phb_gtdb python 01b_archaea_search.py --threads 30 --limit 100  # 测试
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
    "phaZ_Pst":     "BAA32541.1",     # Stutzerimonas stutzeri — PMID: PMC91002
    "phaZ_Cte":     "AAA87070.1",     # Comamonas testosteroni — U16275.1
    "phaZ_Asp":     "BAA35137.1",     # Acidovorax sp. TP4 — AB015309.1
    "phaZ_Dac":     "BAA19791.1",     # Delftia acidovorans — AB003186.1
    "phaZ_Sex":     "AAB02914.1",     # Streptomyces exfoliatus — U58990.1
    "phaZ_Ple3":    "P52090.1",       # Paucimonas lemoignei — Swiss-Prot
    "phaZ_Ple1":    "WP_243656647.1", # P. lemoignei — Type I
    "phaZ_Ple2":    "WP_207907290.1", # P. lemoignei — Type I
    # === 胞内型 PhaZ (无 lipase box) ===
    "phaZ_Cne1":    "BAA33394.1",     # C. necator PhaZ1 — PMID: PMC94854
    "phaZ_Cne2":    "CAJ93939.1",     # C. necator PhaZ2
    "phaZ_Cne5":    "CAJ95805.1",     # C. necator PhaZ5
    # === 其他 ===
    "phaZ_Rpi":     "WKZ88401.1",     # Ralstonia pickettii
    "phaZ_Rpi2":    "UCA14981.1",     # Ralstonia pickettii
    "phaZ_Bsp":     "WP_128854079.1", # Bacillus sp. CDB3 — PMID: PMC1636284
}


def download_phaz_references(output_dir: Path, logger: logging.Logger):
    """下载 PhaZ 参考序列并构建 DIAMOND DB。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "phaz_references.fasta"
    diamond_db = output_dir / "phaz_db.dmnd"

    if diamond_db.exists():
        logger.info("DIAMOND 数据库已存在，跳过下载")
        return str(diamond_db)

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

    # DIAMOND DB
    run_cmd(f"diamond makedb --in {fasta_path} -d {diamond_db} --quiet",
            "diamond makedb", logger=logger)
    return str(diamond_db)


def get_archaeal_genomes(logger: logging.Logger) -> List[Path]:
    """获取古菌基因组列表。"""
    logger.info("加载 GTDB 分类学...")
    taxonomy = load_gtdb_taxonomy()

    domain_map = {}
    for _, row in tqdm(taxonomy.iterrows(), total=len(taxonomy), desc="解析分类"):
        gid = str(row.get("accession", "")).replace("RS_", "").replace("GB_", "")
        tax_str = str(row.get("taxonomy", ""))
        if gid and tax_str:
            tax = parse_gtdb_taxonomy(tax_str)
            domain_map[gid] = tax.get("domain", "Unknown")

    all_genomes = []
    for db in ["GCA", "GCF"]:
        db_dir = GTDB_GENOMES / db
        if db_dir.exists():
            all_genomes.extend(db_dir.rglob("*_genomic.fna.gz"))

    logger.info(f"基因组总数: {len(all_genomes)}")

    archaea = []
    n_bacteria = 0
    n_unknown = 0

    for g in all_genomes:
        gid = g.name.replace("_genomic.fna.gz", "")
        domain = domain_map.get(gid, "Unknown")

        if domain == "Archaea":
            archaea.append(g)
        elif domain == "Bacteria":
            n_bacteria += 1
        else:
            n_unknown += 1

    logger.info(f"过滤: 古菌 {len(archaea)}, 细菌(跳过) {n_bacteria}, 未知 {n_unknown}")
    return sorted(archaea)


def process_genome_pyrodigal(args: Tuple[Path, str, Path]) -> Dict:
    """Pyrodigal → DIAMOND blastp → 返回命中。"""
    genome_path, diamond_db, tmp_base = args
    genome_id = genome_path.name.replace("_genomic.fna.gz", "")

    try:
        tmp_dir = tmp_base / genome_id
        tmp_dir.mkdir(parents=True, exist_ok=True)

        protein_file = tmp_dir / f"{genome_id}.faa"

        # Pyrodigal 基因预测
        if not protein_file.exists() or protein_file.stat().st_size < 100:
            with gzip.open(genome_path, "rb") as fin:
                dna_seq = fin.read()

            import pyrodigal
            orf_finder = pyrodigal.GeneFinder(meta=True)
            genes = orf_finder.find_genes(dna_seq)

            with open(protein_file, "w") as f:
                for i, gene in enumerate(genes):
                    seq = gene.sequence()
                    if seq and len(seq) >= 30:
                        f.write(f">{genome_id}_{i+1}\n{seq}\n")

            del dna_seq

        if not protein_file.exists() or protein_file.stat().st_size < 100:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return {}

        # DIAMOND 搜索
        result_file = tmp_dir / f"{genome_id}.hits.tsv"
        run_cmd(
            f"diamond blastp -q {protein_file} -d {diamond_db} "
            f"-o {result_file} --outfmt 6 qseqid sseqid pident evalue "
            f"-e 1e-10 --id 30 --query-cover 50 "
            f"--threads 1 --quiet --max-target-seqs 5",
            f"DIAMOND:{genome_id}", timeout=300, check=False
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
                                best_evalue = min(best_evalue, evalue)
                                best_pident = max(best_pident, pident)
                        except (ValueError, IndexError):
                            continue

        shutil.rmtree(tmp_dir, ignore_errors=True)

        if n_hits > 0:
            return {
                "genome_id": genome_id,
                "domain": "Archaea",
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
    parser = argparse.ArgumentParser(description="古菌 PhaZ 降解基因搜索 (Pyrodigal+DIAMOND)")
    parser.add_argument("--threads", type=int, default=MAX_THREADS)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    logger = setup_logging(LOGS_DIR, "01b_archaea_search")

    logger.info("=" * 60)
    logger.info("古菌 PhaZ 降解基因搜索 (Pyrodigal + DIAMOND)")
    logger.info("验证古菌中是否存在 phaZ 基因")
    logger.info("=" * 60)

    # Step 1: 参考序列 DB
    diamond_db = download_phaz_references(EXTERNAL_DIR / "phb_references", logger)

    # Step 2: 古菌基因组列表
    logger.info("=" * 60)
    logger.info("获取古菌基因组")
    logger.info("=" * 60)
    archaea_genomes = get_archaeal_genomes(logger)
    if args.limit > 0:
        archaea_genomes = archaea_genomes[:args.limit]

    logger.info(f"古菌基因组总数: {len(archaea_genomes)}")

    if len(archaea_genomes) == 0:
        logger.info("未找到古菌基因组，检查 GTDB 分类数据")
        return

    # Step 3: 断点续传
    result_file = PROCESSED_DIR / "archaea_phb_search_results.tsv"
    completed = set()
    if result_file.exists() and result_file.stat().st_size > 0:
        try:
            existing = pd.read_csv(result_file, sep="\t")
            if "genome_id" in existing.columns:
                completed = set(existing["genome_id"])
        except Exception:
            pass

    to_process = [g for g in archaea_genomes
                  if g.name.replace("_genomic.fna.gz", "") not in completed]
    logger.info(f"已完成: {len(completed)}, 待处理: {len(to_process)}")

    if not to_process:
        logger.info("全部完成！")
        if result_file.exists():
            df = pd.read_csv(result_file, sep="\t")
            n_genomes = len(df)
            n_phaz = int(df["phaZ_count"].sum())
            if n_phaz == 0:
                logger.info("✓ 结论: 古菌中未发现 phaZ 基因 — 与预期一致")
            else:
                logger.warning(f"⚠ 古菌中发现 {n_phaz} 个 phaZ 基因 — 需进一步验证")
        return

    # Step 4: 并行搜索
    logger.info("=" * 60)
    logger.info(f"并行搜索 ({args.threads} 线程, chunksize=1)")
    logger.info("=" * 60)

    tmp_base = PROCESSED_DIR / "tmp_archaea"
    tmp_base.mkdir(parents=True, exist_ok=True)
    batch = []
    batch_size = 500

    with Pool(args.threads) as pool:
        tasks = [(g, diamond_db, tmp_base) for g in to_process]
        it = pool.imap_unordered(process_genome_pyrodigal, tasks, chunksize=1)
        for result in tqdm(it, total=len(to_process),
                           desc="Archaea PhaZ search", ncols=100):
            if result:
                batch.append(result)

            if len(batch) >= batch_size:
                new_df = pd.DataFrame(batch)
                new_df.to_csv(result_file, sep="\t", mode="a",
                             header=not result_file.exists(), index=False)
                logger.info(f"保存 {len(batch)} 条结果 (累计已完成)")
                batch = []

    # 保存剩余
    if batch:
        pd.DataFrame(batch).to_csv(result_file, sep="\t", mode="a",
                                   header=not result_file.exists(), index=False)

    # 汇总
    if result_file.exists():
        df = pd.read_csv(result_file, sep="\t").drop_duplicates(subset="genome_id")
        df.to_csv(result_file, sep="\t", index=False)
        n_genomes = len(df)
        n_phaz = int(df["phaZ_count"].sum())

        logger.info("=" * 60)
        logger.info(f"最终结果: 古菌中 {n_genomes} 个基因组含 PhaZ 同源序列, 共 {n_phaz} 个")
        if n_phaz == 0:
            logger.info("✓ 结论: 所有古菌基因组中未发现 phaZ 基因")
            logger.info("  与前期研究结论一致 — 古菌不携带 PHB depolymerase")
        else:
            logger.warning("⚠ 发现潜在 phaZ 同源序列 — 建议人工核查")
            for _, row in df.head(20).iterrows():
                logger.info(f"  {row['genome_id']}: {int(row['phaZ_count'])} hit(s), "
                          f"best_pident={row['best_pident']:.1f}%")

    logger.info("=" * 60)
    logger.info("古菌搜索完成！")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
