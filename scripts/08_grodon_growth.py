#!/usr/bin/env python3
"""Run gRodon2 growth predictions for matched phaZ+ and phaZ- genomes.

This script is designed for T141. It:
1. builds a same-genus matched manifest;
2. predicts CDS/proteins with Pyrodigal;
3. marks ribosomal-protein-like genes with HMMER or a gRodon seed DIAMOND database;
4. calls gRodon::predictGrowth for each genome;
5. writes growth predictions and a same-genus summary table.

The default HMMER route is intended for the formal run. The DIAMOND route is
kept for quick pilot checks against the gRodon example ribosomal seed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import pyrodigal


PROJECT_ROOT = Path(os.environ.get("PHB_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
GTDB_ROOT = Path(os.environ.get("GTDB_ROOT", "/home/data/haoyu/GTDB"))
GTDB_GENOMES = GTDB_ROOT / "gtdb_genomes_reps_r232" / "database"
GTDB_TAXONOMY = GTDB_ROOT / "taxonomy" / "bac120_taxonomy_r232.tsv"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
TABLES_DIR = RESULTS_DIR / "tables"
LOGS_DIR = RESULTS_DIR / "logs"
EXTERNAL_DIR = PROJECT_ROOT / "data" / "external" / "grodon"

PHB_PYTHON_ENV = Path(os.environ.get("PHB_PYTHON_ENV", "/home/data/haoyu/miniconda3/envs/phb_gtdb"))
GRODON_ENV = Path(os.environ.get("GRODON_ENV", "/home/data/haoyu/miniconda3/envs/grodon2"))
DIAMOND = PHB_PYTHON_ENV / "bin" / "diamond"
HMMSEARCH = PHB_PYTHON_ENV / "bin" / "hmmsearch"
RSCRIPT = GRODON_ENV / "bin" / "Rscript"
RUN_GRODON_ONE = PROJECT_ROOT / "scripts" / "08_run_grodon_one.R"

SUBTYPE_BY_REF = {
    "BAA33394.1": "intracellular.cupriavidus_like",
    "CAJ93939.1": "intracellular.cupriavidus_like",
    "CAJ95805.1": "intracellular.cupriavidus_like",
    "UCA14981.1": "intracellular.ralstonia_like",
    "WKZ88401.1": "intracellular.ralstonia_like",
    "BAA19791.1": "extracellular.general",
    "AAA87070.1": "extracellular.general",
    "BAA35137.1": "extracellular.general",
    "BAA32541.1": "extracellular.general",
    "AAB02914.1": "extracellular.general",
    "P52090.1": "extracellular.lemoignei_like",
    "WP_207907290.1": "extracellular.lemoignei_like",
    "WP_243656647.1": "extracellular.lemoignei_like",
    "WP_128854079.1": "intracellular.bacillus_like",
}


def parse_taxonomy_string(tax: str) -> dict[str, str]:
    out = {}
    ranks = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }
    for part in tax.split(";"):
        if "__" not in part:
            continue
        prefix, value = part.split("__", 1)
        key = ranks.get(prefix)
        if key:
            out[key] = value or "Unknown"
    return out


def normalize_gid(raw: str) -> str:
    return raw.replace("RS_", "").replace("GB_", "")


def genome_id_from_path(path: Path) -> str:
    return path.name.replace("_genomic.fna.gz", "")


def genome_path_from_id(genome_id: str) -> Path:
    match = re.match(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+$", genome_id)
    if not match:
        raise ValueError(f"Cannot construct GTDB genome path from accession: {genome_id}")
    db, part1, part2, part3 = match.groups()
    return GTDB_GENOMES / db / part1 / part2 / part3 / f"{genome_id}_genomic.fna.gz"


def build_genome_path_map(genome_ids) -> dict[str, str]:
    path_map = {}
    for gid in genome_ids:
        try:
            path_map[gid] = str(genome_path_from_id(gid))
        except ValueError:
            continue
    return path_map


def take_existing_genomes(candidates: list[str], path_map: dict[str, str], limit: int) -> list[str]:
    selected = []
    for gid in candidates:
        genome_path = path_map.get(gid)
        if genome_path and Path(genome_path).exists():
            selected.append(gid)
            if len(selected) >= limit:
                break
    return selected


def load_taxonomy() -> dict[str, dict[str, str]]:
    tax_by_gid = {}
    with open(GTDB_TAXONOMY, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gid = normalize_gid(parts[0])
            tax = parse_taxonomy_string(parts[1])
            if tax.get("domain") == "Bacteria":
                tax_by_gid[gid] = tax
    return tax_by_gid


def load_phaz_positive() -> tuple[dict[str, int], dict[str, Counter]]:
    fasta = PROCESSED_DIR / "phaz_proteins_validated.fasta"
    pos_counts: Counter[str] = Counter()
    subtype_counts: dict[str, Counter] = defaultdict(Counter)
    ref_re = re.compile(r"ref=([^\s]+)")
    with open(fasta, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            gid = header.split("|", 1)[0]
            ref_match = ref_re.search(header)
            subtype = "unknown"
            if ref_match:
                ref_raw = ref_match.group(1)
                for acc, label in SUBTYPE_BY_REF.items():
                    if acc in ref_raw:
                        subtype = label
                        break
            pos_counts[gid] += 1
            subtype_counts[gid][subtype] += 1
    return dict(pos_counts), subtype_counts


def load_initial_hit_genomes() -> set[str]:
    p = PROCESSED_DIR / "phb_search_results.tsv"
    if not p.exists():
        return set()
    hits = set()
    with open(p, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hits.add(row["genome_id"])
    return hits


def build_manifest(
    max_per_genus: int,
    pilot: int,
    include_initial_hits_as_negative: bool,
    seed: int = 42,
) -> list[dict[str, str]]:
    random.seed(seed)
    tax_by_gid = load_taxonomy()
    path_map = build_genome_path_map(tax_by_gid.keys())
    pos_counts, subtype_counts = load_phaz_positive()
    initial_hit_genomes = load_initial_hit_genomes()

    by_genus: dict[str, dict[str, list[str]]] = defaultdict(lambda: {"pos": [], "neg": []})
    for gid, tax in tax_by_gid.items():
        genus = tax.get("genus", "Unknown")
        if genus in ("", "Unknown"):
            continue
        if gid not in path_map:
            continue
        if gid in pos_counts:
            by_genus[genus]["pos"].append(gid)
        elif include_initial_hits_as_negative or gid not in initial_hit_genomes:
            by_genus[genus]["neg"].append(gid)

    rows: list[dict[str, str]] = []
    for genus, groups in sorted(by_genus.items()):
        pos_candidates = sorted(groups["pos"])
        neg_candidates = sorted(groups["neg"])
        if not pos_candidates or not neg_candidates:
            continue
        random.shuffle(pos_candidates)
        random.shuffle(neg_candidates)
        pos = take_existing_genomes(pos_candidates, path_map, max_per_genus)
        if not pos:
            continue
        neg = take_existing_genomes(neg_candidates, path_map, min(len(pos), max_per_genus))
        n = min(len(pos), len(neg), max_per_genus)
        if n == 0:
            continue
        selected = [(gid, "phaZ_positive") for gid in pos[:n]]
        selected += [(gid, "phaZ_negative") for gid in neg[:n]]
        for gid, status in selected:
            tax = tax_by_gid[gid]
            counts = subtype_counts.get(gid, Counter())
            major_subtype = "none"
            if counts:
                major_subtype = counts.most_common(1)[0][0]
            rows.append(
                {
                    "genome_id": gid,
                    "phaZ_status": status,
                    "phaZ_validated_count": str(pos_counts.get(gid, 0)),
                    "major_subtype": major_subtype,
                    "domain": tax.get("domain", ""),
                    "phylum": tax.get("phylum", ""),
                    "class": tax.get("class", ""),
                    "order": tax.get("order", ""),
                    "family": tax.get("family", ""),
                    "genus": genus,
                    "species": tax.get("species", ""),
                    "genome_path": path_map[gid],
                }
            )

    if pilot > 0:
        rows = rows[:pilot]
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_done(output_file: Path) -> set[str]:
    if not output_file.exists() or output_file.stat().st_size == 0:
        return set()
    done = set()
    with open(output_file, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status") == "ok":
                done.add(row["genome_id"])
    return done


def parse_fasta_gz(path: str):
    name = None
    chunks = []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks).upper()


def predict_genes(genome_id: str, genome_path: str, work_dir: Path) -> tuple[Path, Path, int]:
    cds_file = work_dir / f"{genome_id}.cds.fna"
    protein_file = work_dir / f"{genome_id}.faa"
    finder = pyrodigal.GeneFinder(meta=True)
    n_genes = 0
    with open(cds_file, "w", encoding="utf-8") as cds_out, open(protein_file, "w", encoding="utf-8") as prot_out:
        for contig_id, seq in parse_fasta_gz(genome_path):
            if len(seq) < 90:
                continue
            genes = finder.find_genes(seq.encode())
            for gene in genes:
                protein = str(gene.translate()).rstrip("*")
                if len(protein) < 30:
                    continue
                cds = str(gene.sequence()).upper()
                if len(cds) < 90:
                    continue
                n_genes += 1
                gene_id = f"{genome_id}_{n_genes}"
                cds_out.write(f">{gene_id}\n{cds}\n")
                prot_out.write(f">{gene_id}\n{protein}\n")
    return cds_file, protein_file, n_genes


def run_diamond_ribosomal(protein_file: Path, he_file: Path, ribo_db: Path, threads: int = 1) -> int:
    hits_file = protein_file.with_suffix(".ribo_hits.tsv")
    cmd = [
        str(DIAMOND),
        "blastp",
        "-q",
        str(protein_file),
        "-d",
        str(ribo_db),
        "-o",
        str(hits_file),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "evalue",
        "qcovhsp",
        "--evalue",
        "1e-5",
        "--id",
        "25",
        "--query-cover",
        "35",
        "--max-target-seqs",
        "1",
        "--threads",
        str(threads),
        "--quiet",
    ]
    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    he_ids = set()
    if hits_file.exists():
        with open(hits_file, "r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    try:
                        if float(parts[2]) >= 25 and float(parts[3]) <= 1e-5:
                            he_ids.add(parts[0])
                    except ValueError:
                        continue
    with open(he_file, "w", encoding="utf-8") as handle:
        for gene_id in sorted(he_ids):
            handle.write(gene_id + "\n")
    return len(he_ids)


def run_hmm_ribosomal(
    protein_file: Path,
    he_file: Path,
    ribosomal_hmm: Path,
    threads: int = 1,
    use_trusted_cutoffs: bool = True,
    evalue: float = 1e-5,
) -> int:
    hits_file = protein_file.with_suffix(".ribo_hmm_hits.tblout")
    cmd = [
        str(HMMSEARCH),
        "--noali",
        "--cpu",
        str(threads),
        "--tblout",
        str(hits_file),
    ]
    if use_trusted_cutoffs:
        cmd.append("--cut_ga")
    else:
        cmd.extend(["-E", str(evalue), "--domE", str(evalue)])
    cmd.extend([str(ribosomal_hmm), str(protein_file)])
    proc = subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        msg = (proc.stderr or "").strip().replace("\n", " ")[:300]
        raise RuntimeError(f"hmmsearch_failed:{msg}")

    he_ids = set()
    if hits_file.exists():
        with open(hits_file, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    full_evalue = float(parts[4])
                except ValueError:
                    continue
                if use_trusted_cutoffs or full_evalue <= evalue:
                    he_ids.add(parts[0])
    with open(he_file, "w", encoding="utf-8") as handle:
        for gene_id in sorted(he_ids):
            handle.write(gene_id + "\n")
    return len(he_ids)


def run_grodon(
    row: dict[str, str],
    heg_method: str,
    marker_db: str,
    tmp_base: str,
    marker_threads: int,
    hmm_use_trusted_cutoffs: bool,
    hmm_evalue: float,
) -> dict[str, str]:
    genome_id = row["genome_id"]
    work_dir = Path(tmp_base) / genome_id
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        cds_file, protein_file, n_genes = predict_genes(genome_id, row["genome_path"], work_dir)
        he_file = work_dir / f"{genome_id}.he_ids.txt"
        if heg_method == "diamond":
            n_he = run_diamond_ribosomal(protein_file, he_file, Path(marker_db), threads=marker_threads)
        else:
            n_he = run_hmm_ribosomal(
                protein_file,
                he_file,
                Path(marker_db),
                threads=marker_threads,
                use_trusted_cutoffs=hmm_use_trusted_cutoffs,
                evalue=hmm_evalue,
            )
        if n_he < 10:
            raise RuntimeError(f"too_few_ribosomal_hits:{n_he}")

        cmd = [str(RSCRIPT), str(RUN_GRODON_ONE), genome_id, str(cds_file), str(he_file)]
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=240)
        result_line = None
        for line in proc.stdout.splitlines():
            if line.startswith("RESULT\t"):
                result_line = line
                break
        if proc.returncode != 0 or result_line is None:
            msg = (proc.stderr or proc.stdout or "").strip().replace("\n", " ")[:300]
            raise RuntimeError(f"gRodon_failed:{msg}")
        parts = result_line.split("\t")
        return {
            **row,
            "heg_method": heg_method,
            "status": "ok",
            "n_cds": parts[2],
            "n_highly_expressed": parts[3],
            "CUBHE": parts[4],
            "GC": parts[5],
            "GCdiv": parts[6],
            "ConsistencyHE": parts[7],
            "CUB": parts[8],
            "CPB": parts[9],
            "FilteredSequences": parts[10],
            "doubling_time_h": parts[11],
            "lower_ci_h": parts[12],
            "upper_ci_h": parts[13],
            "growth_rate_per_h": str(0.6931471805599453 / float(parts[11])),
            "error": "",
        }
    except Exception as exc:
        return {
            **row,
            "heg_method": heg_method,
            "status": "failed",
            "n_cds": "",
            "n_highly_expressed": "",
            "CUBHE": "",
            "GC": "",
            "GCdiv": "",
            "ConsistencyHE": "",
            "CUB": "",
            "CPB": "",
            "FilteredSequences": "",
            "doubling_time_h": "",
            "lower_ci_h": "",
            "upper_ci_h": "",
            "growth_rate_per_h": "",
            "error": str(exc)[:500],
        }
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def run_grodon_task(args: tuple[dict[str, str], str, str, str, int, bool, float]) -> dict[str, str]:
    return run_grodon(*args)


def summarize(output_file: Path, summary_file: Path) -> None:
    rows = []
    with open(output_file, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status") == "ok":
                rows.append(row)
    grouped = defaultdict(lambda: {"pos": [], "neg": []})
    for row in rows:
        value = float(row["growth_rate_per_h"])
        key = row["genus"]
        if row["phaZ_status"] == "phaZ_positive":
            grouped[key]["pos"].append(value)
        else:
            grouped[key]["neg"].append(value)
    out_rows = []
    for genus, values in sorted(grouped.items()):
        pos = values["pos"]
        neg = values["neg"]
        if not pos or not neg:
            continue
        out_rows.append(
            {
                "genus": genus,
                "n_phaZ_positive": len(pos),
                "n_phaZ_negative": len(neg),
                "mean_growth_phaZ_positive": sum(pos) / len(pos),
                "mean_growth_phaZ_negative": sum(neg) / len(neg),
                "delta_positive_minus_negative": sum(pos) / len(pos) - sum(neg) / len(neg),
            }
        )
    write_tsv(
        summary_file,
        out_rows,
        [
            "genus",
            "n_phaZ_positive",
            "n_phaZ_negative",
            "mean_growth_phaZ_positive",
            "mean_growth_phaZ_negative",
            "delta_positive_minus_negative",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run same-genus gRodon growth analysis")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--marker-threads", type=int, default=1, help="HMMER/DIAMOND threads per worker")
    parser.add_argument("--max-per-genus", type=int, default=50, help="Use 0 for all matched genomes per genus")
    parser.add_argument("--pilot", type=int, default=0, help="Run first N manifest rows only")
    parser.add_argument("--include-initial-hits-as-negative", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--heg-method", choices=["hmm", "diamond"], default="hmm")
    parser.add_argument("--ribosomal-hmm", default=str(EXTERNAL_DIR / "ribosomal_pfam.hmm"))
    parser.add_argument("--ribosomal-diamond-db", default=str(EXTERNAL_DIR / "grodon_ribosomal_seed.dmnd"))
    parser.add_argument("--hmm-evalue", type=float, default=1e-5)
    parser.add_argument("--hmm-no-trusted-cutoffs", action="store_true")
    parser.add_argument("--output-label", default="", help="Optional suffix for output table names")
    parser.add_argument("--manifest-only", action="store_true", help="Write the matched genome manifest and exit")
    args = parser.parse_args()

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    label = args.output_label.strip() or ("" if args.heg_method == "diamond" else args.heg_method)
    suffix = f"_{label}" if label else ""
    pilot_suffix = "_pilot" if args.pilot else ""
    output_file = TABLES_DIR / f"grodon_growth_predictions{suffix}{pilot_suffix}.tsv"
    manifest_file = TABLES_DIR / f"grodon_growth_manifest{suffix}{pilot_suffix}.tsv"
    summary_file = TABLES_DIR / f"grodon_growth_same_genus_summary{suffix}{pilot_suffix}.tsv"

    marker_db = Path(args.ribosomal_diamond_db if args.heg_method == "diamond" else args.ribosomal_hmm)
    if not marker_db.exists():
        raise SystemExit(f"Missing {args.heg_method} marker database: {marker_db}")
    if args.heg_method == "hmm" and not HMMSEARCH.exists():
        raise SystemExit(f"Missing hmmsearch executable: {HMMSEARCH}")
    if args.heg_method == "diamond" and not DIAMOND.exists():
        raise SystemExit(f"Missing DIAMOND executable: {DIAMOND}")

    manifest = build_manifest(
        max_per_genus=sys.maxsize if args.max_per_genus == 0 else args.max_per_genus,
        pilot=args.pilot,
        include_initial_hits_as_negative=args.include_initial_hits_as_negative,
    )
    if not manifest:
        raise SystemExit("No matched same-genus genomes found")
    write_tsv(manifest_file, manifest, list(manifest[0].keys()))
    if args.manifest_only:
        n_pos = sum(1 for row in manifest if row["phaZ_status"] == "phaZ_positive")
        n_neg = sum(1 for row in manifest if row["phaZ_status"] == "phaZ_negative")
        n_genera = len({row["genus"] for row in manifest})
        print(f"Wrote {manifest_file}")
        print(f"Manifest rows: {len(manifest)} ({n_pos} phaZ+, {n_neg} phaZ-) across {n_genera} genera")
        return

    done = read_done(output_file) if args.resume else set()
    todo = [row for row in manifest if row["genome_id"] not in done]
    output_fields = list(manifest[0].keys()) + [
        "heg_method",
        "status",
        "n_cds",
        "n_highly_expressed",
        "CUBHE",
        "GC",
        "GCdiv",
        "ConsistencyHE",
        "CUB",
        "CPB",
        "FilteredSequences",
        "doubling_time_h",
        "lower_ci_h",
        "upper_ci_h",
        "growth_rate_per_h",
        "error",
    ]

    write_header = not output_file.exists() or not args.resume
    if write_header:
        with open(output_file, "w", encoding="utf-8", newline="") as handle:
            csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t").writeheader()

    tmp_base = tempfile.mkdtemp(prefix="grodon_", dir=str(PROCESSED_DIR))
    try:
        with Pool(args.threads) as pool, open(output_file, "a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
            tasks = [
                (
                    row,
                    args.heg_method,
                    str(marker_db),
                    tmp_base,
                    args.marker_threads,
                    not args.hmm_no_trusted_cutoffs,
                    args.hmm_evalue,
                )
                for row in todo
            ]
            for i, result in enumerate(pool.imap_unordered(run_grodon_task, tasks), start=1):
                writer.writerow(result)
                handle.flush()
                if i % 10 == 0:
                    print(f"Processed {i}/{len(todo)} genomes", flush=True)
    finally:
        shutil.rmtree(tmp_base, ignore_errors=True)

    summarize(output_file, summary_file)
    print(f"Wrote {output_file}")
    print(f"Wrote {summary_file}")


if __name__ == "__main__":
    main()
