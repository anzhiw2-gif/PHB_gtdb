#!/usr/bin/env python3
"""Prepare a ribosomal-protein Pfam HMM subset for gRodon HEG detection.

The formal gRodon workflow needs a reproducible set of highly expressed genes.
Here we use ribosomal-protein profile HMMs selected from Pfam metadata, then
extract only those models from the full Pfam-A HMM file.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DIR = PROJECT_ROOT / "data" / "external" / "grodon"
PFAM_DIR = EXTERNAL_DIR / "pfam"
PFAM_DAT_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def download(url: str, path: Path, force: bool = False) -> None:
    if path.exists() and path.stat().st_size > 0 and not force:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response, open(tmp, "wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(path)


def is_prokaryotic_ribosomal_model(record: dict[str, str], keywords: list[str]) -> bool:
    pfam_id = record.get("id", "")
    haystack = f"{pfam_id} {record.get('de', '')}".lower()
    keyword_lc = [k.lower() for k in keywords]
    if not any(k in haystack for k in keyword_lc):
        return False

    excluded_terms = (
        "40s",
        "60s",
        "39s",
        "39-s",
        "28s",
        "54s",
        "mitoribosomal",
        "mitochondrial",
        "plastid",
        "chloroplast",
        "eukaryotic",
        "archaeal",
        "sirt1",
        "methyltransferase",
        "ribosome biogenesis",
        "domain 2-like",
        "sigma 54 modulation",
        "s27a",
    )
    if any(t in haystack for t in excluded_terms):
        return False
    if re.search(r"ribosomal_[ls]\d+e(?:_|$)", pfam_id.lower()):
        return False
    if re.search(r"ribosomal_[ls](?:25|27|30|41|44)(?:_|$)", pfam_id.lower()):
        return False

    if "30s ribosomal protein" in haystack or "50s ribosomal protein" in haystack:
        return True
    if re.search(r"ribosomal protein [ls]\d+p\b", haystack):
        return True
    if re.search(r"ribosomal protein [ub][ls]\d+\b", haystack):
        return True

    prokaryotic_id_patterns = (
        r"^bL\d",
        r"^Ribosomal_[LS]\d+[A-Z]?$",
        r"^Ribosomal_[LS]\d+_[CN]$",
        r"^Ribosomal_uL\d+",
        r"^Ribosom_S\d+",
        r"^Ribosomal_TL5",
        r"^RL10_C$",
        r"^Thx$",
        r"^HR_L37$",
    )
    if any(re.search(pattern, pfam_id) for pattern in prokaryotic_id_patterns):
        if not re.search(r"(?:^|_)[LS]\d+e(?:_|$)", pfam_id):
            return True

    return False


def parse_pfam_dat(dat_gz: Path, keywords: list[str]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    with gzip.open(dat_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "//":
                if current:
                    if is_prokaryotic_ribosomal_model(current, keywords):
                        records.append(current)
                current = {}
                continue
            if line.startswith("#=GF ID"):
                current["id"] = line.split(maxsplit=2)[2]
            elif line.startswith("#=GF AC"):
                current["accession"] = line.split(maxsplit=2)[2].split(".")[0]
            elif line.startswith("#=GF DE"):
                current["de"] = line.split(maxsplit=2)[2]
    return [r for r in records if r.get("accession")]


def write_accessions(records: list[dict[str, str]], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("accession\tpfam_id\tdescription\n")
        for record in records:
            handle.write(f"{record['accession']}\t{record.get('id', '')}\t{record.get('de', '')}\n")


def extract_hmms(pfam_hmm_gz: Path, accessions: set[str], out_hmm: Path) -> int:
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    keep = False
    block: list[str] = []
    with gzip.open(pfam_hmm_gz, "rt", encoding="utf-8", errors="replace") as handle, open(
        out_hmm, "w", encoding="utf-8"
    ) as out:
        for line in handle:
            if line.startswith("HMMER"):
                keep = False
                block = [line]
                continue
            block.append(line)
            if line.startswith("ACC"):
                acc = line.split()[1].split(".")[0]
                keep = acc in accessions
            if line.startswith("//"):
                if keep:
                    out.writelines(block)
                    n_written += 1
                keep = False
                block = []
    return n_written


def maybe_hmmpress(out_hmm: Path) -> None:
    hmmpress = shutil.which("hmmpress")
    if not hmmpress:
        return
    subprocess.run([hmmpress, "-f", str(out_hmm)], check=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare ribosomal Pfam HMM subset")
    parser.add_argument("--pfam-dir", default=str(PFAM_DIR))
    parser.add_argument("--output-hmm", default=str(EXTERNAL_DIR / "ribosomal_pfam.hmm"))
    parser.add_argument("--keyword", action="append", default=["ribosomal protein"])
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--skip-hmmpress", action="store_true")
    args = parser.parse_args()

    pfam_dir = Path(args.pfam_dir)
    dat_gz = pfam_dir / "Pfam-A.hmm.dat.gz"
    hmm_gz = pfam_dir / "Pfam-A.hmm.gz"
    out_hmm = Path(args.output_hmm)
    accession_table = out_hmm.with_name(out_hmm.stem + "_accessions.tsv")

    download(PFAM_DAT_URL, dat_gz, force=args.force_download)
    records = parse_pfam_dat(dat_gz, args.keyword)
    if not records:
        raise SystemExit("No Pfam ribosomal-protein HMMs matched the selected keywords")
    write_accessions(records, accession_table)
    print(f"Selected {len(records)} Pfam HMM accessions: {accession_table}", flush=True)

    download(PFAM_HMM_URL, hmm_gz, force=args.force_download)
    n_written = extract_hmms(hmm_gz, {r["accession"] for r in records}, out_hmm)
    if n_written < 10:
        raise SystemExit(f"Only extracted {n_written} HMMs, which is unexpectedly low")
    print(f"Wrote {n_written} HMM models: {out_hmm}", flush=True)

    if not args.skip_hmmpress:
        maybe_hmmpress(out_hmm)


if __name__ == "__main__":
    sys.exit(main())
