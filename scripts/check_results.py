#!/usr/bin/env python3
"""Validate PHB_gtdb result files after a pipeline run."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
RESULTS = ROOT / "results"


EXPECTED_FASTA = {
    "phaz_proteins_all.fasta": 8768,
    "phaz_proteins_c95.fasta": 8731,
    "phaz_proteins_filtered.fasta": 7478,
    "phaz_proteins_validated.fasta": 6532,
    "phaz_bacillus_type.fasta": 60,
    "phaz_extracellular.fasta": 596,
    "phaz_extracellular_lemoignei.fasta": 434,
    "phaz_intracellular.fasta": 4693,
    "phaz_ralstonia.fasta": 2948,
    "phaz_bacillus_type_trim.fasta": 57,
    "phaz_extracellular_trim.fasta": 509,
    "phaz_extracellular_lemoignei_trim.fasta": 412,
    "phaz_intracellular_trim.fasta": 4424,
    "phaz_ralstonia_trim.fasta": 2776,
}

EXPECTED_TSV = {
    "phb_search_results.tsv": (7068, 16486),
    "archaea_phb_search_results.tsv": (2, 2),
}

EXPECTED_TREES = [
    "phaz_bacillus_type_tree.treefile",
    "phaz_extracellular_tree.treefile",
    "phaz_extracellular_lemoignei_tree.treefile",
    "phaz_ralstonia_tree.treefile",
    "phaz_ralstonia_tree_ft.treefile",
    "phaz_intracellular_tree_ft.treefile",
]

EXPECTED_FIGURES = [
    "length_distribution.png",
    "lipase_box_barplot.png",
    "phylum_barplot.png",
    "phylum_subtype_heatmap.png",
    "pipeline_timeline.png",
    "reference_matches.png",
    "subtype_pie.png",
    "top_genera.png",
]


def count_fasta(path: Path) -> int:
    with path.open(errors="ignore") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def check_fasta() -> list[str]:
    problems: list[str] = []
    print("FASTA")
    for name, expected in EXPECTED_FASTA.items():
        path = DATA / name
        if not path.exists():
            problems.append(f"missing FASTA: {path}")
            print(f"  FAIL {name}: missing")
            continue
        if path.stat().st_size == 0:
            problems.append(f"empty FASTA: {path}")
        actual = count_fasta(path)
        status = "OK" if actual == expected else "FAIL"
        if actual != expected:
            problems.append(f"{name}: expected {expected}, got {actual}")
        print(f"  {status} {name}: {actual} sequences, {path.stat().st_size} bytes")
    return problems


def check_tsv() -> list[str]:
    problems: list[str] = []
    print("TSV")
    for name, (expected_rows, expected_count) in EXPECTED_TSV.items():
        path = DATA / name
        if not path.exists():
            problems.append(f"missing TSV: {path}")
            print(f"  FAIL {name}: missing")
            continue
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        total = sum(int(row.get("phaZ_count", "0") or 0) for row in rows)
        status = "OK" if (len(rows), total) == (expected_rows, expected_count) else "FAIL"
        if status == "FAIL":
            problems.append(
                f"{name}: expected rows/count {(expected_rows, expected_count)}, "
                f"got {(len(rows), total)}"
            )
        print(f"  {status} {name}: rows={len(rows)}, phaZ_count_sum={total}")
    return problems


def check_outputs() -> list[str]:
    problems: list[str] = []
    print("TREES")
    for name in EXPECTED_TREES:
        path = RESULTS / name
        status = "OK" if path.exists() and path.stat().st_size > 0 else "FAIL"
        if status == "FAIL":
            problems.append(f"missing/empty tree: {path}")
        print(f"  {status} {name}")

    print("FIGURES")
    for name in EXPECTED_FIGURES:
        path = RESULTS / "figures" / name
        status = "OK" if path.exists() and path.stat().st_size > 0 else "FAIL"
        if status == "FAIL":
            problems.append(f"missing/empty figure: {path}")
        print(f"  {status} {name}")
    return problems


def main() -> int:
    print(f"Checking {ROOT}")
    problems = []
    problems.extend(check_fasta())
    problems.extend(check_tsv())
    problems.extend(check_outputs())

    if problems:
        print("\nPROBLEMS")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print("\nAll expected result files passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
