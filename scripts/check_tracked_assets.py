#!/usr/bin/env python3
"""Validate tracked assets in the GitHub repository (no server data)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Tracked figure data
EXPECTED_TABLES = {
    "phaz_validated_genome_protein_summary.tsv": ["metric", "value"],
    "phaz_validated_copy_number.tsv": ["copies", "n_genomes", "pct"],
    "phaz_validated_subtype_count.tsv": ["subtype", "sequences", "pct", "genomes"],
    "phaz_validated_phylum_count.tsv": ["phylum", "sequences"],
    "phaz_signal_peptide_summary.tsv": ["subtype", "total", "signal_positive", "pct"],
    "phaz_phac_cross_reference.tsv": ["genus", "in_phaz_set", "phaz_count", "known_phac"],
    "phaz_metagenome_seeds.tsv": ["seq_id", "subtype", "length"],
}

# Nature figures (8 total)
EXPECTED_FIGURES = [
    "figure1_workflow_funnel",
    "figure2_phylum_heatmap",
    "figure3_subtype_lipase",
    "figure4_genera_phylogeny",
    "figure5_grodon_growth_comparison_hmm_allmatched",
    "figure26_copy_number_phac",
    "figure27_structures_conservation",
    "figure28_signal_seeds",
]

# Tree files in results/tree/
EXPECTED_TREES = {
    "bacillus_type": "phaz_bacillus_type_tree.treefile",
    "extracellular": "phaz_extracellular_tree.treefile",
    "extracellular_lemoignei": "phaz_extracellular_lemoignei_tree.treefile",
    "intracellular": "phaz_intracellular_tree_ft.treefile",
    "ralstonia": "phaz_ralstonia_tree_ft.treefile",
}

# PDB structures
EXPECTED_PDB_MIN = 20  # minimum expected

# Docs
EXPECTED_DOCS = [
    "../README.md", "../RUN_MANIFEST.md", "SAMPLING_STRATEGY.md",
    "GENOME_VS_PROTEIN_COUNTS.md", "DATA_ACCOUNTING.md",
    "PROJECT_STATUS.md", "PRESENTATION_SUMMARY.md",
    "METHODS.md", "RESULTS.md", "FINAL_REPORT.md",
    "FIGURE_CAPTIONS.md", "PIPELINE.md", "SCRIPT_INDEX.md",
]

problems = []

# Check tables
print("=== Tracked Tables (figure_data/) ===")
for name, expected_cols in EXPECTED_TABLES.items():
    path = ROOT / "figure_data" / name
    if not path.exists():
        problems.append(f"Missing: {name}")
        print(f"  FAIL {name}: missing")
        continue
    header = path.read_text().split("\n")[0].split("\t")
    missing_cols = [c for c in expected_cols if c not in header]
    if missing_cols:
        problems.append(f"{name}: missing columns {missing_cols}")
        print(f"  FAIL {name}: missing columns {missing_cols}")
    else:
        lines = len([l for l in path.read_text().split("\n") if l.strip()]) - 1
        print(f"  OK   {name}: {lines} data rows")

# Check figures
print("\n=== Nature Figures (figures/nature/) ===")
for name in EXPECTED_FIGURES:
    missing = []
    for ext in ["pdf", "png", "svg"]:
        path = ROOT / "figures" / "nature" / f"{name}.{ext}"
        if not path.exists():
            missing.append(ext)
    if missing:
        problems.append(f"Figure {name}: missing {missing}")
        print(f"  FAIL {name}: missing formats {missing}")
    else:
        print(f"  OK   {name}: pdf+png+svg")

# Check trees
print("\n=== Phylogenetic Trees (results/tree/) ===")
for subtype, fname in EXPECTED_TREES.items():
    path = ROOT / "results" / "tree" / subtype / fname
    if path.exists() and path.stat().st_size > 0:
        print(f"  OK   {subtype}: {fname} ({path.stat().st_size} bytes)")
    else:
        problems.append(f"Tree missing: {subtype}/{fname}")
        print(f"  FAIL {subtype}: {fname} missing")

# Check PDBs
print("\n=== 3D Structures (results/structures/) ===")
pdb_dir = ROOT / "results" / "structures"
n_pdb = len(list(pdb_dir.glob("*.pdb")))
if n_pdb >= EXPECTED_PDB_MIN:
    print(f"  OK   {n_pdb} PDB files (min {EXPECTED_PDB_MIN})")
else:
    problems.append(f"Only {n_pdb} PDBs, expected >= {EXPECTED_PDB_MIN}")
    print(f"  FAIL {n_pdb} PDBs (need >= {EXPECTED_PDB_MIN})")

# Check docs
print("\n=== Documentation (docs/) ===")
for name in EXPECTED_DOCS:
    path = ROOT / "docs" / name
    if path.exists():
        print(f"  OK   {name}")
    else:
        if name not in ["PROJECT_STATUS.md"]:  # PROJECT_STATUS is in docs, copied from memory
            problems.append(f"Doc missing: {name}")
            print(f"  FAIL {name}: missing")

# Summary
print(f"\n{'='*50}")
if problems:
    print(f"PROBLEMS ({len(problems)}):")
    for p in problems:
        print(f"  - {p}")
else:
    print("All tracked assets OK.")
print(f"\nNOTE: For full server-side validation, run check_results.py on T141.")
