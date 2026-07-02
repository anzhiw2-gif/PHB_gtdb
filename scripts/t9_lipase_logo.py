#!/usr/bin/env python3
"""#9: Lipase box sequence logo analysis per PhaZ subtype.
Extracts G-X-S-X-G motif region from trimmed alignments."""
import numpy as np
from pathlib import Path
from collections import Counter

PHB = Path("/home/data/haoyu/PHB_gtdb")
ALIGN_DIR = PHB / "data/processed"
OUT_DIR = PHB / "results/tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "intracellular.cupriavidus_like": "phaz_intracellular_trim.fasta",
    "intracellular.ralstonia_like": "phaz_ralstonia_trim.fasta",
    "intracellular.bacillus_like": "phaz_bacillus_type_trim.fasta",
    "extracellular.general": "phaz_extracellular_trim.fasta",
    "extracellular.lemoignei_like": "phaz_extracellular_lemoignei_trim.fasta",
}

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY-"
LIPASE_MOTIFS = ["GWSMG", "GLSMG", "GSSMG", "GHSMG", "GASMG", "GLSAG", "GLSGG",
                 "GFSMG", "GMSMG", "GISSG", "GLSSG", "GQSMG", "GWSTG"]

def read_alignment(path):
    seqs = []
    cur = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else: cur.append(line.strip())
        if cur: seqs.append("".join(cur))
    return seqs

print("=== #9: Lipase Box Sequence Logo Analysis ===\n")

all_motif_data = []
for st, filename in FILES.items():
    path = ALIGN_DIR / filename
    if not path.exists():
        print(f"  SKIP {st}: no alignment")
        continue
    seqs = read_alignment(path)
    n = len(seqs)
    length = max(len(s) for s in seqs)
    padded = [s.ljust(length, "-") for s in seqs]
    
    # Search for lipase box motifs (G-X-S-X-G pattern)
    # Scan each sequence for the pentapeptide pattern
    motif_starts = Counter()
    motif_seqs_found = []
    for si, seq in enumerate(seqs):
        for pos in range(len(seq) - 4):
            pent = seq[pos:pos+5]
            # G-X-S-X-G pattern (or close variants)
            if pent[0] == 'G' and pent[2] == 'S' and pent[4] == 'G':
                motif_starts[pos] += 1
                motif_seqs_found.append({"seq_idx": si, "pos": pos, "motif": pent})
    
    n_with_motif = len(set(m["seq_idx"] for m in motif_seqs_found))
    pct = n_with_motif / n * 100 if n else 0
    
    # Find most common motif
    motif_counter = Counter(m["motif"] for m in motif_seqs_found)
    top_motifs = motif_counter.most_common(5)
    
    print(f"  {st}:")
    print(f"    Sequences: {n}, with G-X-S-X-G: {n_with_motif} ({pct:.1f}%)")
    print(f"    Top motifs: {', '.join(f'{m}({c})' for m,c in top_motifs)}")
    
    # Most common start position
    if motif_starts:
        best_pos = motif_starts.most_common(1)[0]
        print(f"    Best motif position: {best_pos[0]} ({best_pos[1]} seqs)")
    
    # Extract aligned motif region
    if motif_starts:
        best_start = motif_starts.most_common(1)[0][0]
        # Extract 5-mer at that position for each sequence
        motif_col = []
        for s in padded:
            m = s[best_start:best_start+5] if len(s) > best_start+4 else "-----"
            motif_col.append(m)
        
        # Position-wise AA frequency
        logo_data = []
        for pi in range(5):
            aa_count = Counter(m[pi] for m in motif_col)
            row = {"pos": pi + 1}
            for aa in AA_ORDER:
                row[aa] = aa_count.get(aa, 0) / len(motif_col) * 100
            logo_data.append(row)
        
        all_motif_data.append({"subtype": st, "pct": pct, "top_motif": top_motifs[0][0],
                                "n_total": n, "logo": logo_data})
    else:
        all_motif_data.append({"subtype": st, "pct": 0, "top_motif": "-",
                                "n_total": n, "logo": []})
    
    # Save per-subtype logo CSV
    import csv
    logo_out = OUT_DIR / f"lipase_logo_{st.replace('.','_')}.tsv"
    if logo_data:
        with open(logo_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pos"] + list(AA_ORDER), delimiter="\t")
            w.writeheader()
            w.writerows(logo_data)

# Summary
print(f"\n=== Subtype Comparison ===")
print(f"  {'Subtype':35s} {'N':>5s} {'G-X-S-X-G%':>10s} {'Top Motif':>10s}")
for d in all_motif_data:
    print(f"  {d['subtype']:35s} {d['n_total']:5d} {d['pct']:9.1f}% {d['top_motif']:>10s}")

with open(OUT_DIR / "lipase_box_logo_summary.tsv", "w") as f:
    f.write("subtype\tn_total\tn_with_motif\tpct\ttop_motif\n")
    for d in all_motif_data:
        f.write(f"{d['subtype']}\t{d['n_total']}\t{int(d['n_total']*d['pct']/100)}\t{d['pct']:.1f}\t{d['top_motif']}\n")
print(f"\nSaved: {OUT_DIR}/lipase_box_logo_summary.tsv")
