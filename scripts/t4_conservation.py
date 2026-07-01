#!/usr/bin/env python3
"""T4: Catalytic domain conservation analysis across 5 PhaZ subtypes."""
import numpy as np
from pathlib import Path
from collections import Counter

ALIGN_DIR = Path("/home/data/haoyu/PHB_gtdb/data/processed")
OUT_DIR = Path("/home/data/haoyu/PHB_gtdb/results/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUBTYPES = {
    "intracellular_cupriavidus": "phaz_intracellular_trim.fasta",
    "intracellular_ralstonia": "phaz_ralstonia_trim.fasta",
    "intracellular_bacillus": "phaz_bacillus_type_trim.fasta",
    "extracellular_general": "phaz_extracellular_trim.fasta",
    "extracellular_lemoignei": "phaz_extracellular_lemoignei_trim.fasta",
}

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY-"

def read_alignment(path):
    seqs = []
    current = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                current = []
            else:
                current.append(line.strip())
        if current:
            seqs.append("".join(current))
    return seqs

def position_conservation(seqs):
    """Compute per-position identity to consensus."""
    n = len(seqs)
    length = max(len(s) for s in seqs)
    # Pad sequences to same length
    padded = [s.ljust(length, "-") for s in seqs]
    
    results = []
    for pos in range(length):
        chars = [s[pos] for s in padded]
        counts = Counter(chars)
        consensus = counts.most_common(1)[0][0]
        consensus_freq = counts[consensus] / n
        gap_freq = counts.get("-", 0) / n
        results.append({
            "pos": pos + 1,
            "consensus": consensus,
            "consensus_pct": round(consensus_freq * 100, 1),
            "gap_pct": round(gap_freq * 100, 1),
            "entropy": round(-sum((c/n)*np.log2(c/n) for c in counts.values() if c > 0), 3),
        })
    return results

print("=== T4: Catalytic Domain Conservation Analysis ===\n")

all_blocks = []
for name, filename in SUBTYPES.items():
    path = ALIGN_DIR / filename
    if not path.exists():
        print(f"  SKIP {name}: file not found")
        continue
    
    seqs = read_alignment(path)
    cons = position_conservation(seqs)
    
    # Find conserved blocks: contiguous positions with consensus > 60%
    blocks = []
    in_block = False
    block_start = 0
    for r in cons:
        if r["consensus_pct"] >= 60 and r["gap_pct"] < 40:
            if not in_block:
                block_start = r["pos"]
                in_block = True
        else:
            if in_block:
                blocks.append((block_start, r["pos"] - 1))
                in_block = False
    if in_block:
        blocks.append((block_start, cons[-1]["pos"]))
    
    # Top 5 most conserved single positions
    top_positions = sorted(cons, key=lambda x: x["consensus_pct"], reverse=True)[:10]
    
    print(f"--- {name} ---")
    print(f"  Sequences: {len(seqs)}")
    print(f"  Alignment length: {len(cons)}")
    print(f"  Conserved blocks (>60% consensus, <40% gaps): {len(blocks)}")
    for i, (start, end) in enumerate(blocks[:8]):
        length = end - start + 1
        avg_cons = np.mean([r["consensus_pct"] for r in cons[start-1:end]])
        print(f"    Block {i+1}: pos {start}-{end} ({length}aa, avg consensus {avg_cons:.0f}%)")
    
    print(f"  Top 10 most conserved positions:")
    for r in top_positions[:5]:
        print(f"    Pos {r['pos']:4d}: {r['consensus']} ({r['consensus_pct']:.0f}%)")
    
    # Save detailed per-position CSV
    import csv
    out_path = OUT_DIR / f"conservation_{name}.tsv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pos", "consensus", "consensus_pct", "gap_pct", "entropy"], delimiter="\t")
        w.writeheader()
        w.writerows(cons)
    
    all_blocks.append({"subtype": name, "n_seqs": len(seqs), "n_blocks": len(blocks),
                        "align_len": len(cons), "max_consensus": max(r["consensus_pct"] for r in cons)})
    print()

# Summary
print("=== Summary ===")
for b in all_blocks:
    print(f"  {b['subtype']:35s}: {b['n_seqs']:5d} seqs, {b['n_blocks']:2d} conserved blocks, max identity {b['max_consensus']:.0f}%")
print(f"\nSaved per-position CSVs to {OUT_DIR}/")
