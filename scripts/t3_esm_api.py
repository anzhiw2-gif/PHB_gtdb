#!/usr/bin/env python3
"""T3: Predict 3D structures via ESM Atlas API (Meta).
No local installation needed — POST sequences to esmatlas.com."""
import time, json, requests
from pathlib import Path
from Bio import SeqIO

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

PHB = Path("/home/data/haoyu/PHB_gtdb")
IN_FASTA = PHB / "data/processed/phaz_3d_representatives.fasta"
OUT_DIR = PHB / "results/structures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

records = list(SeqIO.parse(str(IN_FASTA), "fasta"))
print(f"Sequences to fold: {len(records)}")
print("=" * 55)

for rec in records:
    label = rec.id.split("|")[0]
    seq = str(rec.seq).replace("*", "").replace("X", "")
    # ESM Atlas limit ~400aa; lemoignei needs shorter truncation
    if "lemoignei" in label.lower() or "Lemoignei" in label:
        seq = seq[:200]
    elif len(seq) > 300:
        seq = seq[:300]
    out_pdb = OUT_DIR / f"{label.replace(' ','_')}.pdb"

    if out_pdb.exists():
        print(f"  {label}: SKIP (already exists)")
        continue

    print(f"  {label} ({len(seq)}aa): folding via ESM Atlas ...", end=" ", flush=True)
    t0 = time.time()

    try:
        resp = requests.post(API_URL, data=seq, timeout=300)
        if resp.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(resp.text)
            elapsed = time.time() - t0
            print(f"done ({elapsed:.0f}s, {len(resp.text)//1024}KB)")
        else:
            print(f"FAILED: HTTP {resp.status_code}")
    except Exception as e:
        print(f"FAILED: {e}")

# Verify results
print(f"\n=== Results ===")
for rec in records:
    label = rec.id.split("|")[0]
    out_pdb = OUT_DIR / f"{label.replace(' ','_')}.pdb"
    if out_pdb.exists():
        size_kb = out_pdb.stat().st_size // 1024
        # Count CA atoms
        ca = sum(1 for line in open(out_pdb) if line.startswith("ATOM") and "CA" in line)
        print(f"  {label}: {ca} residues, {size_kb}KB")
    else:
        print(f"  {label}: MISSING")

print(f"\nDone. PDBs saved to {OUT_DIR}/")
