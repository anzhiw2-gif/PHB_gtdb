#!/usr/bin/env python3
"""T3 extended: Predict 5 representatives per subtype (25 total) via ESM Atlas API."""
import time, requests
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LEN = 300  # catalytic domain truncation

PHB = Path("/home/data/haoyu/PHB_gtdb")
VALIDATED = PHB / "data/processed/phaz_proteins_validated.fasta"
OUT_DIR = PHB / "results/structures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

REF_MAP = {
    "BAA33394.1":"intra_cupriavidus","CAJ93939.1":"intra_cupriavidus",
    "CAJ95805.1":"intra_cupriavidus","UCA14981.1":"intra_ralstonia",
    "WKZ88401.1":"intra_ralstonia","WP_128854079.1":"intra_bacillus",
    "P52090.1":"extra_lemoignei","WP_207907290.1":"extra_lemoignei",
    "WP_243656647.1":"extra_lemoignei","BAA19791.1":"extra_general",
    "AAA87070.1":"extra_general","BAA35137.1":"extra_general",
    "BAA32541.1":"extra_general","AAB02914.1":"extra_general",
}

# Group sequences by subtype
subtype_seqs = defaultdict(list)
with open(VALIDATED) as f:
    for rec in SeqIO.parse(f, "fasta"):
        ref = "unknown"
        for k in REF_MAP:
            if k in rec.description: ref = k; break
        st = REF_MAP.get(ref, "unknown")
        subtype_seqs[st].append(rec)

# Select top 5 per subtype (longest = most complete)
PER_SUBTYPE = 5
selected = []
for st in ["intra_cupriavidus","intra_ralstonia","intra_bacillus",
            "extra_general","extra_lemoignei"]:
    seqs = subtype_seqs.get(st, [])
    seqs.sort(key=lambda x: len(x.seq), reverse=True)
    for i, rec in enumerate(seqs[:PER_SUBTYPE]):
        label = f"{st}_{i+1}"
        rec.id = f"{label}|{rec.id}"
        rec.description = f"{rec.description} [rep={i+1}/{st}]"
        selected.append(rec)
    print(f"  {st}: selected {min(PER_SUBTYPE, len(seqs))} from {len(seqs)}")

print(f"\nTotal to predict: {len(selected)}")
print("=" * 55)

n_new, n_skip, n_fail = 0, 0, 0
for rec in selected:
    label = rec.id.split("|")[0]
    seq = str(rec.seq).replace("*", "").replace("X", "")
    if len(seq) > MAX_LEN:
        seq = seq[:MAX_LEN]
    out_pdb = OUT_DIR / f"{label}.pdb"

    if out_pdb.exists():
        n_skip += 1
        continue

    print(f"  {label} ({len(seq)}aa): ", end="", flush=True)
    t0 = time.time()
    try:
        resp = requests.post(API_URL, data=seq, timeout=120)
        if resp.status_code == 200:
            with open(out_pdb, "w") as f:
                f.write(resp.text)
            elapsed = time.time() - t0
            ca = resp.text.count("CA  ")
            print(f"OK ({elapsed:.0f}s, {ca}CA, {len(resp.text)//1024}KB)")
            n_new += 1
        else:
            print(f"HTTP {resp.status_code}")
            n_fail += 1
    except Exception as e:
        print(f"ERR: {e}")
        n_fail += 1
    time.sleep(0.5)  # rate limit

print(f"\n=== Summary ===")
print(f"  New: {n_new}  Skipped: {n_skip}  Failed: {n_fail}")
print(f"  Total PDBs: {len(list(OUT_DIR.glob('*.pdb')))}")
