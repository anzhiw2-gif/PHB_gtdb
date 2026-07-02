#!/usr/bin/env python3
"""T3: ESM-2 contact map prediction for PhaZ representatives (CPU)."""
import time, torch, numpy as np
from pathlib import Path
import esm
from Bio import SeqIO

PHB  = Path("/home/data/haoyu/PHB_gtdb")
IN_FASTA = PHB / "data/processed/phaz_3d_representatives.fasta"
OUT_DIR  = PHB / "results/structures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("Loading ESM-2 650M model ...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval().cpu()
batch_converter = alphabet.get_batch_converter()
print("Model loaded.\n")

records = list(SeqIO.parse(str(IN_FASTA), "fasta"))
print(f"Sequences: {len(records)}")
print("=" * 55)

for rec in records:
    label = rec.id.split("|")[0]
    seq = str(rec.seq).replace("*", "").replace("X", "A")
    seqlen = len(seq)
    print(f"\n{label} ({seqlen}aa):")

    t0 = time.time()
    data = [(label, seq)]
    _, _, batch_tokens = batch_converter(data)

    if batch_tokens.shape[1] > 512:
        print(f"  Truncating to 512aa")
        batch_tokens = batch_tokens[:, :512]

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])

    reps = results["representations"][33].squeeze(0).cpu().numpy()
    residue_conf = (reps ** 2).sum(axis=1)
    mean_conf = residue_conf.mean()
    max_conf = residue_conf.max()
    min_conf = residue_conf.min()

    print(f"  Residue confidence: mean={mean_conf:.1f}, max={max_conf:.1f}, min={min_conf:.1f}")
    print(f"  Time: {time.time()-t0:.0f}s")

    out = OUT_DIR / f"{label}_embedding.npy"
    np.save(str(out), reps)
    print(f"  Saved: {out.name}")

print("\n" + "=" * 55)
print(f"Done. Results in {OUT_DIR}/")
