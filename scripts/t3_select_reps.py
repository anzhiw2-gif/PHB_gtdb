#!/usr/bin/env python3
"""T3a: Select representative PhaZ sequences for 3D structure prediction.
Strategy: pick 1 EXTRA and 1 INTRA representative per subtype (10 total),
or 2 extremes if CPU time is limited."""
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

VALIDATED = Path("/home/data/haoyu/PHB_gtdb/data/processed/phaz_proteins_validated.fasta")
OUT_FASTA = Path("/home/data/haoyu/PHB_gtdb/data/processed/phaz_3d_representatives.fasta")

REF_MAP = {
    "BAA33394.1":"intracellular.cupriavidus_like","CAJ93939.1":"intracellular.cupriavidus_like",
    "CAJ95805.1":"intracellular.cupriavidus_like","UCA14981.1":"intracellular.ralstonia_like",
    "WKZ88401.1":"intracellular.ralstonia_like","WP_128854079.1":"intracellular.bacillus_like",
    "P52090.1":"extracellular.lemoignei_like","WP_207907290.1":"extracellular.lemoignei_like",
    "WP_243656647.1":"extracellular.lemoignei_like","BAA19791.1":"extracellular.general",
    "AAA87070.1":"extracellular.general","BAA35137.1":"extracellular.general",
    "BAA32541.1":"extracellular.general","AAB02914.1":"extracellular.general",
}

# Parse and group
subtype_seqs = defaultdict(list)
with open(VALIDATED) as f:
    for rec in SeqIO.parse(f, "fasta"):
        ref = "unknown"
        for k in REF_MAP:
            if k in rec.description: ref = k; break
        st = REF_MAP.get(ref, "unknown")
        subtype_seqs[st].append(rec)

# Selection strategy:
# Phase 1 (CPU pilot): 2 extremes — Cupriavidus intracellular + Lemoignei extracellular
# Phase 2 (if time): add 1 per remaining subtype
selected = []

# Phase 1: Pick the longest, highest-confidence (best pident) from each extreme
for st, label in [
    ("intracellular.cupriavidus_like", "Intra_Cupriavidus"),
    ("extracellular.lemoignei_like", "Extra_Lemoignei"),
]:
    seqs = subtype_seqs.get(st, [])
    if not seqs: continue
    # Sort by sequence length (longer = more complete structure)
    seqs.sort(key=lambda x: len(x.seq), reverse=True)
    best = seqs[0]
    # Extract pident from header
    import re
    pident_match = re.search(r'pident=([\d.]+)%', best.description)
    pident = float(pident_match.group(1)) if pident_match else 0
    best.id = f"{label}|{best.id}"
    best.description = f"{best.description} [selected_for_3D]"
    selected.append(best)
    print(f"  {label:25s}: {best.id} | {len(best.seq)}aa | pident={pident:.1f}%")

# Phase 2: Add 1 representative from each remaining subtype
for st, label in [
    ("intracellular.ralstonia_like", "Intra_Ralstonia"),
    ("intracellular.bacillus_like", "Intra_Bacillus"),
    ("extracellular.general", "Extra_General"),
]:
    seqs = subtype_seqs.get(st, [])
    if not seqs: continue
    seqs.sort(key=lambda x: len(x.seq), reverse=True)
    best = seqs[0]
    best.id = f"{label}|{best.id}"
    selected.append(best)
    print(f"  {label:25s}: {best.id} | {len(best.seq)}aa")

# Write
with open(OUT_FASTA, "w") as f:
    SeqIO.write(selected, f, "fasta")

print(f"\n  Total sequences selected: {len(selected)}")
print(f"  Saved: {OUT_FASTA}")
print(f"\n  Phase 1 (run first): {selected[0].id}, {selected[1].id}")
print(f"  Phase 2 (if CPU time permits): remaining {len(selected)-2} sequences")
