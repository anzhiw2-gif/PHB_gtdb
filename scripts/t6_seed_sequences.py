#!/usr/bin/env python3
"""T6: Select representative seed sequences for environmental metagenome screening."""
import subprocess
from pathlib import Path
from collections import defaultdict
from io import StringIO
from Bio import SeqIO

PHB_DIR = Path("/home/data/haoyu/PHB_gtdb")
VALIDATED = PHB_DIR / "data/processed/phaz_proteins_validated.fasta"
OUT_FASTA = PHB_DIR / "data/processed/phaz_metagenome_seeds.fasta"
OUT_TABLE = PHB_DIR / "results/tables/phaz_metagenome_seeds.tsv"

REF_MAP = {
    "BAA33394.1":"intracellular.cupriavidus_like","CAJ93939.1":"intracellular.cupriavidus_like",
    "CAJ95805.1":"intracellular.cupriavidus_like","UCA14981.1":"intracellular.ralstonia_like",
    "WKZ88401.1":"intracellular.ralstonia_like","WP_128854079.1":"intracellular.bacillus_like",
    "P52090.1":"extracellular.lemoignei_like","WP_207907290.1":"extracellular.lemoignei_like",
    "WP_243656647.1":"extracellular.lemoignei_like","BAA19791.1":"extracellular.general",
    "AAA87070.1":"extracellular.general","BAA35137.1":"extracellular.general",
    "BAA32541.1":"extracellular.general","AAB02914.1":"extracellular.general",
}

# Step 1: Parse validated sequences by subtype
subtype_seqs = defaultdict(list)
with open(VALIDATED) as f:
    for rec in SeqIO.parse(f, "fasta"):
        ref = "unknown"
        for k in REF_MAP:
            if k in rec.description: ref = k; break
        st = REF_MAP.get(ref, "unknown")
        subtype_seqs[st].append(rec)

print("=== T6: Metagenome Seed Sequence Selection ===")
print(f"  Input: {VALIDATED}")
print(f"  Total validated sequences: {sum(len(v) for v in subtype_seqs.values()):,}")

# Step 2: CD-HIT c70 per subtype (within each subtype)
# CD-HIT reduces to ~70% identity, preserving functional diversity
all_seeds = []
for st in ["intracellular.cupriavidus_like","intracellular.ralstonia_like",
            "intracellular.bacillus_like","extracellular.general","extracellular.lemoignei_like"]:
    seqs = subtype_seqs.get(st, [])
    if not seqs:
        continue
    
    # Write subtype FASTA
    subtype_fa = PHB_DIR / f"data/processed/tmp_{st.replace('.','_')}.fasta"
    with open(subtype_fa, "w") as f:
        SeqIO.write(seqs, f, "fasta")
    
    # CD-HIT c70
    out_fa = PHB_DIR / f"data/processed/tmp_{st.replace('.','_')}_c70.fasta"
    cmd = f"cd-hit -i {subtype_fa} -o {out_fa} -c 0.70 -n 4 -M 8000 -T 4 -d 0"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Count results
    n_input = len(seqs)
    n_output = sum(1 for _ in open(out_fa) if _.startswith(">")) if out_fa.exists() else 0
    
    # If n_output > 15, take top 15 longest
    cluster_seqs = list(SeqIO.parse(out_fa, "fasta")) if out_fa.exists() else seqs[:10]
    if len(cluster_seqs) > 15:
        cluster_seqs.sort(key=lambda x: len(x.seq), reverse=True)
        cluster_seqs = cluster_seqs[:15]
    
    for rec in cluster_seqs:
        rec.description = f"{rec.description} [subtype={st}]"
        all_seeds.append(rec)
    
    print(f"  {st:35s}: {n_input:5d} -> {n_output:5d} (c70) -> {len(cluster_seqs):3d} (seeds)")
    
    # Cleanup
    subtype_fa.unlink(missing_ok=True)
    out_fa.unlink(missing_ok=True)
    (Path(str(out_fa) + ".clstr")).unlink(missing_ok=True)

# Step 3: Write seed FASTA + table
with open(OUT_FASTA, "w") as f:
    SeqIO.write(all_seeds, f, "fasta")

with open(OUT_TABLE, "w") as f:
    f.write("seq_id\tsubtype\tlength\tgenome\tphylum\tgenus\n")
    for rec in all_seeds:
        parts = rec.description.split()[0].split("|")
        gid = parts[0] if len(parts)>0 else "?"
        phy = parts[1] if len(parts)>1 else "?"
        gen = parts[2] if len(parts)>2 else "?"
        st  = "unknown"
        for k in REF_MAP:
            if k in rec.description: st = REF_MAP.get(k,"unknown"); break
        f.write(f"{rec.id}\t{st}\t{len(rec.seq)}\t{gid}\t{phy}\t{gen}\n")

print(f"\n  Total seed sequences: {len(all_seeds)}")
print(f"  Saved: {OUT_FASTA}")
print(f"  Saved: {OUT_TABLE}")
print(f"\n  Next: These seeds can be used as DIAMOND/MMseqs2 query")
print(f"  for Chapter 4 environmental metagenome screening.")
print("=" * 55)
