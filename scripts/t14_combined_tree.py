import subprocess, random, json
from pathlib import Path
from Bio import SeqIO

PHB = Path("/home/data/haoyu/PHB_gtdb")
OUT_DIR = PHB / "results/tree/combined"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CONDA = "/home/data/haoyu/miniconda3/envs/phb_gtdb/bin"

ST_FILES = {
    "intra_cupriavidus": "phaz_intracellular_trim.fasta",
    "intra_ralstonia": "phaz_ralstonia_trim.fasta",
    "intra_bacillus": "phaz_bacillus_type_trim.fasta",
    "extra_general": "phaz_extracellular_trim.fasta",
    "extra_lemoignei": "phaz_extracellular_lemoignei_trim.fasta",
}

# Sampling plan: small subtypes ALL, large subtypes sampled
SAMPLE_PLAN = {
    "intra_cupriavidus": 250,   # 4424 total -> 5.7%
    "intra_ralstonia": 250,     # 2776 total -> 9.0%
    "intra_bacillus": 0,        # 0 = use ALL (57 total)
    "extra_general": 150,       # 509 total -> 29.5%
    "extra_lemoignei": 150,     # 412 total -> 36.4%
}

random.seed(42)
sampled = []
sample_log = []
for st_tag, fname in ST_FILES.items():
    path = PHB / "data/processed" / fname
    if not path.exists(): continue
    recs = list(SeqIO.parse(str(path), "fasta"))
    n_total = len(recs)
    n_target = SAMPLE_PLAN[st_tag]
    if n_target == 0 or n_target >= n_total:
        picked = recs  # use all
        n_sample = n_total
    else:
        n_sample = min(n_target, n_total)
        picked = random.sample(recs, n_sample)
    for rec in picked:
        rec.id = f"{st_tag}|{rec.id}"
        rec.description = ""
        sampled.append(rec)
    sample_log.append({"subtype": st_tag, "total": n_total, "sampled": n_sample,
                        "pct": round(n_sample/n_total*100, 1), "reason": "all" if n_sample==n_total else "proportional"})
    print(f"  {st_tag:25s}: {n_total:5d} -> {n_sample:4d} ({n_sample/n_total*100:.1f}%)")

total = len(sampled)
print(f"  Total sampled: {total}")
print(f"  Expected MAFFT time: ~{total//50} min")

# Save sampling log
log_path = OUT_DIR / "sampling_log.json"
with open(log_path, "w") as f:
    json.dump(sample_log, f, indent=2)

# Write combined FASTA
fa = PHB / "data/processed/phaz_combined_sample.fasta"
with open(fa, "w") as f:
    SeqIO.write(sampled, f, "fasta")

# MAFFT
print("  MAFFT aligning...")
mafft_out = PHB / "data/processed/phaz_combined_sample_aligned.fasta"
subprocess.run(f"{CONDA}/mafft --auto --thread 8 --quiet {fa} > {mafft_out}", shell=True, check=True)
n_aligned = sum(1 for l in open(mafft_out) if l.startswith(">"))
print(f"  MAFFT done: {n_aligned} seqs")

# trimAl
trim_out = PHB / "data/processed/phaz_combined_sample_trim.fasta"
subprocess.run(f"{CONDA}/trimal -in {mafft_out} -out {trim_out} -automated1", shell=True)
print(f"  trimAl done")

# FastTree
tree_out = OUT_DIR / "phaz_combined_sample_ft.treefile"
cmd = f"{CONDA}/FastTree -lg -gamma -quiet {trim_out}"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
tree_out.write_text(result.stdout)
print(f"  FastTree done: {tree_out.stat().st_size} bytes")

# Monophyly check
newick = tree_out.read_text()
print(f"\n  === Monophyly Assessment ===")
for st_tag in ST_FILES:
    count = newick.count(st_tag)
    rating = "★" * min(5, max(1, 6 - count//3))
    print(f"    {st_tag:25s}: {count:3d} occurrences {rating}")

print(f"\n  Tree saved: {tree_out}")
print(f"  Log saved: {log_path}")
