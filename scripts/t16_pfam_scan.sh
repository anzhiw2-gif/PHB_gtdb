#!/bin/bash
# T16: Pfam domain annotation of validated PhaZ via HMMER
set -e
PFAM_DIR=/home/data/haoyu/GTDB/pfam
PFAM_HMM=$PFAM_DIR/Pfam-A.hmm
PHB_DIR=/home/data/haoyu/PHB_gtdb
INPUT=$PHB_DIR/data/processed/phaz_proteins_validated.fasta
OUTPUT=$PHB_DIR/results/tables/phaz_pfam_domains.tsv
TMP=$PHB_DIR/data/processed/phaz_pfam_scan.tsv

mkdir -p $PFAM_DIR $PHB_DIR/results/tables

# Step 1: Download Pfam-A if not present
if [ ! -f "$PFAM_HMM" ]; then
    echo "Downloading Pfam-A HMM..."
    wget -q --show-progress -O $PFAM_HMM.gz \
        ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip $PFAM_HMM.gz
else
    echo "Pfam-A HMM exists"
fi

# Step 2: Press HMM database
if [ ! -f "$PFAM_HMM.h3m" ]; then
    echo "Pressing HMM database..."
    /home/data/haoyu/miniconda3/envs/phb_gtdb/bin/hmmpress $PFAM_HMM
fi

# Step 3: Run hmmscan
echo "Running hmmscan on $(grep -c '>' $INPUT) sequences (16 threads)..."
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/hmmscan \
    --cpu 16 --domtblout $TMP -E 1e-5 $PFAM_HMM $INPUT

echo "Scan done: $(wc -l < $TMP) lines"

# Step 4: Parse
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python3 << 'PYEOF'
from collections import defaultdict, Counter
domains = defaultdict(list)
with open("/home/data/haoyu/PHB_gtdb/data/processed/phaz_pfam_scan.tsv") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.split()
        if len(parts) < 23: continue
        seq_id = parts[3]
        pfam_acc = parts[1]
        pfam_name = parts[0]
        evalue = float(parts[12])
        domains[seq_id].append((pfam_acc, pfam_name, evalue))

pfam_counts = Counter()
for seq_id, dom_list in domains.items():
    for acc, name, e in dom_list:
        pfam_counts[(acc, name)] += 1

print(f"Seqs with Pfam hits: {len(domains)}/6532")
print(f"Unique Pfam families: {len(pfam_counts)}")
print("\nTop 20 Pfam domains:")
for (acc, name), cnt in pfam_counts.most_common(20):
    pct = cnt / 6532 * 100
    print(f"  {acc:15s} {name:35s} {cnt:5d} ({pct:4.1f}%)")

outpath = "/home/data/haoyu/PHB_gtdb/results/tables/phaz_pfam_domains.tsv"
with open(outpath, "w") as out:
    out.write("seq_id\tpfam_acc\tpfam_name\tevalue\n")
    for seq_id, dom_list in sorted(domains.items()):
        for acc, name, e in dom_list:
            out.write(f"{seq_id}\t{acc}\t{name}\t{e:.1e}\n")
print(f"\nSaved: {outpath}")
PYEOF

echo "T16 complete!"
