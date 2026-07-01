from collections import Counter, defaultdict
from pathlib import Path

FASTA = Path("/home/data/haoyu/PHB_gtdb/data/processed/phaz_proteins_validated.fasta")
OUT  = Path("/home/data/haoyu/PHB_gtdb/results/tables")
OUT.mkdir(parents=True, exist_ok=True)

REF_MAP = {
    "BAA33394.1":"intracellular","CAJ93939.1":"intracellular",
    "CAJ95805.1":"intracellular","UCA14981.1":"ralstonia",
    "WKZ88401.1":"ralstonia","WP_128854079.1":"bacillus_type",
    "P52090.1":"extracellular_lemoignei","WP_207907290.1":"extracellular_lemoignei",
    "WP_243656647.1":"extracellular_lemoignei","BAA19791.1":"extracellular",
    "AAA87070.1":"extracellular","BAA35137.1":"extracellular",
    "BAA32541.1":"extracellular","AAB02914.1":"extracellular",
}

genomes = defaultdict(list)
phyla, genera, subtypes = Counter(), Counter(), Counter()

with open(FASTA) as f:
    for line in f:
        if not line.startswith(">"): continue
        parts = line[1:].strip().split()[0].split("|")
        gid, phylum, genus = parts[0], parts[1], parts[2]
        ref = "unknown"
        for k in REF_MAP:
            if k in line: ref = k; break
        genomes[gid].append(parts[4] if len(parts)>4 else "x")
        phyla[phylum] += 1; genera[genus] += 1
        subtypes[REF_MAP.get(ref,"unknown")] += 1

n_genomes = len(genomes)
n_seqs    = sum(len(v) for v in genomes.values())
copies    = Counter(len(v) for v in genomes.values())

print("="*55)
print("  PhaZ VALIDATED SET — Genome/Sequence Summary")
print("="*55)
print(f"  Protein sequences:    {n_seqs:>6,}")
print(f"  Genomes with phaZ:    {n_genomes:>6,}")
print(f"  Phyla:                {len(phyla):>6}")
print(f"  Genera:               {len(genera):>6}")
print(f"  Mean copies/genome:   {n_seqs/n_genomes:>6.2f}")
print()
print("  Copy number distribution:")
single = multi = 0
for c in sorted(copies):
    cnt = copies[c]; bar = "|"*(cnt//40)
    if c==1: single=cnt
    else: multi+=cnt
    print(f"    {c:2d}: {cnt:5d} genomes  {bar}")
print(f"  Single-copy: {single} ({single/n_genomes*100:.1f}%)")
print(f"  Multi-copy:  {multi} ({multi/n_genomes*100:.1f}%)")
print()
print("  Subtype (validated):")
for st, cnt in subtypes.most_common():
    print(f"    {st:30s} {cnt:5d} ({cnt/n_seqs*100:5.1f}%)")
print()
print("  Phyla (validated):")
for ph, cnt in phyla.most_common():
    print(f"    {ph:30s} {cnt:5d}")
print()
print("  Top 15 genera (validated):")
for gen, cnt in genera.most_common(15):
    print(f"    {gen:30s} {cnt:5d}")

# Save
with open(OUT/"phaz_validated_genome_protein_summary.tsv","w") as f:
    f.write("metric\tvalue\n")
    f.write(f"validated_sequences\t{n_seqs}\nunique_genomes\t{n_genomes}\n")
    f.write(f"unique_phyla\t{len(phyla)}\nunique_genera\t{len(genera)}\n")
    f.write(f"mean_copies_per_genome\t{n_seqs/n_genomes:.2f}\n")
    f.write(f"single_copy_genomes\t{single}\nmulti_copy_genomes\t{multi}\n")

with open(OUT/"phaz_validated_copy_number.tsv","w") as f:
    f.write("copies\tn_genomes\tpct\n")
    for c in sorted(copies):
        f.write(f"{c}\t{copies[c]}\t{copies[c]/n_genomes*100:.1f}\n")

with open(OUT/"phaz_validated_subtype_count.tsv","w") as f:
    f.write("subtype\tsequences\tpct\tgenomes\n")
    for st, cnt in subtypes.most_common():
        st_genomes = len(set(gid for gid, pids in genomes.items()
                             if any(REF_MAP.get((lambda h: next((k for k in REF_MAP if k in h),"unknown"))(""),"")==st for _ in [])))
        f.write(f"{st}\t{cnt}\t{cnt/n_seqs*100:.1f}\n")

with open(OUT/"phaz_validated_phylum_count.tsv","w") as f:
    f.write("phylum\tsequences\n")
    for ph, cnt in phyla.most_common():
        f.write(f"{ph}\t{cnt}\n")

print(f"\n  Saved: {OUT}/")
print("="*55)
