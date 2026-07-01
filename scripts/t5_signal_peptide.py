#!/usr/bin/env python3
"""T5: Signal peptide & N-terminal hydrophobicity analysis on validated PhaZ set.
Uses Kyte-Doolittle hydrophobicity scale — well-established in literature."""
from collections import defaultdict
from pathlib import Path

FASTA = Path("/home/data/haoyu/PHB_gtdb/data/processed/phaz_proteins_validated.fasta")
OUT  = Path("/home/data/haoyu/PHB_gtdb/results/tables")
OUT.mkdir(parents=True, exist_ok=True)

# Kyte-Doolittle hydrophobicity scale
KD = {
    'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,
    'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,
    'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3,'X':0,'-':0
}

REF_MAP = {
    "BAA33394.1":"intracellular","CAJ93939.1":"intracellular",
    "CAJ95805.1":"intracellular","UCA14981.1":"ralstonia",
    "WKZ88401.1":"ralstonia","WP_128854079.1":"bacillus_type",
    "P52090.1":"extracellular_lemoignei","WP_207907290.1":"extracellular_lemoignei",
    "WP_243656647.1":"extracellular_lemoignei","BAA19791.1":"extracellular",
    "AAA87070.1":"extracellular","BAA35137.1":"extracellular",
    "BAA32541.1":"extracellular","AAB02914.1":"extracellular",
}

def signal_peptide_score(seq, window=9):
    """Detect putative Sec signal peptide: N-region (pos charged) + H-region (hydrophobic).
    Returns (has_signal, n_charged, max_hydrophobicity)."""
    n_term = seq[:30]  # signal peptides typically in first 20-30 aa
    
    # N-region: count positive charges in first 10 aa
    n_charged = sum(1 for aa in n_term[:10] if aa in 'KR')
    
    # H-region: sliding window hydrophobicity in residues 5-25
    max_hydro = -999
    for i in range(4, min(25, len(n_term) - window)):
        window_seq = n_term[i:i+window]
        hydro = sum(KD.get(aa, 0) for aa in window_seq) / window
        max_hydro = max(max_hydro, hydro)
    
    # Classical Sec signal: N-region has 1+ basic residues, H-region > 1.5 KD score
    has_signal = (n_charged >= 1 and max_hydro > 1.5)
    return has_signal, n_charged, max_hydro

# Parse sequences
subtypes = defaultdict(list)
seqs = defaultdict(list)
with open(FASTA) as f:
    seq_data = []
    for line in f:
        if line.startswith(">"):
            if seq_data:
                seqs[gid].append({"header": header, "seq": "".join(seq_data)})
            header = line[1:].strip()
            parts = header.split()[0].split("|")
            gid = parts[0]
            seq_data = []
        else:
            seq_data.append(line.strip())
    if seq_data:
        seqs[gid].append({"header": header, "seq": "".join(seq_data)})

# Analyze each sequence
results = []
for gid, proteins in seqs.items():
    for p in proteins:
        ref = "unknown"
        for key in REF_MAP:
            if key in p["header"]:
                ref = key
                break
        st = REF_MAP.get(ref, "unknown")
        has_sig, n_pos, max_hydro = signal_peptide_score(p["seq"])
        results.append({
            "genome": gid,
            "subtype": st,
            "has_signal": has_sig,
            "n_positive_first10": n_pos,
            "max_hydrophobicity": round(max_hydro, 2),
            "length": len(p["seq"]),
            "n_term_30": p["seq"][:30],
        })

# By subtype summary
from collections import Counter
st_summary = defaultdict(lambda: {"total":0, "signal":0, "hydro_scores":[]})
for r in results:
    st = r["subtype"]
    st_summary[st]["total"] += 1
    if r["has_signal"]:
        st_summary[st]["signal"] += 1
    st_summary[st]["hydro_scores"].append(r["max_hydrophobicity"])

print("=" * 60)
print("  T5: Signal Peptide Prediction (Kyte-Doolittle)")
print("=" * 60)
print(f"  Total validated sequences analyzed: {len(results):,}")
print()
print(f"  {'Subtype':30s} {'Total':>6s} {'Signal+':>8s} {'Pct':>6s} {'Mean Hydro':>10s}")
print(f"  {'-'*30} {'-'*6} {'-'*8} {'-'*6} {'-'*10}")

for st in ["intracellular","ralstonia","bacillus_type","extracellular","extracellular_lemoignei"]:
    s = st_summary[st]
    mean_hydro = sum(s["hydro_scores"])/len(s["hydro_scores"]) if s["hydro_scores"] else 0
    pct = s["signal"]/s["total"]*100 if s["total"] else 0
    print(f"  {st:30s} {s['total']:6d} {s['signal']:8d} {pct:5.1f}% {mean_hydro:10.2f}")

# Extra vs intra comparison
intra_st = ["intracellular","ralstonia","bacillus_type"]
extra_st = ["extracellular","extracellular_lemoignei"]
intra_total = sum(st_summary[s]["total"] for s in intra_st)
intra_sig   = sum(st_summary[s]["signal"] for s in intra_st)
extra_total = sum(st_summary[s]["total"] for s in extra_st)
extra_sig   = sum(st_summary[s]["signal"] for s in extra_st)

print(f"\n  Intracellular combined: {intra_sig}/{intra_total} ({intra_sig/intra_total*100:.1f}%)")
print(f"  Extracellular combined: {extra_sig}/{extra_total} ({extra_sig/extra_total*100:.1f}%)")
print(f"  Fold enrichment:        {extra_sig/extra_total/(intra_sig/intra_total) if intra_sig else 'N/A':.1f}x")
print()

# Save
with open(OUT / "phaz_signal_peptide_summary.tsv", "w") as f:
    f.write("subtype\ttotal\tsignal_positive\tpct\tmean_max_hydrophobicity\n")
    for st in ["intracellular","ralstonia","bacillus_type","extracellular","extracellular_lemoignei"]:
        s = st_summary[st]
        mean_hydro = sum(s["hydro_scores"])/len(s["hydro_scores"])
        f.write(f"{st}\t{s['total']}\t{s['signal']}\t{s['signal']/s['total']*100:.1f}\t{mean_hydro:.2f}\n")

with open(OUT / "phaz_signal_peptide_detail.tsv", "w") as f:
    f.write("genome\tsubtype\thas_signal\tn_positive_first10\tmax_hydrophobicity\tlength\tn_term_30\n")
    for r in results:
        f.write(f"{r['genome']}\t{r['subtype']}\t{r['has_signal']}\t{r['n_positive_first10']}\t{r['max_hydrophobicity']}\t{r['length']}\t{r['n_term_30']}\n")

print(f"  Saved: {OUT}/phaz_signal_peptide_summary.tsv")
print(f"  Saved: {OUT}/phaz_signal_peptide_detail.tsv")
print("=" * 60)
