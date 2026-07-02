#!/usr/bin/env python3
"""Generate Nature-style supplement figures (Fig 2.6, 2.7, 2.8) for Chapter 2."""
from __future__ import annotations
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib import gridspec

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "figure_data"
TABLES_DIR = ROOT / "results" / "tables"
STRUCT_DIR = ROOT / "results" / "structures"
OUT_DIR = ROOT / "figures" / "nature"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

PAL = {
    "blue": "#2B5C8C", "blue2": "#5B92C5", "blue3": "#94BFE5",
    "teal": "#3FA7A0", "purple": "#8F6AB8",
    "light": "#D8D8D8", "mid": "#8F8F8F", "dark": "#4D4D4D", "black": "#272727",
    "bg": "#EEF4FA",
}
ST_COLORS = {
    "intracellular.cupriavidus_like": "#2B5C8C",
    "intracellular.ralstonia_like": "#5B92C5",
    "intracellular.bacillus_like": "#94BFE5",
    "extracellular.general": "#3FA7A0",
    "extracellular.lemoignei_like": "#8F6AB8",
}
ST_ORDER = ["intracellular.cupriavidus_like","intracellular.ralstonia_like",
            "intracellular.bacillus_like","extracellular.general","extracellular.lemoignei_like"]
ST_LABELS = {
    "intracellular.cupriavidus_like": "Intra.\nCupriavidus",
    "intracellular.ralstonia_like": "Intra.\nRalstonia",
    "intracellular.bacillus_like": "Intra.\nBacillus",
    "extracellular.general": "Extra.\ngeneral",
    "extracellular.lemoignei_like": "Extra.\nlemoignei",
}


def apply_style(fs=7):
    matplotlib.rcParams.update({
        "font.size": fs, "axes.spines.right": False, "axes.spines.top": False,
        "axes.linewidth": 0.65, "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "legend.frameon": False,
        "pdf.fonttype": 42, "savefig.dpi": 600,
    })


def panel_label(ax, label, x=-0.08, y=1.04):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="left", va="bottom", color=PAL["black"])


def save(fig, name):
    base = OUT_DIR / name
    for ext in ("svg", "pdf", "png"):
        kw = {"bbox_inches": "tight"}
        if ext == "png": kw["dpi"] = 600
        fig.savefig(base.with_suffix(f".{ext}"), **kw)
    plt.close(fig)


# ===================== Fig 2.6: Copy number + PhaC co-occurrence =====================
def figure26():
    apply_style()
    fig = plt.figure(figsize=(7.2, 3.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Panel A: Copy number distribution
    panel_label(ax1, "a")
    cn = pd.read_csv(DATA_DIR / "phaz_validated_copy_number.tsv", sep="\t")
    copies, n_genomes = cn["copies"].values, cn["n_genomes"].values
    colors_a = [PAL["blue"] if c == 1 else PAL["mid"] for c in copies]
    bars = ax1.bar(copies, n_genomes, color=colors_a, edgecolor="white", width=0.7)
    ax1.set_xlabel("PhaZ copies per genome")
    ax1.set_ylabel("Genomes")
    ax1.set_title("Copy number distribution\n(validated, n=5,531)", loc="left", fontsize=8, pad=4)
    for c, n in zip(copies, n_genomes):
        ax1.text(c, n + 80, str(n), ha="center", fontsize=6, color=PAL["dark"])
    # Annotation
    ax1.text(0.95, 0.92, "84.7% single-copy", transform=ax1.transAxes, fontsize=7,
             ha="right", color=PAL["blue"], fontweight="bold")
    ax1.text(0.95, 0.84, f"Mean: 1.18 copies/genome", transform=ax1.transAxes, fontsize=6.5,
             ha="right", color=PAL["dark"])
    ax1.set_xticks(copies)
    ax1.grid(axis="y", color="#E5E5E5", linewidth=0.5)
    ax1.set_ylim(0, n_genomes.max() * 1.18)

    # Panel B: PhaC-phaZ co-occurrence
    panel_label(ax2, "b")
    phac = pd.read_csv(DATA_DIR / "phaz_phac_cross_reference.tsv", sep="\t")
    in_phaz = phac[phac["in_phaz_set"] == True]
    known_phac_with_phaz = phac[phac["known_phac"] == True]
    has_both = known_phac_with_phaz[known_phac_with_phaz["phaz_count"] > 0]
    has_both_sorted = has_both.nlargest(12, "phaz_count")

    y = np.arange(len(has_both_sorted))[::-1]
    ax2.barh(y, has_both_sorted["phaz_count"].values, color=PAL["teal"], height=0.65, edgecolor="white")
    ax2.set_yticks(y)
    ax2.set_yticklabels(has_both_sorted["genus"].values, fontsize=6.3, fontstyle="italic")
    ax2.set_xlabel("Validated PhaZ sequences")
    ax2.set_title("Genera with both phaZ & PhaC\n(12 known PhaC genera)", loc="left", fontsize=8, pad=4)
    for i, (_, row) in enumerate(has_both_sorted.iterrows()):
        ax2.text(row["phaz_count"] + 5, y[i], str(int(row["phaz_count"])), va="center", fontsize=5.8)
    ax2.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax2.tick_params(axis="y", length=0)

    # Annotation: genera without phaZ
    no_phaz = known_phac_with_phaz[known_phac_with_phaz["phaz_count"] == 0]
    ax2.text(0.98, 0.08,
             f"7 PhaC genera lack phaZ:\n{', '.join(no_phaz['genus'].values[:5])}...",
             transform=ax2.transAxes, fontsize=5.8, ha="right", va="bottom", color=PAL["dark"],
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.6))

    save(fig, "figure26_copy_number_phac")


# ===================== Fig 2.7: 3D structures + Conservation =====================
def figure27():
    apply_style()
    fig = plt.figure(figsize=(7.2, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.28)
    ax1 = fig.add_subplot(gs[0, :])   # PDB summary — spans full width top
    ax2 = fig.add_subplot(gs[1, 0])   # Conservation profile — Cupriavidus
    ax3 = fig.add_subplot(gs[1, 1])   # Conservation profile — Lemoignei

    # Panel A: PDB structure summary
    panel_label(ax1, "a")
    pdb_dir = STRUCT_DIR
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    # Count residues per file
    pdb_data = []
    for pf in pdb_files:
        ca = sum(1 for line in open(pf) if line.startswith("ATOM") and "CA  " in line)
        name = pf.stem
        # Map to subtype
        if name.startswith("intra_cupriavidus"): st = "intracellular.cupriavidus_like"
        elif name.startswith("intra_ralstonia"): st = "intracellular.ralstonia_like"
        elif name.startswith("intra_bacillus"): st = "intracellular.bacillus_like"
        elif name.startswith("extra_general"): st = "extracellular.general"
        elif name.startswith("extra_lemoignei"): st = "extracellular.lemoignei_like"
        elif "Cupriavidus" in name: st = "intracellular.cupriavidus_like"
        elif "Ralstonia" in name: st = "intracellular.ralstonia_like"
        elif "Bacillus" in name: st = "intracellular.bacillus_like"
        elif "General" in name: st = "extracellular.general"
        elif "Lemoignei" in name: st = "extracellular.lemoignei_like"
        else: st = "unknown"
        pdb_data.append({"subtype": st, "name": name, "ca_atoms": ca})
    pdb_df = pd.DataFrame(pdb_data)

    # Grouped bar by subtype
    st_list = [s for s in ST_ORDER if s in pdb_df["subtype"].values]
    positions = []
    labels_pos = []
    for i, st in enumerate(st_list):
        vals = pdb_df[pdb_df["subtype"] == st]["ca_atoms"].values
        x_pos = np.arange(i * 6, i * 6 + len(vals))
        ax1.bar(x_pos, vals, color=ST_COLORS.get(st, PAL["mid"]), edgecolor="white", width=0.8)
        # Mean line
        ax1.axhline(vals.mean(), xmin=(i*6-0.5)/30, xmax=(i*6+len(vals)-0.5)/30,
                    color=PAL["black"], linewidth=0.8, linestyle="--")
        positions.extend(x_pos)
        labels_pos.append(i * 6 + (len(vals) - 1) / 2)

    ax1.set_xticks(labels_pos)
    ax1.set_xticklabels([ST_LABELS.get(s, s) for s in st_list], fontsize=6.2, ha="center")
    ax1.set_ylabel("CA atoms (residues)")
    ax1.set_title("PhaZ 3D structures: resolved residues per PDB (ESMFold v1)", loc="left", fontsize=8, pad=4)
    ax1.text(0.99, 0.92, f"28 structures, 5 subtypes", transform=ax1.transAxes,
             fontsize=6.5, ha="right", color=PAL["dark"])
    ax1.grid(axis="y", color="#E5E5E5", linewidth=0.5)

    # Panel B: Conservation — Cupriavidus (intracellular)
    panel_label(ax2, "b")
    cons_file = TABLES_DIR / "conservation_intracellular_cupriavidus.tsv"
    if cons_file.exists():
        cons = pd.read_csv(cons_file, sep="\t")
        ax2.fill_between(cons["pos"], cons["consensus_pct"], alpha=0.3, color=PAL["blue"])
        ax2.plot(cons["pos"], cons["consensus_pct"], color=PAL["blue"], linewidth=0.8)
        ax2.axhline(50, color=PAL["mid"], linestyle="--", linewidth=0.6)
        ax2.set_xlabel("Alignment position")
        ax2.set_ylabel("Consensus (%)")
        ax2.set_title("Intra. Cupriavidus conservation\n(max 58%)", loc="left", fontsize=8, pad=4)
        ax2.set_ylim(0, 65)
        # Mark top conserved region
        top_idx = cons["consensus_pct"].idxmax()
        ax2.annotate(f"Pos {cons['pos'][top_idx]}: {cons['consensus_pct'][top_idx]:.0f}%",
                     (cons["pos"][top_idx], cons["consensus_pct"][top_idx]),
                     textcoords="offset points", xytext=(0, 10), fontsize=5.8, ha="center",
                     color=PAL["blue"])
        ax2.grid(color="#E5E5E5", linewidth=0.5)
    else:
        ax2.text(0.5, 0.5, "Conservation data not found", ha="center", va="center", fontsize=8)
        ax2.set_title("Conservation (data missing)", loc="left", fontsize=8, pad=4)

    # Panel C: Conservation — Lemoignei (extracellular, highest conservation)
    panel_label(ax3, "c")
    cons_lemo_file = TABLES_DIR / "conservation_extracellular_lemoignei.tsv"
    if cons_lemo_file.exists():
        cons_l = pd.read_csv(cons_lemo_file, sep="\t")
        ax3.fill_between(cons_l["pos"], cons_l["consensus_pct"], alpha=0.3, color=PAL["purple"])
        ax3.plot(cons_l["pos"], cons_l["consensus_pct"], color=PAL["purple"], linewidth=0.8)
        ax3.axhline(50, color=PAL["mid"], linestyle="--", linewidth=0.6)
        ax3.set_xlabel("Alignment position")
        ax3.set_ylabel("Consensus (%)")
        ax3.set_title("Extra. Lemoignei conservation\n(max 79%)", loc="left", fontsize=8, pad=4)
        ax3.set_ylim(0, 85)
        top_idx_l = cons_l["consensus_pct"].idxmax()
        ax3.annotate(f"G at pos {cons_l['pos'][top_idx_l]}: {cons_l['consensus_pct'][top_idx_l]:.0f}%",
                     (cons_l["pos"][top_idx_l], cons_l["consensus_pct"][top_idx_l]),
                     textcoords="offset points", xytext=(0, 10), fontsize=5.8, ha="center",
                     color=PAL["purple"])
        ax3.grid(color="#E5E5E5", linewidth=0.5)
    else:
        ax3.text(0.5, 0.5, "Conservation data not found", ha="center", va="center", fontsize=8)
        ax3.set_title("Conservation (data missing)", loc="left", fontsize=8, pad=4)

    save(fig, "figure27_structures_conservation")


# ===================== Fig 2.8: Signal peptide + Seed sequences =====================
def figure28():
    apply_style()
    fig = plt.figure(figsize=(7.2, 3.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Panel A: Signal peptide
    panel_label(ax1, "a")
    sp = pd.read_csv(DATA_DIR / "phaz_signal_peptide_summary.tsv", sep="\t")
    # Map to ST_ORDER
    sp_map = {"intracellular": "intracellular.cupriavidus_like", "ralstonia": "intracellular.ralstonia_like",
              "bacillus_type": "intracellular.bacillus_like", "extracellular": "extracellular.general",
              "extracellular_lemoignei": "extracellular.lemoignei_like"}
    sp["subtype_clean"] = sp["subtype"].map(sp_map)
    sp_ordered = sp.set_index("subtype_clean").reindex(ST_ORDER).dropna()

    y = np.arange(len(sp_ordered))[::-1]
    colors_sp = [ST_COLORS.get(s, PAL["mid"]) for s in sp_ordered.index]
    ax1.barh(y, sp_ordered["pct"].values, color=colors_sp, height=0.65, edgecolor="white")
    ax1.set_yticks(y)
    ax1.set_yticklabels([ST_LABELS.get(s, s) for s in sp_ordered.index], fontsize=6.3)
    ax1.set_xlabel("Sequences with putative signal peptide (%)")
    ax1.set_title("Signal peptide prediction\n(Kyte-Doolittle hydrophobicity)", loc="left", fontsize=8, pad=4)
    for i, (_, row) in enumerate(sp_ordered.iterrows()):
        ax1.text(row["pct"] + 0.5, i, f"{row['pct']:.1f}%", va="center", fontsize=6.2)
    ax1.set_xlim(0, 20)
    # Add extra vs intra summary
    intra_pct = sp_ordered.iloc[:3]["signal_positive"].sum() / sp_ordered.iloc[:3]["total"].sum() * 100
    extra_pct = sp_ordered.iloc[3:]["signal_positive"].sum() / sp_ordered.iloc[3:]["total"].sum() * 100
    ax1.text(0.98, 0.92, f"Extra: {extra_pct:.1f}%  |  Intra: {intra_pct:.1f}%\nFold enrichment: {extra_pct/intra_pct:.1f}x",
             transform=ax1.transAxes, fontsize=6.2, ha="right", va="top", color=PAL["dark"],
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.6))
    ax1.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax1.tick_params(axis="y", length=0)

    # Panel B: Seed sequences
    panel_label(ax2, "b")
    seeds = pd.read_csv(DATA_DIR / "phaz_metagenome_seeds.tsv", sep="\t")
    # Group by subtype
    seed_counts = seeds["subtype"].value_counts()
    seed_ordered = [seed_counts.get(s, 0) for s in ST_ORDER]
    y2 = np.arange(len(ST_ORDER))[::-1]
    colors_sd = [ST_COLORS.get(s, PAL["mid"]) for s in ST_ORDER]
    ax2.barh(y2, seed_ordered, color=colors_sd, height=0.65, edgecolor="white")
    ax2.set_yticks(y2)
    ax2.set_yticklabels([ST_LABELS.get(s, s) for s in ST_ORDER], fontsize=6.3)
    ax2.set_xlabel("Representative seed sequences")
    ax2.set_title("Metagenome screening seeds\n(CD-HIT c70 + manual curation)", loc="left", fontsize=8, pad=4)
    for i, val in enumerate(seed_ordered):
        ax2.text(val + 0.2, i, str(val), va="center", fontsize=6.2)
    ax2.text(0.98, 0.92, f"Total: {len(seeds)} seeds\nFor Chapter 4\nmetagenome screening",
             transform=ax2.transAxes, fontsize=6.2, ha="right", va="top", color=PAL["dark"],
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.6))
    ax2.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax2.tick_params(axis="y", length=0)
    ax2.set_xlim(0, max(seed_ordered) * 1.4)

    save(fig, "figure28_signal_seeds")


def main():
    print("Generating Fig 2.6 ...")
    figure26()
    print("Generating Fig 2.7 ...")
    figure27()
    print("Generating Fig 2.8 ...")
    figure28()
    print(f"Done. Output: {OUT_DIR}/")


if __name__ == "__main__":
    main()
