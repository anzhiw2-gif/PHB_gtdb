#!/usr/bin/env python3
"""Generate Nature-style summary figures for the PHB_gtdb project.

Inputs are copied from the T141 runtime into ./figure_data. The script writes
editable SVG/PDF plus high-resolution PNG outputs under ./figures/nature.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib import gridspec, patches


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "figure_data"
OUT_DIR = ROOT / "figures" / "nature"
SOURCE_DIR = OUT_DIR / "source_data"


# Mandatory editable SVG text settings from nature-figure.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"


PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "red_strong": "#B64342",
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#8F8F8F",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "bg_blue": "#EEF4FA",
    "bg_teal": "#EAF6F5",
    "bg_red": "#F9ECEA",
}

SUBTYPE_ORDER = [
    "intracellular.cupriavidus_like",
    "intracellular.ralstonia_like",
    "intracellular.bacillus_like",
    "extracellular.general",
    "extracellular.lemoignei_like",
]

SUBTYPE_LABELS = {
    "intracellular.cupriavidus_like": "Intra.\n(Cupriavidus)",
    "intracellular.ralstonia_like": "Intra.\n(Ralstonia)",
    "intracellular.bacillus_like": "Intra.\n(Bacillus)",
    "extracellular.general": "Extra.\n(general)",
    "extracellular.lemoignei_like": "Extra.\n(lemoignei)",
}

SUBTYPE_COLORS = {
    "intracellular.cupriavidus_like": "#2B5C8C",
    "intracellular.ralstonia_like": "#5B92C5",
    "intracellular.bacillus_like": "#94BFE5",
    "extracellular.general": "#3FA7A0",
    "extracellular.lemoignei_like": "#8F6AB8",
}

REF_TO_SUBTYPE = {
    "BAA33394.1": "intracellular.cupriavidus_like",
    "CAJ93939.1": "intracellular.cupriavidus_like",
    "CAJ95805.1": "intracellular.cupriavidus_like",
    "UCA14981.1": "intracellular.ralstonia_like",
    "WKZ88401.1": "intracellular.ralstonia_like",
    "BAA19791.1": "extracellular.general",
    "AAA87070.1": "extracellular.general",
    "BAA35137.1": "extracellular.general",
    "BAA32541.1": "extracellular.general",
    "AAB02914.1": "extracellular.general",
    "P52090.1": "extracellular.lemoignei_like",
    "WP_207907290.1": "extracellular.lemoignei_like",
    "WP_243656647.1": "extracellular.lemoignei_like",
    "WP_128854079.1": "intracellular.bacillus_like",
}

LIPASE_BOX_PERCENT = {
    "intracellular.cupriavidus_like": 10.1,
    "intracellular.ralstonia_like": 9.8,
    "intracellular.bacillus_like": 50.0,
    "extracellular.general": 28.6,
    "extracellular.lemoignei_like": 53.2,
}


@dataclass(frozen=True)
class PhaZRecord:
    genome_id: str
    phylum: str
    genus: str
    species: str
    protein_id: str
    ref: str
    subtype: str


def apply_style(font_size: float = 7.0) -> None:
    mpl.rcParams.update(
        {
            "font.size": font_size,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def add_panel_label(ax, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=PALETTE["neutral_black"],
    )


def save_figure(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base = OUT_DIR / name
    for ext in ("svg", "pdf", "png"):
        kwargs = {"bbox_inches": "tight"}
        if ext == "png":
            kwargs["dpi"] = 600
        fig.savefig(base.with_suffix(f".{ext}"), **kwargs)
    plt.close(fig)


def normalize_ref(raw_ref: str) -> str:
    for accession in REF_TO_SUBTYPE:
        if accession in raw_ref:
            return accession
    return raw_ref


def parse_validated_fasta(path: Path) -> pd.DataFrame:
    records: list[PhaZRecord] = []
    header_re = re.compile(r"ref=([^\s]+)")
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            meta = line[1:].strip()
            token = meta.split()[0]
            parts = token.split("|")
            if len(parts) < 5:
                raise ValueError(f"Unexpected FASTA header: {line[:120]}")
            match = header_re.search(meta)
            if not match:
                raise ValueError(f"Missing ref= field: {line[:120]}")
            ref = normalize_ref(match.group(1))
            subtype = REF_TO_SUBTYPE.get(ref)
            if subtype is None:
                raise ValueError(f"Unmapped reference accession: {ref}")
            records.append(
                PhaZRecord(
                    genome_id=parts[0],
                    phylum=parts[1],
                    genus=parts[2],
                    species=parts[3],
                    protein_id=parts[4],
                    ref=ref,
                    subtype=subtype,
                )
            )
    df = pd.DataFrame([r.__dict__ for r in records])
    if len(df) != 6532:
        raise ValueError(f"Expected 6532 validated records, observed {len(df)}")
    return df


def prepare_source_data(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)

    subtype_counts = (
        df["subtype"]
        .value_counts()
        .reindex(SUBTYPE_ORDER)
        .rename_axis("subtype")
        .reset_index(name="count")
    )
    subtype_counts["percent"] = subtype_counts["count"] / subtype_counts["count"].sum() * 100
    subtype_counts["lipase_box_percent"] = subtype_counts["subtype"].map(LIPASE_BOX_PERCENT)

    phylum_subtype = pd.crosstab(df["phylum"], df["subtype"]).reindex(columns=SUBTYPE_ORDER, fill_value=0)
    phylum_subtype["total"] = phylum_subtype.sum(axis=1)
    phylum_subtype = phylum_subtype.sort_values("total", ascending=False)

    genus_subtype = pd.crosstab(df["genus"], df["subtype"]).reindex(columns=SUBTYPE_ORDER, fill_value=0)
    genus_subtype["total"] = genus_subtype.sum(axis=1)
    genus_subtype = genus_subtype.sort_values("total", ascending=False)

    tree_summary = pd.DataFrame(
        [
            {
                "subtype": "intracellular.bacillus_like",
                "tips": 57,
                "method": "IQ-TREE",
                "treefile": "phaz_intracellular_bacillus_like_tree.treefile",
            },
            {
                "subtype": "extracellular.general",
                "tips": 509,
                "method": "IQ-TREE",
                "treefile": "phaz_extracellular_general_tree.treefile",
            },
            {
                "subtype": "extracellular.lemoignei_like",
                "tips": 412,
                "method": "IQ-TREE",
                "treefile": "phaz_extracellular_lemoignei_like_tree.treefile",
            },
            {
                "subtype": "intracellular.ralstonia_like",
                "tips": 2776,
                "method": "FastTree",
                "treefile": "phaz_intracellular_ralstonia_like_tree_ft.treefile",
            },
            {
                "subtype": "intracellular.cupriavidus_like",
                "tips": 4424,
                "method": "FastTree",
                "treefile": "phaz_intracellular_cupriavidus_like_tree_ft.treefile",
            },
        ]
    )

    funnel = pd.DataFrame(
        [
            ("Bacterial DIAMOND hits", 16486),
            ("Extracted candidates", 8768),
            ("CD-HIT c95", 8731),
            (">=100 aa", 7478),
            ("HMMer verified", 6033),
            ("Final validated union", 6532),
        ],
        columns=["stage", "count"],
    )

    phylum_subtype.to_csv(SOURCE_DIR / "figure2_phylum_subtype_matrix.csv")
    subtype_counts.to_csv(SOURCE_DIR / "figure3_subtype_lipase.csv", index=False)
    genus_subtype.head(20).to_csv(SOURCE_DIR / "figure4_top_genera_subtype.csv")
    tree_summary.to_csv(SOURCE_DIR / "figure4_tree_summary.csv", index=False)
    funnel.to_csv(SOURCE_DIR / "figure1_funnel.csv", index=False)

    return {
        "subtype_counts": subtype_counts,
        "phylum_subtype": phylum_subtype,
        "genus_subtype": genus_subtype,
        "tree_summary": tree_summary,
        "funnel": funnel,
    }


def draw_arrow(ax, x1, y1, x2, y2, color=PALETTE["neutral_mid"]) -> None:
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=color, shrinkA=4, shrinkB=4),
    )


def figure1(stats: dict[str, pd.DataFrame]) -> None:
    apply_style()
    fig = plt.figure(figsize=(7.2, 4.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.26)
    ax_flow = fig.add_subplot(gs[0, 0])
    ax_fun = fig.add_subplot(gs[0, 1])

    add_panel_label(ax_flow, "a")
    ax_flow.set_axis_off()
    ax_flow.set_title("Genome-to-PhaZ workflow", loc="left", fontsize=8, pad=4)

    steps = [
        ("1", "GTDB R232", "189,801 bacteria\n10,122 archaea"),
        ("2", "Protein calling", "Pyrodigal\ntranslated ORFs"),
        ("3", "PhaZ search", "14 curated refs\nDIAMOND blastp"),
        ("4", "Candidate QC", "CD-HIT c95\nlength >=100 aa"),
        ("5", "Validation", "HMMer + DIAMOND\nhigh-conf. union"),
        ("6", "Subtype trees", "MAFFT / trimAl\nIQ-TREE / FastTree"),
    ]
    xs = [0.16, 0.50, 0.84, 0.84, 0.50, 0.16]
    ys = [0.74, 0.74, 0.74, 0.34, 0.34, 0.34]
    box_w, box_h = 0.27, 0.21
    box_colors = [
        PALETTE["bg_blue"],
        "#F7F7F7",
        PALETTE["bg_teal"],
        "#F7F7F7",
        PALETTE["bg_red"],
        "#F7F7F7",
    ]
    for i, ((number, title, subtitle), x, y, color) in enumerate(zip(steps, xs, ys, box_colors)):
        rect = patches.FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2),
            box_w,
            box_h,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            linewidth=0.8,
            edgecolor=PALETTE["neutral_dark"],
            facecolor=color,
        )
        ax_flow.add_patch(rect)
        badge_x = x - box_w / 2 + 0.034
        badge_y = y + box_h / 2 - 0.034
        ax_flow.add_patch(
            patches.Circle(
                (badge_x, badge_y),
                0.020,
                facecolor="white",
                edgecolor=PALETTE["neutral_mid"],
                linewidth=0.65,
                zorder=3,
            )
        )
        ax_flow.text(
            badge_x,
            badge_y,
            number,
            ha="center",
            va="center",
            fontsize=5.6,
            fontweight="bold",
            color=PALETTE["neutral_dark"],
            zorder=4,
        )
        ax_flow.text(x, y + 0.030, title, ha="center", va="center", fontsize=6.6, fontweight="bold")
        ax_flow.text(
            x,
            y - 0.044,
            subtitle,
            ha="center",
            va="center",
            fontsize=5.7,
            color=PALETTE["neutral_dark"],
            linespacing=1.15,
        )

    connectors = [
        (xs[0] + box_w / 2, ys[0], xs[1] - box_w / 2, ys[1]),
        (xs[1] + box_w / 2, ys[1], xs[2] - box_w / 2, ys[2]),
        (xs[2], ys[2] - box_h / 2, xs[3], ys[3] + box_h / 2),
        (xs[3] - box_w / 2, ys[3], xs[4] + box_w / 2, ys[4]),
        (xs[4] - box_w / 2, ys[4], xs[5] + box_w / 2, ys[5]),
    ]
    for x1, y1, x2, y2 in connectors:
        draw_arrow(ax_flow, x1, y1, x2, y2)
    ax_flow.text(
        0.03,
        0.08,
        "Archaeal screen: 2 short DIAMOND hits -> 0 validated PhaZ",
        ha="left",
        va="center",
        fontsize=6.2,
        color=PALETTE["red_strong"],
    )
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)

    add_panel_label(ax_fun, "b")
    funnel = stats["funnel"].copy()
    y = np.arange(len(funnel))[::-1]
    widths = funnel["count"].to_numpy()
    max_width = widths.max()
    colors = [PALETTE["blue_secondary"]] * (len(funnel) - 1) + [PALETTE["blue_main"]]
    ax_fun.barh(y, widths, color=colors, edgecolor="white", height=0.68)
    ax_fun.set_yticks(y)
    ax_fun.set_yticklabels(funnel["stage"], fontsize=6.4)
    ax_fun.set_xlabel("PhaZ candidates / proteins", fontsize=7)
    ax_fun.set_xlim(0, max_width * 1.18)
    ax_fun.set_title("Screening funnel", loc="left", fontsize=8, pad=4)
    for yi, count in zip(y, widths):
        ax_fun.text(count + max_width * 0.025, yi, f"{count:,}", va="center", fontsize=6.5)
    ax_fun.text(
        max_width * 0.68,
        0.55,
        "+499 high-confidence\nDIAMOND rescue",
        ha="left",
        va="center",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.6),
    )
    ax_fun.spines["left"].set_visible(False)
    ax_fun.tick_params(axis="y", length=0)
    ax_fun.grid(axis="x", color="#E5E5E5", linewidth=0.5)

    save_figure(fig, "figure1_workflow_funnel")


def figure2(stats: dict[str, pd.DataFrame]) -> None:
    apply_style()
    phylum_subtype = stats["phylum_subtype"]
    top_rows = phylum_subtype.head(13).copy()
    matrix = top_rows[SUBTYPE_ORDER].to_numpy()
    totals = top_rows["total"]

    fig = plt.figure(figsize=(7.2, 4.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.38], wspace=0.18)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])

    add_panel_label(ax_bar, "a")
    y = np.arange(len(top_rows))[::-1]
    phyla = list(top_rows.index)
    bar_colors = [PALETTE["blue_main"] if p == "Pseudomonadota" else PALETTE["neutral_mid"] for p in phyla]
    ax_bar.barh(y, totals.values, color=bar_colors, height=0.68)
    ax_bar.set_xscale("log")
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(phyla, fontsize=6.3)
    ax_bar.set_xlabel("Validated PhaZ proteins (log scale)")
    ax_bar.set_title("Phylum-level distribution", loc="left", fontsize=8, pad=4)
    for yi, value in zip(y, totals.values):
        ax_bar.text(value * 1.12, yi, f"{value:,}", va="center", fontsize=6)
    ax_bar.grid(axis="x", which="major", color="#E5E5E5", linewidth=0.5)
    ax_bar.tick_params(axis="y", length=0)

    add_panel_label(ax_heat, "b")
    log_matrix = np.log10(matrix + 1)
    im = ax_heat.imshow(log_matrix, cmap="Blues", aspect="auto", vmin=0)
    ax_heat.set_xticks(np.arange(len(SUBTYPE_ORDER)))
    ax_heat.set_xticklabels([SUBTYPE_LABELS[s] for s in SUBTYPE_ORDER], rotation=30, ha="right", fontsize=6.3)
    ax_heat.set_yticks(np.arange(len(phyla)))
    ax_heat.set_yticklabels(phyla, fontsize=6.3)
    ax_heat.set_title("Phylum x subtype matrix", loc="left", fontsize=8, pad=4)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix[i, j])
            if value == 0:
                continue
            text_color = "white" if log_matrix[i, j] > log_matrix.max() * 0.55 else PALETTE["neutral_black"]
            ax_heat.text(j, i, str(value), ha="center", va="center", fontsize=5.5, color=text_color)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.tick_params(length=0)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.035, pad=0.02)
    cbar.set_label("log10(count + 1)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)

    save_figure(fig, "figure2_phylum_heatmap")


def figure3(stats: dict[str, pd.DataFrame]) -> None:
    apply_style()
    subtype_counts = stats["subtype_counts"].copy()

    fig = plt.figure(figsize=(7.2, 3.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.28)
    ax_comp = fig.add_subplot(gs[0, 0])
    ax_lip = fig.add_subplot(gs[0, 1])

    add_panel_label(ax_comp, "a")
    y = np.arange(len(subtype_counts))[::-1]
    colors = [SUBTYPE_COLORS[s] for s in subtype_counts["subtype"]]
    ax_comp.barh(y, subtype_counts["percent"], color=colors, height=0.68, edgecolor="white")
    ax_comp.set_yticks(y)
    ax_comp.set_yticklabels([SUBTYPE_LABELS[s] for s in subtype_counts["subtype"]], fontsize=6.8)
    ax_comp.set_xlabel("Share of validated PhaZ (%)")
    ax_comp.set_xlim(0, 62)
    ax_comp.set_title("Subtype composition (n=6,532)", loc="left", fontsize=8, pad=4)
    for yi, row in zip(y, subtype_counts.itertuples()):
        ax_comp.text(
            row.percent + 1.1,
            yi,
            f"{int(row.count):,} ({row.percent:.1f}%)",
            va="center",
            fontsize=6.4,
        )
    ax_comp.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax_comp.tick_params(axis="y", length=0)

    add_panel_label(ax_lip, "b")
    lip_values = subtype_counts["lipase_box_percent"].to_numpy()
    ax_lip.barh(y, lip_values, color=colors, height=0.68, edgecolor="white")
    ax_lip.set_yticks(y)
    ax_lip.set_yticklabels([SUBTYPE_LABELS[s] for s in subtype_counts["subtype"]], fontsize=6.8)
    ax_lip.set_xlabel("Sequences with lipase box (%)")
    ax_lip.set_xlim(0, 60)
    ax_lip.set_title("Catalytic motif support", loc="left", fontsize=8, pad=4)
    ax_lip.axvline(10, color=PALETTE["neutral_mid"], linestyle="--", linewidth=0.7)
    ax_lip.text(10.8, -0.48, "intracellular baseline\n(~10%)", fontsize=5.8, color=PALETTE["neutral_dark"])
    for yi, value in zip(y, lip_values):
        ax_lip.text(value + 1.0, yi, f"{value:.1f}%", va="center", fontsize=6.4)
    ax_lip.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax_lip.tick_params(axis="y", length=0)

    save_figure(fig, "figure3_subtype_lipase")


def draw_tree_silhouette(ax, tree_path: Path, color: str) -> None:
    tree = Phylo.read(str(tree_path), "newick")
    try:
        tree.ladderize()
    except Exception:
        pass
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda _: None)
    for line in ax.lines:
        line.set_color(color)
        line.set_linewidth(0.25)
        line.set_alpha(0.85)
    for collection in ax.collections:
        collection.set_color(color)
        collection.set_linewidth(0.25)
        collection.set_alpha(0.85)
    ax.set_axis_off()


def figure4(stats: dict[str, pd.DataFrame]) -> None:
    apply_style()
    genus_subtype = stats["genus_subtype"].head(15).copy()
    tree_summary = stats["tree_summary"].copy()

    fig = plt.figure(figsize=(7.2, 5.7), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.08, 0.92], wspace=0.24)
    ax_genus = fig.add_subplot(gs[0, 0])
    right = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[0, 1], hspace=0.08)
    tree_axes = [fig.add_subplot(right[i, 0]) for i in range(5)]

    add_panel_label(ax_genus, "a")
    genera = list(genus_subtype.index)
    y = np.arange(len(genera))[::-1]
    left = np.zeros(len(genera))
    for subtype in SUBTYPE_ORDER:
        values = genus_subtype[subtype].to_numpy()
        ax_genus.barh(
            y,
            values,
            left=left,
            color=SUBTYPE_COLORS[subtype],
            edgecolor="white",
            height=0.68,
            label=SUBTYPE_LABELS[subtype],
        )
        left += values
    ax_genus.set_yticks(y)
    ax_genus.set_yticklabels([f"{g}" for g in genera], fontsize=6.3, fontstyle="italic")
    ax_genus.set_xlabel("Validated PhaZ proteins")
    ax_genus.set_title("Top genera by subtype", loc="left", fontsize=8, pad=4)
    for yi, total in zip(y, genus_subtype["total"]):
        ax_genus.text(total + 4, yi, f"{int(total)}", va="center", fontsize=5.8)
    ax_genus.set_xlim(0, genus_subtype["total"].max() * 1.25)
    ax_genus.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax_genus.tick_params(axis="y", length=0)
    ax_genus.legend(
        ncol=2,
        loc="lower right",
        bbox_to_anchor=(1.0, -0.18),
        fontsize=5.7,
        handlelength=1.0,
        columnspacing=0.8,
    )

    ordered_trees = tree_summary.set_index("subtype").loc[SUBTYPE_ORDER]
    for idx, (subtype, row) in enumerate(ordered_trees.iterrows()):
        ax = tree_axes[idx]
        if idx == 0:
            add_panel_label(ax, "b", x=-0.10, y=1.08)
            ax.set_title("Subtype tree overview", loc="left", fontsize=8, pad=2)
        draw_tree_silhouette(ax, DATA_DIR / row["treefile"], SUBTYPE_COLORS[subtype])
        ax.text(
            0.01,
            0.92,
            f"{SUBTYPE_LABELS[subtype]}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.4,
            fontweight="bold",
            color=SUBTYPE_COLORS[subtype],
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.8),
        )
        ax.text(
            0.99,
            0.08,
            f"{int(row['tips']):,} tips, {row['method']}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.6,
            color=PALETTE["neutral_dark"],
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.8),
        )

    save_figure(fig, "figure4_genera_phylogeny")


def main() -> None:
    df = parse_validated_fasta(DATA_DIR / "phaz_proteins_validated.fasta")
    stats = prepare_source_data(df)
    figure1(stats)
    figure2(stats)
    figure3(stats)
    figure4(stats)
    print(f"Generated Nature-style figures in {OUT_DIR}")


if __name__ == "__main__":
    main()
