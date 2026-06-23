#!/usr/bin/env python3
"""
Statistical tests and publication-style figure for gRodon2 growth results.

Primary analysis unit: genus. Within each genus, phaZ-positive and
phaZ-negative genomes are first balanced to the same count, then the genus
mean growth-rate difference is tested against zero.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".mplconfig"))

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    }
)


PALETTE = {
    "positive": "#2C7FB8",
    "negative": "#D95F02",
    "neutral": "#4D4D4D",
    "light": "#E6E6E6",
    "very_light": "#F5F5F5",
    "accent": "#2E8B57",
    "warn": "#B2182B",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run genus-level gRodon2 statistics and draw Figure 5."
    )
    parser.add_argument(
        "--balanced",
        type=Path,
        default=PROJECT_ROOT / "results/tables/grodon_growth_balanced_by_genus_hmm_allmatched.tsv",
        help="Balanced genome-level gRodon2 table.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=PROJECT_ROOT
        / "results/tables/grodon_growth_balanced_by_genus_summary_hmm_allmatched.tsv",
        help="Balanced genus-level summary table.",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=PROJECT_ROOT
        / "results/tables/grodon_growth_balanced_genus_selected_counts_hmm_allmatched.tsv",
        help="Per-genus selected count table.",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=PROJECT_ROOT / "results/tables/grodon_growth_balance_audit_hmm_allmatched.tsv",
        help="Balance audit table.",
    )
    parser.add_argument(
        "--failed",
        type=Path,
        default=PROJECT_ROOT / "results/tables/grodon_failed_genomes_hmm_allmatched.tsv",
        help="Failed gRodon2 genome table.",
    )
    parser.add_argument("--label", default="hmm_allmatched")
    parser.add_argument(
        "--total-genomes",
        type=int,
        default=8788,
        help="Total genomes submitted to the gRodon2 run.",
    )
    parser.add_argument("--seed", type=int, default=20260623)
    parser.add_argument("--n-bootstrap", type=int, default=10000)
    parser.add_argument("--n-permutations", type=int, default=5000)
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=PROJECT_ROOT / "results/tables",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=PROJECT_ROOT / "figures/nature",
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=PROJECT_ROOT / "figures/nature/source_data",
    )
    return parser.parse_args()


def normal_p_two_sided(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def exact_sign_test_p(k_positive: int, n_nonzero: int) -> float:
    k = min(k_positive, n_nonzero - k_positive)
    prob = sum(math.comb(n_nonzero, i) for i in range(k + 1)) / (2**n_nonzero)
    return min(1.0, 2.0 * prob)


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and sorted_values[j] == sorted_values[i]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def wilcoxon_signed_rank_normal(delta: np.ndarray) -> dict[str, float]:
    nonzero = delta[delta != 0]
    n = len(nonzero)
    abs_delta = np.abs(nonzero)
    ranks = average_ranks(abs_delta)
    w_pos = float(ranks[nonzero > 0].sum())
    w_neg = float(ranks[nonzero < 0].sum())
    expected = n * (n + 1) / 4.0

    _, tie_counts = np.unique(abs_delta, return_counts=True)
    tie_term = float(np.sum(tie_counts**3 - tie_counts))
    variance = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term / 48.0

    z = (w_pos - expected) / math.sqrt(variance)
    p = normal_p_two_sided(z)
    return {
        "n_nonzero": n,
        "w_positive": w_pos,
        "w_negative": w_neg,
        "z": z,
        "p_value": p,
        "effect_r": z / math.sqrt(n),
    }


def bootstrap_ci(
    values: np.ndarray,
    reducer,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> tuple[float, float]:
    values = np.asarray(values)
    estimates = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        sample = values[rng.integers(0, len(values), size=len(values))]
        estimates[i] = reducer(sample)
    return tuple(np.percentile(estimates, [2.5, 97.5]))


def stratified_permutation(
    balanced: pd.DataFrame,
    rng: np.random.Generator,
    n_permutations: int,
) -> tuple[float, float, float, float]:
    groups = []
    observed_unweighted = []
    observed_weighted_num = 0.0
    observed_weighted_den = 0
    for _, group in balanced.groupby("genus", sort=True):
        pos_values = group.loc[
            group["phaZ_status"] == "phaZ_positive", "growth_rate_per_h"
        ].to_numpy(float)
        neg_values = group.loc[
            group["phaZ_status"] == "phaZ_negative", "growth_rate_per_h"
        ].to_numpy(float)
        n_pos = len(pos_values)
        n_neg = len(neg_values)
        if n_pos != n_neg:
            raise ValueError("Balanced table contains unequal within-genus counts.")
        delta = float(pos_values.mean() - neg_values.mean())
        observed_unweighted.append(delta)
        observed_weighted_num += delta * n_pos
        observed_weighted_den += n_pos
        values = np.concatenate([pos_values, neg_values])
        groups.append((values, n_pos))

    observed_unweighted = float(np.mean(observed_unweighted))
    observed_weighted = float(observed_weighted_num / observed_weighted_den)

    null_unweighted = np.empty(n_permutations, dtype=float)
    null_weighted = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        deltas = []
        weighted_num = 0.0
        weighted_den = 0
        for values, n_pos in groups:
            perm = rng.permutation(values)
            delta = float(perm[:n_pos].mean() - perm[n_pos:].mean())
            deltas.append(delta)
            weighted_num += delta * n_pos
            weighted_den += n_pos
        null_unweighted[i] = np.mean(deltas)
        null_weighted[i] = weighted_num / weighted_den

    p_unweighted = (np.sum(np.abs(null_unweighted) >= abs(observed_unweighted)) + 1) / (
        n_permutations + 1
    )
    p_weighted = (np.sum(np.abs(null_weighted) >= abs(observed_weighted)) + 1) / (
        n_permutations + 1
    )
    return observed_unweighted, p_unweighted, observed_weighted, p_weighted


def prepare_tables(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    balanced = pd.read_csv(args.balanced, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    counts = pd.read_csv(args.counts, sep="\t")

    balanced = balanced.copy()
    balanced["growth_rate_per_h"] = pd.to_numeric(
        balanced["growth_rate_per_h"], errors="coerce"
    )
    balanced = balanced[balanced["growth_rate_per_h"].notna()].copy()

    effects = summary.merge(
        counts[
            [
                "genus",
                "selected_phaZ_positive",
                "selected_phaZ_negative",
                "total_selected",
                "ok_phaZ_positive",
                "ok_phaZ_negative",
            ]
        ],
        on="genus",
        how="left",
    )
    effects["relative_delta_percent"] = (
        effects["delta_positive_minus_negative"]
        / effects["mean_growth_phaZ_negative"].replace(0, np.nan)
        * 100.0
    )
    effects["direction"] = np.where(
        effects["delta_positive_minus_negative"] > 0,
        "phaZ_positive_higher",
        "phaZ_negative_higher",
    )
    effects = effects.sort_values(
        ["n_phaZ_positive", "delta_positive_minus_negative"],
        ascending=[False, False],
    )
    return balanced, summary, effects


def run_statistics(
    args: argparse.Namespace,
    balanced: pd.DataFrame,
    summary: pd.DataFrame,
    effects: pd.DataFrame,
) -> pd.DataFrame:
    rng = np.random.default_rng(args.seed)
    delta = summary["delta_positive_minus_negative"].to_numpy(float)
    weights = summary["n_phaZ_positive"].to_numpy(float)

    mean_delta = float(np.mean(delta))
    median_delta = float(np.median(delta))
    weighted_mean = float(np.average(delta, weights=weights))
    mean_ci = bootstrap_ci(delta, np.mean, rng, args.n_bootstrap)
    median_ci = bootstrap_ci(delta, np.median, rng, args.n_bootstrap)
    weighted_ci = bootstrap_ci(
        np.arange(len(delta)),
        lambda idx: np.average(delta[idx.astype(int)], weights=weights[idx.astype(int)]),
        rng,
        args.n_bootstrap,
    )

    positive = int(np.sum(delta > 0))
    negative = int(np.sum(delta < 0))
    sign_p = exact_sign_test_p(positive, positive + negative)
    wilcoxon = wilcoxon_signed_rank_normal(delta)
    perm_unweighted, perm_p_unweighted, perm_weighted, perm_p_weighted = stratified_permutation(
        balanced, rng, args.n_permutations
    )

    failed_count = 0
    failed_genera = 0
    if args.failed.exists():
        failed = pd.read_csv(args.failed, sep="\t")
        failed_count = len(failed)
        failed_genera = failed["genus"].nunique()

    audit_excluded = 0
    if args.audit.exists():
        audit = pd.read_csv(args.audit, sep="\t")
        audit_excluded = int(
            audit["n_excluded_positive"].sum() + audit["n_excluded_negative"].sum()
        )
    missing_growth = args.total_genomes - failed_count - len(balanced) - audit_excluded

    rows = [
        ("n_grodon2_input_genomes", args.total_genomes, "", ""),
        ("n_balanced_genomes", len(balanced), "", ""),
        ("n_balanced_phaZ_positive", int((balanced["phaZ_status"] == "phaZ_positive").sum()), "", ""),
        ("n_balanced_phaZ_negative", int((balanced["phaZ_status"] == "phaZ_negative").sum()), "", ""),
        ("n_balanced_genera", summary["genus"].nunique(), "", ""),
        ("n_failed_genomes", failed_count, "", ""),
        ("n_failed_genera", failed_genera, "", ""),
        ("n_status_ok_but_missing_growth_rate", missing_growth, "", ""),
        ("n_successful_excluded_by_balance", audit_excluded, "", ""),
        ("genus_mean_delta_growth_rate_per_h", mean_delta, f"{mean_ci[0]},{mean_ci[1]}", "bootstrap_95_ci"),
        ("genus_median_delta_growth_rate_per_h", median_delta, f"{median_ci[0]},{median_ci[1]}", "bootstrap_95_ci"),
        (
            "genus_weighted_mean_delta_growth_rate_per_h",
            weighted_mean,
            f"{weighted_ci[0]},{weighted_ci[1]}",
            "bootstrap_95_ci_weighted_by_selected_pairs",
        ),
        ("n_genera_delta_positive", positive, "", ""),
        ("n_genera_delta_negative", negative, "", ""),
        ("exact_sign_test_two_sided_p", sign_p, "", "genus_delta_signs"),
        ("wilcoxon_signed_rank_z", wilcoxon["z"], "", "normal_approximation"),
        ("wilcoxon_signed_rank_two_sided_p", wilcoxon["p_value"], "", "normal_approximation"),
        ("wilcoxon_effect_r", wilcoxon["effect_r"], "", "z/sqrt(n_genera)"),
        (
            "stratified_permutation_unweighted_mean_delta",
            perm_unweighted,
            "",
            f"{args.n_permutations}_within_genus_permutations",
        ),
        (
            "stratified_permutation_unweighted_p",
            perm_p_unweighted,
            "",
            f"{args.n_permutations}_within_genus_permutations",
        ),
        (
            "stratified_permutation_weighted_mean_delta",
            perm_weighted,
            "",
            f"{args.n_permutations}_within_genus_permutations",
        ),
        (
            "stratified_permutation_weighted_p",
            perm_p_weighted,
            "",
            f"{args.n_permutations}_within_genus_permutations",
        ),
        (
            "genome_level_mean_growth_phaZ_positive",
            float(
                balanced.loc[
                    balanced["phaZ_status"] == "phaZ_positive", "growth_rate_per_h"
                ].mean()
            ),
            "",
            "descriptive_only_not_independent",
        ),
        (
            "genome_level_mean_growth_phaZ_negative",
            float(
                balanced.loc[
                    balanced["phaZ_status"] == "phaZ_negative", "growth_rate_per_h"
                ].mean()
            ),
            "",
            "descriptive_only_not_independent",
        ),
        (
            "genome_level_median_growth_phaZ_positive",
            float(
                balanced.loc[
                    balanced["phaZ_status"] == "phaZ_positive", "growth_rate_per_h"
                ].median()
            ),
            "",
            "descriptive_only_not_independent",
        ),
        (
            "genome_level_median_growth_phaZ_negative",
            float(
                balanced.loc[
                    balanced["phaZ_status"] == "phaZ_negative", "growth_rate_per_h"
                ].median()
            ),
            "",
            "descriptive_only_not_independent",
        ),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "interval", "note"])


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )


def draw_figure(
    args: argparse.Namespace,
    balanced: pd.DataFrame,
    summary: pd.DataFrame,
    effects: pd.DataFrame,
    tests: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(7.2, 6.6), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.55, 1.15, 1.05],
        width_ratios=[1.0, 1.05],
        hspace=0.55,
        wspace=0.38,
    )

    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])

    # Panel A: accounting.
    panel_label(ax0, "a")
    ax0.axis("off")
    n_genera = int(summary["genus"].nunique())
    n_pos = int((balanced["phaZ_status"] == "phaZ_positive").sum())
    n_neg = int((balanced["phaZ_status"] == "phaZ_negative").sum())
    total_run = int(tests.loc[tests["metric"] == "n_grodon2_input_genomes", "value"].iloc[0])
    failed = int(tests.loc[tests["metric"] == "n_failed_genomes", "value"].iloc[0])
    missing_growth = int(
        tests.loc[tests["metric"] == "n_status_ok_but_missing_growth_rate", "value"].iloc[0]
    )
    excluded = int(
        tests.loc[tests["metric"] == "n_successful_excluded_by_balance", "value"].iloc[0]
    )
    boxes = [
        ("gRodon2 run", f"{total_run:,} genomes"),
        ("Usable growth", f"{n_pos + n_neg + excluded:,} genomes"),
        ("Strict balance", f"{n_pos:,} phaZ+ / {n_neg:,} phaZ-"),
        ("Genus tests", f"{n_genera:,} genera"),
    ]
    x0 = 0.04
    for i, (title, value) in enumerate(boxes):
        x = x0 + i * 0.235
        ax0.add_patch(
            mpl.patches.FancyBboxPatch(
                (x, 0.28),
                0.19,
                0.48,
                boxstyle="round,pad=0.015,rounding_size=0.018",
                facecolor=PALETTE["very_light"],
                edgecolor="#BDBDBD",
                linewidth=0.8,
                transform=ax0.transAxes,
            )
        )
        ax0.text(x + 0.095, 0.60, title, ha="center", va="center", fontsize=7)
        ax0.text(
            x + 0.095,
            0.43,
            value,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
        if i < len(boxes) - 1:
            ax0.annotate(
                "",
                xy=(x + 0.228, 0.52),
                xytext=(x + 0.202, 0.52),
                xycoords=ax0.transAxes,
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=PALETTE["neutral"]),
            )
    ax0.text(
        0.04,
        0.12,
        f"Failed predictions: {failed}; missing growth rate: {missing_growth}; records removed only for balance: {excluded}.",
        transform=ax0.transAxes,
        ha="left",
        va="center",
        color=PALETTE["neutral"],
    )

    # Panel B: genus mean scatter.
    panel_label(ax1, "b")
    sizes = np.clip(summary["n_phaZ_positive"].to_numpy(float), 1, 80)
    sizes = 10 + 1.8 * np.sqrt(sizes)
    delta = summary["delta_positive_minus_negative"].to_numpy(float)
    sc = ax1.scatter(
        summary["mean_growth_phaZ_negative"],
        summary["mean_growth_phaZ_positive"],
        c=delta,
        cmap="RdBu_r",
        s=sizes,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.25,
    )
    upper = float(
        max(summary["mean_growth_phaZ_negative"].max(), summary["mean_growth_phaZ_positive"].max())
    )
    upper = math.ceil(upper * 10) / 10
    ax1.plot([0, upper], [0, upper], color="#8C8C8C", lw=0.8, ls="--")
    ax1.set_xlim(0, upper)
    ax1.set_ylim(0, upper)
    ax1.set_xlabel("Mean growth rate, phaZ- (h$^{-1}$)")
    ax1.set_ylabel("Mean growth rate, phaZ+ (h$^{-1}$)")
    cbar = fig.colorbar(sc, ax=ax1, fraction=0.046, pad=0.03)
    cbar.set_label("Delta (phaZ+ - phaZ-)")

    # Panel C: delta distribution.
    panel_label(ax2, "c")
    ax2.axvline(0, color="#8C8C8C", lw=0.8, ls="--")
    bins = np.linspace(
        np.percentile(delta, 1),
        np.percentile(delta, 99),
        34,
    )
    ax2.hist(delta, bins=bins, color="#9ECAE1", edgecolor="white", linewidth=0.3)
    median_delta = float(np.median(delta))
    ax2.axvline(median_delta, color=PALETTE["positive"], lw=1.4)
    ax2.set_xlabel("Genus-level delta growth rate (h$^{-1}$)")
    ax2.set_ylabel("Number of genera")
    p_wilcox = float(
        tests.loc[tests["metric"] == "wilcoxon_signed_rank_two_sided_p", "value"].iloc[0]
    )
    p_perm = float(
        tests.loc[tests["metric"] == "stratified_permutation_unweighted_p", "value"].iloc[0]
    )
    ax2.text(
        0.98,
        0.94,
        f"median = {median_delta:+.4f} h$^{{-1}}$\nWilcoxon P = {p_wilcox:.3g}\nPermutation P = {p_perm:.3g}",
        transform=ax2.transAxes,
        ha="right",
        va="top",
        fontsize=7,
    )

    # Panel D: sample-size distribution.
    panel_label(ax3, "d")
    n_each = summary["n_phaZ_positive"].astype(int)
    bins = np.arange(1, min(25, n_each.max()) + 2) - 0.5
    ax3.hist(n_each.clip(upper=25), bins=bins, color="#BDBDBD", edgecolor="white")
    ax3.set_xlabel("Selected genomes per status within genus")
    ax3.set_ylabel("Number of genera")
    n_one = int((n_each == 1).sum())
    n_two = int((n_each == 2).sum())
    ax3.text(
        0.97,
        0.94,
        f"1 vs 1: {n_one} genera\n2 vs 2: {n_two} genera",
        transform=ax3.transAxes,
        ha="right",
        va="top",
        fontsize=7,
    )

    # Panel E: top-sample genera paired means.
    panel_label(ax4, "e")
    top = effects.sort_values("n_phaZ_positive", ascending=False).head(16).iloc[::-1]
    y = np.arange(len(top))
    ax4.hlines(
        y,
        top["mean_growth_phaZ_negative"],
        top["mean_growth_phaZ_positive"],
        color="#BDBDBD",
        lw=0.9,
    )
    ax4.scatter(
        top["mean_growth_phaZ_negative"],
        y,
        color=PALETTE["negative"],
        s=14,
        label="phaZ-",
        zorder=3,
    )
    ax4.scatter(
        top["mean_growth_phaZ_positive"],
        y,
        color=PALETTE["positive"],
        s=14,
        label="phaZ+",
        zorder=3,
    )
    ax4.set_yticks(y)
    ax4.set_yticklabels(top["genus"])
    ax4.set_xlabel("Mean growth rate (h$^{-1}$)")
    ax4.legend(loc="lower right", handletextpad=0.3)

    fig.suptitle(
        "Genus-balanced gRodon2 comparison of phaZ-positive and phaZ-negative genomes",
        fontsize=9,
        y=0.985,
    )

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    stem = args.figure_dir / f"figure5_grodon_growth_comparison_{args.label}"
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_outputs(
    args: argparse.Namespace,
    balanced: pd.DataFrame,
    effects: pd.DataFrame,
    tests: pd.DataFrame,
) -> None:
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.source_dir.mkdir(parents=True, exist_ok=True)

    effects_path = args.tables_dir / f"grodon_growth_genus_effects_{args.label}.tsv"
    tests_path = args.tables_dir / f"grodon_growth_statistical_tests_{args.label}.tsv"
    effects.to_csv(effects_path, sep="\t", index=False)
    tests.to_csv(tests_path, sep="\t", index=False)

    # Lightweight source data for the publication figure.
    effects.to_csv(
        args.source_dir / f"figure5_grodon_genus_effects_{args.label}.tsv",
        sep="\t",
        index=False,
    )
    tests.to_csv(
        args.source_dir / f"figure5_grodon_statistical_tests_{args.label}.tsv",
        sep="\t",
        index=False,
    )
    balanced[
        [
            "genome_id",
            "phaZ_status",
            "major_subtype",
            "phylum",
            "genus",
            "species",
            "growth_rate_per_h",
            "doubling_time_h",
        ]
    ].to_csv(
        args.source_dir / f"figure5_grodon_balanced_genomes_{args.label}.tsv",
        sep="\t",
        index=False,
    )


def main() -> None:
    args = parse_args()
    balanced, summary, effects = prepare_tables(args)
    tests = run_statistics(args, balanced, summary, effects)
    write_outputs(args, balanced, effects, tests)
    draw_figure(args, balanced, summary, effects, tests)

    print(tests.to_string(index=False))
    print(f"Wrote Figure 5 to {args.figure_dir}")
    print(f"Wrote source data to {args.source_dir}")


if __name__ == "__main__":
    main()
