#!/usr/bin/env python3
"""
Create genus-balanced gRodon2 analysis tables.

For each genus, keep the same number of phaZ-positive and phaZ-negative
genomes among successful gRodon2 predictions. This avoids within-genus
comparisons being driven by unequal sample counts after failed predictions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize gRodon2 failures and create genus-balanced tables."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/tables/grodon_growth_predictions_hmm_allmatched.tsv"),
        help="gRodon2 prediction table.",
    )
    parser.add_argument(
        "--label",
        default="hmm_allmatched",
        help="Output label used in generated table names.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/tables"),
        help="Directory for output tables.",
    )
    parser.add_argument(
        "--sort-by",
        choices=["genome_id", "growth_rate_per_h"],
        default="genome_id",
        help="Deterministic within-genus selection order.",
    )
    return parser.parse_args()


def write_failures(df: pd.DataFrame, outdir: Path, label: str) -> None:
    failed = df[df["status"] != "ok"].copy()

    failure_cols = [
        "genome_id",
        "phaZ_status",
        "phaZ_validated_count",
        "major_subtype",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "status",
        "error",
        "genome_path",
    ]
    failed[failure_cols].to_csv(
        outdir / f"grodon_failed_genomes_{label}.tsv",
        sep="\t",
        index=False,
    )

    failure_summary = (
        failed.groupby(["genus", "phaZ_status", "error"], dropna=False)
        .size()
        .reset_index(name="n_failed")
        .sort_values(["n_failed", "genus", "phaZ_status"], ascending=[False, True, True])
    )
    failure_summary.to_csv(
        outdir / f"grodon_failed_genera_summary_{label}.tsv",
        sep="\t",
        index=False,
    )


def select_balanced(df: pd.DataFrame, sort_by: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    ok = df[df["status"] == "ok"].copy()
    ok["growth_rate_per_h"] = pd.to_numeric(ok["growth_rate_per_h"], errors="coerce")
    ok = ok[ok["growth_rate_per_h"].notna()].copy()

    selected_parts = []
    genus_rows = []

    for genus, genus_df in ok.groupby("genus", sort=True):
        pos = genus_df[genus_df["phaZ_status"] == "phaZ_positive"].copy()
        neg = genus_df[genus_df["phaZ_status"] == "phaZ_negative"].copy()
        n_pos = len(pos)
        n_neg = len(neg)
        n_each = min(n_pos, n_neg)

        genus_rows.append(
            {
                "genus": genus,
                "n_phaZ_positive_ok": n_pos,
                "n_phaZ_negative_ok": n_neg,
                "n_each_selected": n_each,
                "n_excluded_positive": n_pos - n_each,
                "n_excluded_negative": n_neg - n_each,
                "included_in_balanced": n_each > 0,
            }
        )

        if n_each == 0:
            continue

        if sort_by == "growth_rate_per_h":
            pos = pos.sort_values(["growth_rate_per_h", "genome_id"], ascending=[False, True])
            neg = neg.sort_values(["growth_rate_per_h", "genome_id"], ascending=[False, True])
        else:
            pos = pos.sort_values("genome_id")
            neg = neg.sort_values("genome_id")

        selected_parts.append(pos.head(n_each))
        selected_parts.append(neg.head(n_each))

    selected = pd.concat(selected_parts, ignore_index=True) if selected_parts else ok.iloc[0:0]
    balance_audit = pd.DataFrame(genus_rows).sort_values(
        ["included_in_balanced", "n_each_selected", "genus"],
        ascending=[False, False, True],
    )
    return selected, balance_audit


def summarize_balanced(selected: pd.DataFrame) -> pd.DataFrame:
    summary_rows = []
    for genus, genus_df in selected.groupby("genus", sort=True):
        pos = genus_df[genus_df["phaZ_status"] == "phaZ_positive"]
        neg = genus_df[genus_df["phaZ_status"] == "phaZ_negative"]
        summary_rows.append(
            {
                "genus": genus,
                "n_phaZ_positive": len(pos),
                "n_phaZ_negative": len(neg),
                "mean_growth_phaZ_positive": pos["growth_rate_per_h"].mean(),
                "mean_growth_phaZ_negative": neg["growth_rate_per_h"].mean(),
                "delta_positive_minus_negative": (
                    pos["growth_rate_per_h"].mean() - neg["growth_rate_per_h"].mean()
                ),
                "median_growth_phaZ_positive": pos["growth_rate_per_h"].median(),
                "median_growth_phaZ_negative": neg["growth_rate_per_h"].median(),
            }
        )
    return pd.DataFrame(summary_rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    write_failures(df, args.outdir, args.label)

    selected, balance_audit = select_balanced(df, args.sort_by)
    summary = summarize_balanced(selected)

    selected.to_csv(
        args.outdir / f"grodon_growth_balanced_by_genus_{args.label}.tsv",
        sep="\t",
        index=False,
    )
    summary.to_csv(
        args.outdir / f"grodon_growth_balanced_by_genus_summary_{args.label}.tsv",
        sep="\t",
        index=False,
    )
    balance_audit.to_csv(
        args.outdir / f"grodon_growth_balance_audit_{args.label}.tsv",
        sep="\t",
        index=False,
    )

    print(f"input_rows={len(df)}")
    print(f"ok_rows={(df['status'] == 'ok').sum()}")
    print(f"failed_rows={(df['status'] != 'ok').sum()}")
    print(f"balanced_rows={len(selected)}")
    print(f"balanced_genera={summary.shape[0]}")
    print(f"selected_each_total={len(selected) // 2}")


if __name__ == "__main__":
    main()
