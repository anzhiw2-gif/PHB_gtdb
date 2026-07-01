#!/usr/bin/env python3
"""Live monitor for the PHB_gtdb gRodon2 growth prediction run.

The script reads the manifest, prediction table, and nohup log produced by
08_grodon_growth.py. It prints a live terminal dashboard and can optionally
write a small PNG progress plot at each refresh.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("PHB_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
LOGS_DIR = PROJECT_ROOT / "results" / "logs"


@dataclass
class Snapshot:
    timestamp: float
    total: int
    completed: int
    ok: int
    failed: int
    log_completed: int | None
    log_total: int | None
    last_log_lines: list[str]
    output_mtime: float | None


def count_manifest(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        n_lines = sum(1 for _ in handle)
    return max(0, n_lines - 1)


def read_predictions(path: Path) -> tuple[int, int, int]:
    if not path.exists():
        return 0, 0, 0
    completed = 0
    ok = 0
    failed = 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            completed += 1
            status = row.get("status", "")
            if status == "ok":
                ok += 1
            elif status:
                failed += 1
    return completed, ok, failed


def tail_lines(path: Path, n: int = 12) -> list[str]:
    if not path.exists():
        return []
    lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lines.append(line.rstrip("\n"))
            if len(lines) > n:
                lines.pop(0)
    return lines


def parse_log_progress(lines: list[str]) -> tuple[int | None, int | None]:
    pattern = re.compile(r"Processed\s+(\d+)/(\d+)\s+genomes")
    for line in reversed(lines):
        match = pattern.search(line)
        if match:
            return int(match.group(1)), int(match.group(2))
    return None, None


def snapshot(manifest: Path, predictions: Path, log: Path) -> Snapshot:
    total = count_manifest(manifest)
    completed, ok, failed = read_predictions(predictions)
    lines = tail_lines(log, n=14)
    log_completed, log_total = parse_log_progress(lines)
    output_mtime = predictions.stat().st_mtime if predictions.exists() else None
    if total == 0 and log_total:
        total = log_total
    return Snapshot(
        timestamp=time.time(),
        total=total,
        completed=completed,
        ok=ok,
        failed=failed,
        log_completed=log_completed,
        log_total=log_total,
        last_log_lines=lines,
        output_mtime=output_mtime,
    )


def progress_bar(done: int, total: int, width: int = 42) -> str:
    if total <= 0:
        return "[" + "." * width + "]"
    filled = min(width, max(0, round(width * done / total)))
    return "[" + "#" * filled + "." * (width - filled) + "]"


def fmt_seconds(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return "unknown"
    seconds = int(seconds)
    days, rem = divmod(seconds, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, sec = divmod(rem, 60)
    if days:
        return f"{days}d {hours:02d}h {minutes:02d}m"
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {sec:02d}s"
    return f"{sec}s"


def estimate_rate(history: list[Snapshot], window: int = 6) -> float | None:
    usable = [s for s in history if s.completed > 0]
    if len(usable) < 2:
        return None
    recent = usable[-window:]
    first = recent[0]
    last = recent[-1]
    dt = last.timestamp - first.timestamp
    dn = last.completed - first.completed
    if dt <= 0 or dn <= 0:
        return None
    return dn / dt


def append_history(path: Path, snap: Snapshot, rate_per_s: float | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if write_header:
            writer.writerow(
                [
                    "timestamp",
                    "datetime",
                    "total",
                    "completed",
                    "ok",
                    "failed",
                    "percent",
                    "rate_genomes_per_min",
                ]
            )
        percent = (100 * snap.completed / snap.total) if snap.total else 0
        rate_min = rate_per_s * 60 if rate_per_s else ""
        writer.writerow(
            [
                int(snap.timestamp),
                datetime.fromtimestamp(snap.timestamp).isoformat(timespec="seconds"),
                snap.total,
                snap.completed,
                snap.ok,
                snap.failed,
                f"{percent:.4f}",
                f"{rate_min:.4f}" if rate_min != "" else "",
            ]
        )


def write_plot(path: Path, history: list[Snapshot]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except Exception:
        return

    points = [s for s in history if s.total > 0]
    if not points:
        return

    times = [datetime.fromtimestamp(s.timestamp) for s in points]
    percents = [100 * s.completed / s.total for s in points]
    total = points[-1].total
    completed = points[-1].completed

    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=150)
    ax.plot(times, percents, color="#2f6f9f", lw=2.4)
    ax.scatter(times[-1], percents[-1], color="#d94841", s=32, zorder=3)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Completed (%)")
    ax.set_xlabel("Time")
    ax.set_title(f"gRodon2 progress: {completed}/{total} genomes")
    ax.grid(True, axis="y", color="#d9d9d9", lw=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fig.autofmt_xdate()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def clear_screen() -> None:
    sys.stdout.write("\033[2J\033[H")


def render(snap: Snapshot, history: list[Snapshot], args: argparse.Namespace) -> str:
    percent = (100 * snap.completed / snap.total) if snap.total else 0.0
    rate_per_s = estimate_rate(history)
    rate_per_min = rate_per_s * 60 if rate_per_s else None
    remaining = max(0, snap.total - snap.completed)
    eta_seconds = remaining / rate_per_s if rate_per_s else None
    eta_time = datetime.fromtimestamp(snap.timestamp) + timedelta(seconds=eta_seconds) if eta_seconds else None
    mtime = datetime.fromtimestamp(snap.output_mtime).strftime("%F %T") if snap.output_mtime else "not created"

    term_width = shutil.get_terminal_size((100, 24)).columns
    bar_width = max(24, min(56, term_width - 38))
    lines = [
        "PHB_gtdb gRodon2 live monitor",
        "=" * min(term_width, 88),
        f"Project       : {PROJECT_ROOT}",
        f"Manifest      : {args.manifest}",
        f"Predictions   : {args.predictions}",
        f"Log           : {args.log}",
        "",
        f"Progress      : {progress_bar(snap.completed, snap.total, bar_width)} {percent:6.2f}%",
        f"Completed     : {snap.completed}/{snap.total} genomes",
        f"Status        : ok={snap.ok}  failed={snap.failed}  pending={remaining}",
        f"Log checkpoint: {snap.log_completed if snap.log_completed is not None else '-'}"
        f"/{snap.log_total if snap.log_total is not None else '-'}",
        f"Rate          : {rate_per_min:.2f} genomes/min" if rate_per_min else "Rate          : collecting baseline...",
        f"ETA remaining : {fmt_seconds(eta_seconds)}",
        f"ETA finish    : {eta_time.strftime('%F %T') if eta_time else 'unknown'}",
        f"Output updated: {mtime}",
        f"Now           : {datetime.fromtimestamp(snap.timestamp).strftime('%F %T')}",
        "",
        "Recent log:",
    ]
    recent = snap.last_log_lines[-args.log_lines :]
    lines.extend(f"  {line}" for line in recent)
    lines.append("")
    lines.append("Press Ctrl+C to stop monitoring. The gRodon2 job itself will keep running.")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Live monitor for PHB_gtdb gRodon2 progress")
    parser.add_argument("--label", default="hmm_allmatched", help="Output label used by 08_grodon_growth.py")
    parser.add_argument("--interval", type=int, default=60, help="Refresh interval in seconds")
    parser.add_argument("--once", action="store_true", help="Print one snapshot and exit")
    parser.add_argument("--log-lines", type=int, default=10, help="Number of recent log lines to show")
    parser.add_argument("--history", default="", help="Optional TSV file to append monitor history")
    parser.add_argument("--plot", default="", help="Optional PNG path for live progress plot")
    parser.add_argument("--no-clear", action="store_true", help="Do not clear terminal between refreshes")
    parser.add_argument("--manifest", default="")
    parser.add_argument("--predictions", default="")
    parser.add_argument("--log", default="")
    return parser


def resolve_paths(args: argparse.Namespace) -> argparse.Namespace:
    label = args.label
    args.manifest = args.manifest or str(TABLES_DIR / f"grodon_growth_manifest_{label}.tsv")
    args.predictions = args.predictions or str(TABLES_DIR / f"grodon_growth_predictions_{label}.tsv")
    args.log = args.log or str(LOGS_DIR / f"grodon_growth_{label}.out")
    args.manifest = str(Path(args.manifest))
    args.predictions = str(Path(args.predictions))
    args.log = str(Path(args.log))
    return args


def main() -> int:
    args = resolve_paths(build_parser().parse_args())
    history: list[Snapshot] = []
    history_path = Path(args.history) if args.history else None
    plot_path = Path(args.plot) if args.plot else None

    try:
        while True:
            snap = snapshot(Path(args.manifest), Path(args.predictions), Path(args.log))
            history.append(snap)
            rate = estimate_rate(history)
            if history_path:
                append_history(history_path, snap, rate)
            if plot_path:
                write_plot(plot_path, history)
            if not args.no_clear:
                clear_screen()
            print(render(snap, history, args), flush=True)
            if args.once:
                break
            if snap.total > 0 and snap.completed >= snap.total:
                print("\nMonitoring complete: all manifest rows have prediction records.", flush=True)
                break
            time.sleep(max(1, args.interval))
    except KeyboardInterrupt:
        print("\nStopped monitor only. The gRodon2 job was not interrupted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
