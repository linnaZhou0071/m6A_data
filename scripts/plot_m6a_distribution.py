#!/usr/bin/env python3
"""Rebuild 10-percent m6A histograms from the versioned HEK/HeLa union CSVs.

Bins are [0,10), [10,20), ..., [80,90), [90,100]; 100 is included
in the final bin. Run: python scripts/plot_m6a_distribution.py
The original Prism project and JPEGs are historical and are not overwritten.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
PROCESSED = ROOT / "data" / "processed"
EDGES = np.arange(0, 101, 10)
LABELS = [f"{i}–{i + 10}%" for i in range(0, 100, 10)]
SOURCES = (
    ("HEK group (SAC HEK293 + GLORI HEK293T)", "HEK_union.csv", "#39739d"),
    ("HeLa group", "HeLa_union.csv", "#b45e6e"),
)


def histogram(path: Path) -> tuple[np.ndarray, int]:
    frame = pd.read_csv(path, usecols=["pos", "m6a_mean"])
    if frame.pos.isna().any() or frame.pos.duplicated().any():
        raise ValueError(f"Missing or duplicate site ID in {path}")
    values = pd.to_numeric(frame.m6a_mean, errors="raise").to_numpy()
    if not np.isfinite(values).all() or ((values < 0) | (values > 100)).any():
        raise ValueError(f"Invalid percentage in {path}")
    counts, _ = np.histogram(values, bins=EDGES)
    if int(counts.sum()) != len(frame):
        raise AssertionError(f"Histogram lost rows from {path}")
    return counts, len(frame)


def main() -> None:
    rows = []
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout="constrained")
    for ax, (name, filename, color) in zip(axes, SOURCES):
        counts, n_sites = histogram(PROCESSED / filename)
        ax.bar(range(10), counts, width=0.82, color=color, edgecolor="white")
        ax.set_ylabel("Unique sites")
        ax.set_title(f"{name} · n={n_sites:,}", loc="left", fontsize=11)
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)
        for index, count in enumerate(counts):
            rows.append({
                "group": "HEK" if filename.startswith("HEK") else "HeLa",
                "source_table": filename,
                "bin_lower_inclusive_percent": int(EDGES[index]),
                "bin_upper_percent": int(EDGES[index + 1]),
                "count": int(count),
                "total_unique_sites": n_sites,
            })
    axes[-1].set_xticks(range(10), LABELS)
    axes[-1].set_xlabel("Mean m6A level over available source measurements")
    fig.suptitle("Distribution of curated single-base m6A measurements", fontsize=13)
    fig.text(0.99, 0.005, "10% bins; final bin includes 100%. Groups are not identical cell-line labels.",
             ha="right", fontsize=8)
    FIGURES.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(FIGURES / "m6a_distribution_counts.csv", index=False)
    for suffix in ("svg", "png"):
        fig.savefig(FIGURES / f"m6a_distribution_corrected.{suffix}", dpi=220)
    plt.close(fig)
    for name, _, _ in SOURCES:
        subset = [r["count"] for r in rows if r["group"] == name.split()[0]]
        print(name, subset, sum(subset))


if __name__ == "__main__":
    main()
