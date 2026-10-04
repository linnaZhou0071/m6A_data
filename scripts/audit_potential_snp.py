#!/usr/bin/env python3
"""Audit the historical potential_SNP workbook against corrected CDS effects.

This is a read-only comparison. The local workbook and row-level effect CSVs
are excluded from the public repository; the script reports aggregate counts
without exporting their records.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "local/legacy/potential_SNP.xlsx"
EXAMPLE = "chr10_11462944_-"


def site_ids(frame: pd.DataFrame) -> set[str]:
    return set(frame["pos"].astype(str))


def audit(branch: str) -> None:
    workbook = pd.read_excel(WORKBOOK, sheet_name=branch, dtype={"pos": str})
    stem = f"HEK_HeLa_over80_{branch}"
    historical = pd.read_csv(
        ROOT / f"data/intermediate/{stem}_CDS.AtoG.effects.csv",
        usecols=["pos"],
        dtype=str,
    )
    corrected = pd.read_csv(
        ROOT / f"data/reanalysis/{stem}_CDS.AtoG.effects.corrected.csv",
        usecols=["pos", "analysis_transcript", "is_mRNA_A"],
        dtype=str,
    )

    positions = site_ids(workbook)
    site_strand = workbook.drop_duplicates("pos")["pos"].str.rsplit("_", n=1).str[-1]
    corrected_subset = corrected[corrected["pos"].isin(positions)]
    corrected_a = corrected_subset["is_mRNA_A"].str.lower().eq("true")
    old_false = workbook["is_mRNA_A"].astype(str).str.lower().eq("false")
    pair_columns = ["pos", "analysis_transcript"]
    matched_pairs = workbook[pair_columns].drop_duplicates().merge(
        corrected[pair_columns].drop_duplicates(), on=pair_columns,
        how="inner",
    )

    if not old_false.all():
        raise AssertionError(f"{branch}: workbook contains rows not labelled is_mRNA_A=False")
    if positions != site_ids(corrected_subset) or not corrected_a.all():
        raise AssertionError(f"{branch}: corrected effects do not cover all workbook sites as RNA A")

    print(
        f"{branch}: workbook rows={len(workbook):,}; unique sites={len(positions):,}; "
        f"negative-strand sites={int(site_strand.eq('-').sum()):,}; "
        f"old false rows={int(old_false.sum()):,}; "
        f"sites absent from old effects={len(positions - site_ids(historical)):,}; "
        f"sites present in corrected effects={len(site_ids(corrected_subset)):,}; "
        f"corrected true rows={int(corrected_a.sum()):,}; "
        f"same site-transcript pairs={len(matched_pairs):,}/{len(workbook):,}"
    )
    if branch == "intersection":
        example = workbook[workbook["pos"].eq(EXAMPLE)].iloc[0]
        print(
            f"Example {EXAMPLE}: transcript={example['analysis_transcript']}; "
            f"genomic bases in transcript order={example['genomic_codon']}; "
            f"old RNA codon={example['mrna_codon']}; old RNA base={example['mRNA_base']}"
        )


def main() -> None:
    if not WORKBOOK.is_file():
        raise FileNotFoundError(f"Local historical workbook not found: {WORKBOOK}")
    for branch in ("intersection", "union"):
        audit(branch)


if __name__ == "__main__":
    main()
