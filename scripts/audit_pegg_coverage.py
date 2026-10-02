#!/usr/bin/env python3
"""Explain historical no-guide sites using PEGG 2.1.0 PAM/RTT/RHA rules.

This reproduces candidate eligibility, not PEGG scoring or a new design run.
It checks predicted presence/absence against both historical full output CSVs.
"""
import csv
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

from predict_AtoG_coding_effects_v2 import (
    build_refseq_base_index, load_fasta_contigs, resolve_fasta_contig, revcomp,
)

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "GRCh38.fa"
INPUT = ROOT / "data/intermediate/pegg_input_intersection.csv"
OUT = ROOT / "data/reanalysis/historical_pegg_coverage_audit.csv"
RUNS = (
    ("default_like", ROOT / "pegg_results_intersection.csv", (5, 10, 15, 25, 30), 1),
    ("custom", ROOT / "pegg_results_intersection_custom.csv", (20, 22, 25, 27, 30, 32, 35), 6),
)
PAM_RE = re.compile(r"(?=([ACGT]GG))")
CONTEXT = 120
PROTO_SIZE = 19


def guide_counts(path: Path) -> Counter:
    counts = Counter()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            counts[row["Original_Pos_ID"]] += 1
    return counts


def classify(wt: str, rtt_lengths: tuple[int, ...], min_rha: int) -> tuple[int, int, str]:
    # SNP at index 120 in the 241-nt reference context.
    if len(wt) != 2 * CONTEXT + 1:
        raise ValueError("Incomplete sequence context")
    eligible = 0
    valid = 0
    for sequence in (wt, revcomp(wt)):
        left_len = CONTEXT
        low = max(left_len - (max(rtt_lengths) - 1 - 3), PROTO_SIZE + 3)
        high = left_len + 3
        for match in PAM_RE.finditer(sequence):
            pam_start = match.start()
            if not low <= pam_start <= high:
                continue
            eligible += 1
            left_rtt_length = left_len - (pam_start - 3)
            for rtt_length in rtt_lengths:
                rha_length = rtt_length - left_rtt_length - 1
                if rha_length >= min_rha and left_len + 1 + rha_length <= len(sequence):
                    valid += 1
    if eligible == 0:
        return eligible, valid, "no_eligible_NGG_PAM"
    if valid == 0:
        return eligible, valid, "no_RTT_meets_min_RHA"
    return eligible, valid, "candidate_eligible"


def main() -> None:
    inputs = pd.read_csv(INPUT, dtype=str, keep_default_na=False)
    if inputs.Original_Pos_ID.duplicated().any():
        raise ValueError("Duplicate input site")
    corrected = set(pd.read_csv(
        ROOT / "data/reanalysis/pegg_input_intersection.corrected.csv", usecols=["Original_Pos_ID"]
    ).Original_Pos_ID)
    input_sites = set(inputs.Original_Pos_ID)
    if not corrected <= input_sites:
        raise AssertionError("Corrected sites must be a subset of the historical design input")
    contigs = set(load_fasta_contigs(str(FASTA)))
    aliases = build_refseq_base_index(list(contigs))
    observed = {name: guide_counts(path) for name, path, _, _ in RUNS}
    if any(not set(counts) <= input_sites for counts in observed.values()):
        raise AssertionError("Historical guide output contains unexpected site IDs")
    rows = []
    with pysam.FastaFile(str(FASTA)) as fasta:
        for source in inputs.itertuples(index=False):
            site = source.Original_Pos_ID
            chrom = resolve_fasta_contig("chr" + source.Chromosome, contigs, aliases)
            if chrom is None:
                raise ValueError(f"Unmapped chromosome at {site}")
            pos = int(source.Start_Position)
            wt = fasta.fetch(chrom, pos - CONTEXT - 1, pos + CONTEXT).upper()
            if wt[CONTEXT] != source.Reference_Allele:
                raise AssertionError(f"Reference mismatch at {site}")
            for name, _, rtt_lengths, min_rha in RUNS:
                pam_count, valid, reason = classify(wt, rtt_lengths, min_rha)
                count = observed[name][site]
                if (reason == "candidate_eligible") != (count > 0):
                    raise AssertionError(f"Eligibility does not match observed output: {site}, {name}")
                rows.append({
                    "Original_Pos_ID": site,
                    "historical_run": name,
                    "retained_in_corrected_filter": site in corrected,
                    "eligible_NGG_PAMs": pam_count,
                    "valid_PAM_RTT_combinations": valid,
                    "historical_pegRNA_count": count,
                    "coverage_reason": reason,
                })
    OUT.parent.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(pd.DataFrame(rows).groupby(["historical_run", "coverage_reason"]).size().to_string())
    print(f"Wrote {OUT.relative_to(ROOT)}; all eligibility predictions match observed site coverage")


if __name__ == "__main__":
    main()
