#!/usr/bin/env python3
"""Audit corrected MANE-filtered sites and subset historical pegRNAs.

The coding-effect correction changes *site eligibility*, not the PEGG sequence
for a retained site. This script never reruns PEGG or overwrites historical data.
Use --full-guides to write large local CSVs for a possible Zenodo deposit.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from filter_non_mane_nonsynonymous import clean_transcript_id, parse_mane_transcripts

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data/processed"
INTERMEDIATE = ROOT / "data/intermediate"
REANALYSIS = ROOT / "data/reanalysis"
SUMMARY = ROOT / "data/summary"
METADATA = ROOT / "metadata"
GUIDES = (
    ("default", "pegg_results_intersection.csv"),
    ("custom", "pegg_results_intersection_custom.csv"),
)
RANK_COLUMNS = [
    "Original_Pos_ID", "pegRNA_id", "Chromosome", "Start_Position",
    "Reference_Allele", "Tumor_Seq_Allele2", "PAM_strand",
    "Protospacer", "PAM", "RTT", "RTT_length", "PBS",
    "PBS_length", "RTT_PBS", "RHA_size", "PEGG2_Score",
    "RF_Score", "pegRNA_rank", "contains_polyT_terminator",
]


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def site_set(frame: pd.DataFrame, column: str = "pos") -> set[str]:
    if frame[column].eq("").any() or frame[column].duplicated().any():
        raise AssertionError(f"Missing/duplicate {column}")
    return set(frame[column])


def subset_full(source: Path, target: Path, keep: set[str]) -> tuple[int, set[str]]:
    target.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    seen = set()
    with source.open(newline="") as src, target.open("w", newline="") as dst:
        reader = csv.DictReader(src)
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            site = row["Original_Pos_ID"]
            if site in keep:
                writer.writerow(row)
                n += 1
                seen.add(site)
    return n, seen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-guides", action="store_true",
                        help="Write large corrected-filter subsets for later Zenodo review")
    args = parser.parse_args()
    mane = parse_mane_transcripts(ROOT / "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz")
    coverage = read(REANALYSIS / "historical_pegg_coverage_audit.csv")
    coverage_by_run = {
        name: frame.set_index("Original_Pos_ID") for name, frame in coverage.groupby("historical_run")
    }
    report = {"scope": "corrected_MANE_only_candidate_snapshot", "branches": {}}
    for branch in ("intersection", "union"):
        stem = f"HEK_HeLa_over80_{branch}"
        original = read(PROCESSED / f"{stem}.csv")
        annotated = read(INTERMEDIATE / f"{stem}.annotated_v3.csv")
        cds = read(INTERMEDIATE / f"{stem}_CDS.v3fixed.csv")
        effects = read(REANALYSIS / f"{stem}_CDS.AtoG.effects.corrected.csv")
        filtered = read(REANALYSIS / f"{stem}_filtered.corrected.csv")
        historical_filtered = read(INTERMEDIATE / f"{stem}_filtered.csv")
        pegg = read(REANALYSIS / f"pegg_input_{branch}.corrected.csv")
        original_sites = site_set(original)
        cds_sites = site_set(cds)
        filtered_sites = site_set(filtered)
        historical_filtered_sites = site_set(historical_filtered)
        if not filtered_sites <= historical_filtered_sites:
            raise AssertionError(f"{branch}: corrected filter unexpectedly adds historical exclusions")
        pegg_sites = site_set(pegg, "Original_Pos_ID")
        if site_set(annotated) != original_sites:
            raise AssertionError(f"{branch}: annotation does not cover source sites")
        if set(effects.pos) != cds_sites:
            raise AssertionError(f"{branch}: coding effects do not cover every CDS site")
        if not effects.analysis_result.eq("ok").all() or not effects.is_mRNA_A.eq("True").all():
            raise AssertionError(f"{branch}: failed/invalid coding effect")
        mane_mask = effects.analysis_transcript.map(clean_transcript_id).isin(mane)
        mane_sites = set(effects.loc[mane_mask, "pos"])
        excluded = set(effects.loc[mane_mask & effects.effect.eq("nonsynonymous"), "pos"])
        any_nonsyn = set(effects.loc[effects.effect.eq("nonsynonymous"), "pos"])
        if filtered_sites != original_sites - excluded:
            raise AssertionError(f"{branch}: corrected MANE subtraction differs")
        skipped = filtered_sites - pegg_sites
        if skipped != {p for p in filtered_sites if p.startswith("chrM_")}:
            raise AssertionError(f"{branch}: unexplained PEGG input attrition")
        fate = annotated[["pos", "genomic_feature"]].copy()
        fate["in_cds"] = fate.pos.isin(cds_sites)
        fate["cds_effect_rows"] = fate.pos.map(effects.groupby("pos").size()).fillna(0).astype(int)
        fate["has_mane_effect"] = fate.pos.isin(mane_sites)
        fate["nonsynonymous_in_any_transcript"] = fate.pos.isin(any_nonsyn)
        fate["excluded_mane_nonsynonymous"] = fate.pos.isin(excluded)
        fate["retained_after_historical_filter"] = fate.pos.isin(historical_filtered_sites)
        fate["retained_after_corrected_filter"] = fate.pos.isin(filtered_sites)
        fate["newly_excluded_after_strand_correction"] = fate.pos.isin(
            historical_filtered_sites - filtered_sites
        )
        fate["in_corrected_pegg_input"] = fate.pos.isin(pegg_sites)
        fate["corrected_design_status"] = "excluded_before_design"
        fate.loc[fate.pos.isin(skipped), "corrected_design_status"] = "unsupported_chrM"
        if branch == "union":
            fate.loc[fate.in_corrected_pegg_input, "corrected_design_status"] = "not_attempted"
        fate["coding_review_status"] = "non_CDS_by_annotation"
        fate.loc[fate.in_cds & fate.has_mane_effect, "coding_review_status"] = "MANE_synonymous"
        fate.loc[fate.in_cds & ~fate.has_mane_effect, "coding_review_status"] = "CDS_no_MANE_effect"
        fate.loc[fate.excluded_mane_nonsynonymous, "coding_review_status"] = "MANE_nonsynonymous_excluded"
        counts = {
            "input_sites": len(original_sites),
            "cds_sites": len(cds_sites),
            "effect_rows": len(effects),
            "excluded_mane_nonsynonymous_sites": len(excluded),
            "retained_sites": len(filtered_sites),
            "newly_excluded_after_strand_correction": len(historical_filtered_sites - filtered_sites),
            "retained_cds_without_mane_effect": len((cds_sites & filtered_sites) - mane_sites),
            "retained_with_nonmane_nonsynonymous_effect": len(filtered_sites & (any_nonsyn - excluded)),
            "pegg_input_sites": len(pegg_sites),
            "unsupported_chrM_sites": sorted(skipped),
        }
        if branch == "intersection":
            old_input = site_set(read(INTERMEDIATE / "pegg_input_intersection.csv"), "Original_Pos_ID")
            if not pegg_sites <= old_input:
                raise AssertionError("Corrected PEGG input includes a site never designed historically")
            for label, _ in GUIDES:
                run = "default_like" if label == "default" else label
                lookup = coverage_by_run[run]
                if set(lookup.index) != old_input:
                    raise AssertionError(f"{label}: incomplete coverage ledger")
                fate[f"{label}_historical_pegRNA_count"] = (
                    fate.pos.map(lookup.historical_pegRNA_count).fillna("0").astype(int)
                )
                fate[f"{label}_historical_coverage_reason"] = (
                    fate.pos.map(lookup.coverage_reason).fillna("not_in_historical_PEGG_input")
                )
                rank = read(SUMMARY / f"pegg_intersection_rank1_{label}.csv")
                selected = rank[rank.Original_Pos_ID.isin(pegg_sites)].copy()
                expected = {s for s in pegg_sites if int(lookup.loc[s, "historical_pegRNA_count"]) > 0}
                if label == "custom":
                    fate.loc[fate.in_corrected_pegg_input, "corrected_design_status"] = (
                        "no_guide_under_historical_custom_parameters"
                    )
                    fate.loc[fate.pos.isin(expected), "corrected_design_status"] = "historical_guide_reused"
                if site_set(selected, "Original_Pos_ID") != expected:
                    raise AssertionError(f"{label}: rank-1 sites do not match historical guide coverage")
                if not selected.pegRNA_rank.astype(float).eq(1).all():
                    raise AssertionError(f"{label}: non-rank-1 guide in summary")
                selected[RANK_COLUMNS].to_csv(
                    SUMMARY / f"pegg_intersection_rank1_{label}_corrected_filter.csv", index=False
                )
                counts[f"{label}_retained_sites_with_guides"] = len(expected)
                counts[f"{label}_retained_sites_without_guides"] = len(pegg_sites - expected)
                expected_rows = sum(int(lookup.loc[s, "historical_pegRNA_count"]) for s in pegg_sites)
                counts[f"{label}_guide_rows_in_corrected_subset"] = expected_rows
                if args.full_guides:
                    n, sites = subset_full(
                        ROOT / dict(GUIDES)[label],
                        REANALYSIS / "full_guides" / f"pegg_intersection_{label}_corrected_filter.csv",
                        pegg_sites,
                    )
                    if sites != expected or n != expected_rows:
                        raise AssertionError(f"{label}: full-guide subset differs from coverage audit")
                    counts[f"{label}_full_guide_rows"] = n
        fate.to_csv(SUMMARY / f"site_fate_{branch}_corrected.csv", index=False)
        report["branches"][branch] = counts
    METADATA.mkdir(exist_ok=True)
    (METADATA / "corrected_qc_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
