#!/usr/bin/env python3
"""Audit the local m⁶A tables and build small, reviewable release summaries.

Run from any directory: python scripts/build_release.py
This script never changes source/intermediate CSVs or full PEGG outputs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
INTERMEDIATE = ROOT / "data" / "intermediate"
SUMMARY = ROOT / "data" / "summary"
METADATA = ROOT / "metadata"


def table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def positions(df: pd.DataFrame, column: str = "pos") -> set[str]:
    values = df[column]
    if values.eq("").any() or values.duplicated().any():
        raise ValueError(f"{column} must be nonempty and unique")
    return set(values)


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    required_full = [
        ROOT / "pegg_results_intersection.csv",
        ROOT / "pegg_results_intersection_custom.csv",
    ]
    missing_full = [str(path) for path in required_full if not path.is_file()]
    if missing_full:
        raise SystemExit(
            "Cannot rebuild summaries without local full PEGG results: "
            + ", ".join(missing_full)
        )
    SUMMARY.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    hek = table(PROCESSED / "HEK_over_80%.csv")
    hela = table(PROCESSED / "HeLa_over_80%.csv")
    hek_sites, hela_sites = positions(hek), positions(hela)
    report: dict = {
        "threshold": "methylation >= 80%",
        "high_methylation_sites": {"HEK": len(hek_sites), "HeLa": len(hela_sites)},
        "branches": {},
        "release_status": "review_required",
    }
    for branch in ("intersection", "union"):
        stem = f"HEK_HeLa_over80_{branch}"
        original = table(PROCESSED / f"{stem}.csv")
        annotated = table(INTERMEDIATE / f"{stem}.annotated_v3.csv")
        cds = table(INTERMEDIATE / f"{stem}_CDS.v3fixed.csv")
        effects = table(INTERMEDIATE / f"{stem}_CDS.AtoG.effects.csv")
        mane = table(INTERMEDIATE / f"{stem}_CDS.AtoG.effects.MANE_nonsynonymous.csv")
        filtered = table(INTERMEDIATE / f"{stem}_filtered.csv")
        pegg = table(INTERMEDIATE / f"pegg_input_{branch}.csv")
        original_sites = positions(original)
        annotated_sites = positions(annotated)
        cds_sites = positions(cds)
        filtered_sites = positions(filtered)
        pegg_sites = positions(pegg, "Original_Pos_ID")
        effect_sites = set(effects.pos)
        mane_sites = set(mane.pos)
        any_nonsyn = set(effects.loc[effects.effect == "nonsynonymous", "pos"])
        expected = hek_sites & hela_sites if branch == "intersection" else hek_sites | hela_sites
        expected_cds = set(annotated.loc[annotated.genomic_feature == "CDS", "pos"])
        checks = {
            "set_definition": original_sites == expected,
            "annotation_covers_all_sites": annotated_sites == original_sites,
            "cds_matches_annotation": cds_sites == expected_cds,
            "effects_within_cds": effect_sites <= cds_sites,
            "filter_equals_original_minus_mane_nonsyn": filtered_sites == original_sites - mane_sites,
            "pegg_input_within_filtered": pegg_sites <= filtered_sites,
        }
        missing_effects = cds_sites - effect_sites
        counts = {
            "original_sites": len(original_sites),
            "cds_sites": len(cds_sites),
            "effect_rows": len(effects),
            "effect_sites": len(effect_sites),
            "cds_missing_effect_record": len(missing_effects),
            "cds_missing_effect_record_minus_strand": sum(x.endswith("_-") for x in missing_effects),
            "mane_nonsyn_sites_removed": len(mane_sites),
            "filtered_sites": len(filtered_sites),
            "filtered_sites_with_nonsyn_in_any_transcript": len(filtered_sites & any_nonsyn),
            "pegg_input_sites": len(pegg_sites),
            "filtered_not_in_pegg_input": sorted(filtered_sites - pegg_sites),
        }
        fate = annotated[["pos", "genomic_feature"]].copy()
        fate["has_effect_record"] = fate.pos.isin(effect_sites)
        fate["nonsynonymous_in_any_transcript"] = fate.pos.isin(any_nonsyn)
        fate["excluded_mane_nonsynonymous"] = fate.pos.isin(mane_sites)
        fate["retained_after_filter"] = fate.pos.isin(filtered_sites)
        fate["in_pegg_input"] = fate.pos.isin(pegg_sites)
        if branch == "intersection":
            for label, filename in (
                ("default", "pegg_results_intersection.csv"),
                ("custom", "pegg_results_intersection_custom.csv"),
            ):
                full_path = ROOT / filename
                full = table(full_path)
                full_sites = set(full.Original_Pos_ID)
                counts[f"{label}_designs"] = len(full)
                counts[f"{label}_sites_with_designs"] = len(full_sites)
                counts[f"{label}_input_sites_without_designs"] = len(pegg_sites - full_sites)
                checks[f"{label}_design_sites_within_input"] = full_sites <= pegg_sites
                checks[f"{label}_unique_pegRNA_id"] = not full.pegRNA_id.duplicated().any()
                per_site = full.groupby("Original_Pos_ID").size()
                fate[f"{label}_pegRNA_count"] = fate.pos.map(per_site).fillna(0).astype(int)
                rank_one = full.loc[pd.to_numeric(full.pegRNA_rank, errors="coerce").eq(1)]
                checks[f"{label}_one_rank1_per_designed_site"] = (
                    len(rank_one) == len(full_sites)
                    and not rank_one.Original_Pos_ID.duplicated().any()
                )
                selected = [
                    "Original_Pos_ID", "pegRNA_id", "Chromosome", "Start_Position",
                    "Reference_Allele", "Tumor_Seq_Allele2", "PAM_strand",
                    "Protospacer", "PAM", "RTT", "RTT_length", "PBS",
                    "PBS_length", "RTT_PBS", "RHA_size", "PEGG2_Score",
                    "RF_Score", "pegRNA_rank", "contains_polyT_terminator",
                ]
                rank_one[selected].to_csv(
                    SUMMARY / f"pegg_intersection_rank1_{label}.csv", index=False
                )
        fate.to_csv(SUMMARY / f"site_fate_{branch}.csv", index=False)
        report["branches"][branch] = {"counts": counts, "checks": checks}

    report["threshold_qc"] = {
        "HEK_exactly_80_percent": int(pd.to_numeric(hek.m6a_mean, errors="coerce").eq(80).sum())
        if "m6a_mean" in hek else "column not found",
        "HeLa_exactly_80_percent": int(pd.to_numeric(hela.m6a_mean, errors="coerce").eq(80).sum())
        if "m6a_mean" in hela else "column not found",
        "HEK_all_at_least_80_percent": bool(pd.to_numeric(hek.m6a_mean, errors="coerce").ge(80).all()),
        "HeLa_all_at_least_80_percent": bool(pd.to_numeric(hela.m6a_mean, errors="coerce").ge(80).all()),
    }
    # Hash only release files and the two locally retained, oversized PEGG outputs.
    files = [
        p for folder in ("data/processed", "data/intermediate", "data/summary", "figures")
        for p in (ROOT / folder).rglob("*") if p.is_file()
    ]
    files += [ROOT / name for name in (
        "pegg_results_intersection.csv", "pegg_results_intersection_custom.csv"
    ) if (ROOT / name).is_file()]
    with (METADATA / "checksums.sha256").open("w", encoding="utf-8") as handle:
        for path in sorted(files):
            handle.write(f"{file_hash(path)}  {path.relative_to(ROOT)}\n")
    (METADATA / "qc_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
