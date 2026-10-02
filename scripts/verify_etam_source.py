#!/usr/bin/env python3
"""Verify the eTAM-seq TXT against the HeLa union table without editing either file."""
import argparse
from decimal import Decimal, InvalidOperation
from pathlib import Path

import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--raw",
        default=root / "eTAM-seq data" / "GSE211303_hela.polya.wt.ftom.ftop.rep1.deep.hits.txt",
        type=Path,
    )
    parser.add_argument(
        "--union", default=root / "data/processed/HeLa_union.csv", type=Path
    )
    args = parser.parse_args()
    raw = pd.read_csv(args.raw, sep="\t", dtype=str, keep_default_na=False)
    union = pd.read_csv(args.union, dtype=str, keep_default_na=False)
    sourced = union.loc[union.m6a_file3.ne(""), ["pos", "m6a_file3"]]
    if raw.pos.duplicated().any() or sourced.pos.duplicated().any():
        raise SystemExit("FAIL: duplicate pos in source or HeLa union")
    raw_map = dict(zip(raw.pos, raw.methylation))
    union_map = dict(zip(sourced.pos, sourced.m6a_file3))
    raw_sites, union_sites = set(raw_map), set(union_map)
    missing = sorted(raw_sites - union_sites)
    extra = sorted(union_sites - raw_sites)
    differing = []
    for pos in sorted(raw_sites & union_sites):
        try:
            same = Decimal(raw_map[pos]) == Decimal(union_map[pos])
        except InvalidOperation:
            same = False
        if not same:
            differing.append((pos, raw_map[pos], union_map[pos]))
    print(f"raw rows={len(raw)} unique={len(raw_sites)}")
    print(f"union m6a_file3 nonempty={len(sourced)} unique={len(union_sites)}")
    print(f"raw absent from union={len(missing)}; union absent from raw={len(extra)}")
    print(f"methylation mismatches={len(differing)}")
    if missing or extra or differing:
        print("examples:", missing[:3], extra[:3], differing[:3])
        raise SystemExit(1)
    print("PASS: all source positions and methylation values match one-to-one")


if __name__ == "__main__":
    main()
