#!/usr/bin/env python3
"""Check every PEGG input reference allele against the indexed local GRCh38 FASTA."""
import argparse
from pathlib import Path

import pandas as pd
import pysam


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=root / "GRCh38.fa")
    parser.add_argument("--input-dir", type=Path, default=root / "data/intermediate")
    parser.add_argument("--suffix", default="", help="Use .corrected for corrected PEGG inputs")
    parser.add_argument("--all-union-sites", action="store_true",
                        help="Also check every high-methylation union site is genomic A on RNA + or T on RNA -")
    args = parser.parse_args()
    fasta = pysam.FastaFile(str(args.fasta))
    contigs = {}
    for chrom in list(map(str, range(1, 23))) + ["X", "Y"]:
        number = {"X": 23, "Y": 24}.get(chrom, int(chrom) if chrom.isdigit() else 0)
        prefix = f"NC_{number:06d}."
        matches = [name for name in fasta.references if name.startswith(prefix)]
        if len(matches) != 1:
            raise SystemExit(f"Expected one FASTA contig for chr{chrom}; found {matches}")
        contigs[chrom] = matches[0]
    mito = [name for name in fasta.references if name.startswith("NC_012920.")]
    if len(mito) != 1:
        raise SystemExit(f"Expected one mitochondrial FASTA contig; found {mito}")
    contigs["M"] = mito[0]
    total_failures = 0
    for branch in ("intersection", "union"):
        input_path = args.input_dir / f"pegg_input_{branch}{args.suffix}.csv"
        rows = pd.read_csv(input_path, dtype=str, keep_default_na=False)
        failures = []
        for row in rows.itertuples(index=False):
            if row.Chromosome not in contigs:
                failures.append((row.Original_Pos_ID, "unmapped chromosome"))
                continue
            position = int(row.Start_Position)
            observed = fasta.fetch(contigs[row.Chromosome], position - 1, position).upper()
            if observed != row.Reference_Allele:
                failures.append((row.Original_Pos_ID, observed, row.Reference_Allele))
        print(f"{branch}: checked={len(rows)} reference mismatches={len(failures)}")
        if failures:
            print("examples:", failures[:5])
        total_failures += len(failures)
    if args.all_union_sites:
        union = pd.read_csv(root / "data/processed/HEK_HeLa_over80_union.csv",
                            usecols=["pos"], dtype=str, keep_default_na=False)
        if union.pos.eq("").any() or union.pos.duplicated().any():
            raise SystemExit("Union has missing/duplicate site IDs")
        failures = []
        for site in union.pos:
            try:
                chrom, position, strand = site.rsplit("_", 2)
                expected = {"+" : "A", "-": "T"}[strand]
                observed = fasta.fetch(contigs[chrom[3:]], int(position) - 1, int(position)).upper()
                if observed != expected:
                    failures.append((site, observed, expected))
            except (ValueError, KeyError):
                failures.append((site, "unparseable or unmapped"))
        print(f"all union sites: checked={len(union)} RNA-A reference mismatches={len(failures)}")
        if failures:
            print("examples:", failures[:5])
        total_failures += len(failures)
    if total_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
