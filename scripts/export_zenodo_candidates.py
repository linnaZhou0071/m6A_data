"""Export lean, unvalidated candidates; this does not rerun PEGG."""
from __future__ import annotations
import argparse
import csv
import hashlib
import os
import tempfile
from pathlib import Path
from typing import Optional

FIELDS = (
    "Original_Pos_ID", "pegRNA_id", "Chromosome", "Start_Position",
    "Reference_Allele", "Tumor_Seq_Allele2", "PAM_strand", "Protospacer",
    "PAM", "RTT", "RTT_length", "PBS", "PBS_length", "RTT_PBS",
    "RHA_size", "PEGG2_Score", "RF_Score", "pegRNA_rank",
    "contains_polyT_terminator",
)
EXPECTED = {"default": 129128, "custom": 248876}
EXPECTED_SITES = 3564


def export(source: Path, destination: Path, expected_rows: Optional[int] = None,
           expected_sites: Optional[int] = None) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    sites = set()
    guide_ids = set()
    temporary = None
    try:
        with source.open(newline="", encoding="utf-8-sig") as source_handle:
            reader = csv.DictReader(source_handle)
            missing = set(FIELDS) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{source}: missing columns {sorted(missing)}")
            with tempfile.NamedTemporaryFile(
                mode="w", newline="", encoding="utf-8", dir=destination.parent,
                prefix=".zenodo-export-", suffix=".csv", delete=False,
            ) as target_handle:
                temporary = Path(target_handle.name)
                writer = csv.DictWriter(target_handle, fieldnames=FIELDS,
                                        lineterminator="\n")
                writer.writeheader()
                for row in reader:
                    site = row["Original_Pos_ID"]
                    parts = site.rsplit("_", 2)
                    if len(parts) != 3 or parts[2] not in ("+", "-"):
                        raise ValueError(f"Invalid site ID: {site!r}")
                    chrom = parts[0].removeprefix("chr")
                    if chrom != row["Chromosome"] or parts[1] != row["Start_Position"]:
                        raise ValueError(f"Coordinate mismatch: {site!r}")
                    ref, alt = ("A", "G") if parts[2] == "+" else ("T", "C")
                    if (row["Reference_Allele"], row["Tumor_Seq_Allele2"]) != (ref, alt):
                        raise ValueError(f"Allele mismatch: {site!r}")
                    guide_id = row["pegRNA_id"]
                    if not guide_id or guide_id in guide_ids:
                        raise ValueError(f"Empty or duplicate pegRNA_id: {guide_id!r}")
                    if not row["Protospacer"] or not row["RTT"] or not row["PBS"]:
                        raise ValueError(f"Missing guide sequence: {guide_id!r}")
                    guide_ids.add(guide_id)
                    sites.add(site)
                    writer.writerow({field: row[field] for field in FIELDS})
                    count += 1
        if expected_rows is not None and count != expected_rows:
            raise ValueError(f"{source}: expected {expected_rows} rows, observed {count}")
        if expected_sites is not None and len(sites) != expected_sites:
            raise ValueError(f"{source}: expected {expected_sites} sites, observed {len(sites)}")
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    digest = hashlib.sha256()
    with destination.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return dict(file=str(destination), rows=count, sites=len(sites),
                bytes=destination.stat().st_size, sha256=digest.hexdigest())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path,
                        default=Path("data/reanalysis/full_guides"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("data/zenodo_preview"))
    args = parser.parse_args()
    exported = []
    for name, rows in EXPECTED.items():
        source = args.source_dir / f"pegg_intersection_{name}_corrected_filter.csv"
        destination = args.output_dir / f"pegg_intersection_{name}_candidate_v0.2.0.csv"
        result = export(source, destination, rows, EXPECTED_SITES)
        exported.append(result)
        print(result)
    manifest = args.output_dir / "SHA256SUMS.txt"
    manifest.write_text(
        "".join(
            "{}  {}\n".format(item["sha256"], Path(item["file"]).name)
            for item in exported
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
