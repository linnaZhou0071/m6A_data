#!/usr/bin/env python3
"""Hash and classify local candidate files without claiming publication clearance.

Run after the final edit. This inventory deliberately separates original code/
documentation from source-derived data awaiting rights review and large Zenodo
staging files. Source downloads and genome references are listed elsewhere.
"""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "metadata/local_file_inventory.tsv"
HASHES = ROOT / "metadata/local_candidate_checksums.sha256"


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def collect() -> list[tuple[Path, str]]:
    items = []
    for name in (".gitignore", "README.md", "LICENSE", "LICENSE-CC-BY-4.0.md", "CITATION.cff", "environment.yml"):
        items.append((ROOT / name, "original_code_or_documentation"))
    for folder in ("scripts", "tests"):
        items.extend((p, "original_code_or_documentation")
                     for p in (ROOT / folder).rglob("*.py"))
    items.extend((p, "original_code_or_documentation")
                 for p in (ROOT / "metadata").glob("*")
                 if p.is_file() and p.name not in {
                     OUT.name, HASHES.name, "checksums.sha256",
                 })
    items.extend((p, "source_derived_rights_hold")
                 for p in (ROOT / "data").rglob("*.csv")
                 if "full_guides" not in p.parts)
    items.extend((p, "original_code_or_documentation" if p.name == "workflow.svg"
                  else "source_derived_rights_hold")
                 for p in (ROOT / "figures").glob("*")
                 if p.suffix.lower() in {".svg", ".png", ".csv"})
    items.extend((p, "zenodo_candidate_rights_hold")
                 for p in (ROOT / "data/reanalysis/full_guides").glob("*.csv"))
    items.extend((p, "historical_local_not_current")
                 for p in ROOT.glob("pegg_results_*.csv"))
    items.extend((p, "historical_local_not_current")
                 for p in (ROOT / "figures").glob("*.jpg"))
    # The user's renamed local Prism project; it is not published.
    prism = ROOT / "m6A_plotting.pzfx"
    if prism.exists():
        items.append((prism, "historical_local_not_current"))
    if len(items) != len({p for p, _ in items}):
        raise AssertionError("Duplicate inventory path")
    return sorted(items, key=lambda item: str(item[0].relative_to(ROOT)))


def main() -> None:
    items = collect()
    rows = []
    for path, status in items:
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({
            "relative_path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "publication_status": status,
        })
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with HASHES.open("w", encoding="utf-8") as handle:
        for row in rows:
            if row["publication_status"] not in {
                "historical_local_not_current", "zenodo_candidate_rights_hold"
            }:
                handle.write(f'{row["sha256"]}  {row["relative_path"]}\n')
    from collections import Counter
    print(f"Inventory: {len(rows)} files; statuses={dict(Counter(r['publication_status'] for r in rows))}")
    print(f"Written {OUT.relative_to(ROOT)} and {HASHES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
