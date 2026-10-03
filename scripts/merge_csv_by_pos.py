import argparse
import csv
import os
from typing import Dict, List, Set, Tuple
import pandas as pd


def read_table_as_dict(path: str, key_col: str = "pos", exclude_cols=None) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    Read a CSV or spreadsheet as {pos: {column: value}} and return data columns.
    CSV files are read using utf-8-sig.
    Spreadsheets are read from the first sheet using pandas and openpyxl.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        df = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    elif ext in (".xlsx", ".xlsm", ".xls"):
        df = pd.read_excel(path, dtype=str, keep_default_na=False)  # first sheet
    else:
        raise ValueError(f"Unsupported file type: {ext} -> {path}")

    if key_col not in df.columns:
        raise ValueError(f"Input file lacks key column: {key_col} -> {path}")

    # Keep the first occurrence of each position; set operations use unique IDs.
    df = df.drop_duplicates(subset=[key_col])

    # Normalize values to strings so missing values do not become NaN keys.
    df = df.astype(str)
    df = df.replace({"nan": ""})

    cols = [c for c in df.columns if c != key_col and c not in (exclude_cols or set())]
    data: Dict[str, Dict[str, str]] = {}
    for _, row in df.iterrows():
        key = row[key_col]
        data[key] = {c: row[c] for c in cols}
    return data, cols


def write_union_intersection(
    d1: Dict[str, Dict[str, str]],
    d2: Dict[str, Dict[str, str]],
    cols1: List[str],
    cols2: List[str],
    label1: str,
    label2: str,
    out_prefix: str,
    key_col: str = "pos",
) -> Tuple[str, str]:
    """
    Write union and intersection CSVs and return their paths.
    Columns are the position key and source-suffixed columns from both inputs.
    """
    keys1: Set[str] = set(d1.keys())
    keys2: Set[str] = set(d2.keys())
    inter_keys = keys1 & keys2
    union_keys = keys1 | keys2

    # Build the output header.
    header_union = [key_col]
    header_union += [f"{c}_{label1}" for c in cols1]
    header_union += [f"{c}_{label2}" for c in cols2]

    # Write the union.
    union_path = f"{out_prefix}_union.csv"
    with open(union_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header_union)
        writer.writeheader()
        for k in sorted(union_keys):
            row: Dict[str, str] = {key_col: k}
            v1 = d1.get(k, {})
            v2 = d2.get(k, {})
            for c in cols1:
                row[f"{c}_{label1}"] = v1.get(c, "")
            for c in cols2:
                row[f"{c}_{label2}"] = v2.get(c, "")
            writer.writerow(row)

    # Write the intersection.
    inter_path = f"{out_prefix}_intersection.csv"
    with open(inter_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header_union)
        writer.writeheader()
        for k in sorted(inter_keys):
            row = {key_col: k}
            v1 = d1.get(k, {})
            v2 = d2.get(k, {})
            for c in cols1:
                row[f"{c}_{label1}"] = v1.get(c, "")
            for c in cols2:
                row[f"{c}_{label2}"] = v2.get(c, "")
            writer.writerow(row)

    # Report set cardinalities and the union identity check.
    print("====== Merge summary ======")
    print(f"Unique positions in input 1: {len(keys1)}")
    print(f"Unique positions in input 2: {len(keys2)}")
    print(f"Positions only in input 1: {len(keys1 - keys2)}")
    print(f"Positions only in input 2: {len(keys2 - keys1)}")
    print(f"Intersection positions: {len(inter_keys)}")
    print(f"Union positions: {len(union_keys)}")
    # Verify |A union B| = |A| + |B| - |A intersection B|.
    lhs = len(union_keys)
    rhs = len(keys1) + len(keys2) - len(inter_keys)
    print(f"Union identity check: {lhs} == {rhs} -> {'PASS' if lhs == rhs else 'FAIL'}")
    print(f"Union output: {union_path}")
    print(f"Intersection output: {inter_path}")

    return union_path, inter_path


def derive_label_from_path(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    # Use simple suffix characters to avoid spaces in output column names.
    safe = []
    for ch in stem:
        if ch.isalnum() or ch in ("_", "-"):
            safe.append(ch)
    label = "".join(safe) or "file"
    return label


def main():
    parser = argparse.ArgumentParser(
        description="Merge two tables by position into union and intersection CSVs with source-suffixed columns.",
    )
    parser.add_argument("--input1", "-i1", required=True, help="First input CSV or spreadsheet")
    parser.add_argument("--input2", "-i2", required=True, help="Second input CSV or spreadsheet")
    parser.add_argument("--label1", default=None, help="Suffix for input 1 columns (default: filename stem)")
    parser.add_argument("--label2", default=None, help="Suffix for input 2 columns (default: filename stem)")
    parser.add_argument("--out", "-o", default=None, help="Output prefix (default: label1_label2_pos)")
    parser.add_argument("--key", default="pos", help="Position key column (default: pos)")
    parser.add_argument("--exclude-columns", nargs="*", default=[], help="Auxiliary columns to omit from merged output")

    args = parser.parse_args()

    label1 = args.label1 or derive_label_from_path(args.input1)
    label2 = args.label2 or derive_label_from_path(args.input2)
    out_prefix = args.out or f"{label1}_{label2}_pos"

    d1, cols1 = read_table_as_dict(args.input1, key_col=args.key, exclude_cols=set(args.exclude_columns))
    d2, cols2 = read_table_as_dict(args.input2, key_col=args.key, exclude_cols=set(args.exclude_columns))

    write_union_intersection(d1, d2, cols1, cols2, label1, label2, out_prefix, key_col=args.key)


if __name__ == "__main__":
    main()
