import argparse
import csv
import os
from typing import Dict, List, Set, Tuple
import pandas as pd


def read_table_as_dict(path: str, key_col: str = "pos", exclude_cols=None) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    读取CSV或XLSX为 {pos: {col: value}} 的字典，并返回列名列表（不含 key_col）。
    - .csv：按 utf-8-sig 读取
    - .xlsx：使用 pandas + openpyxl 读取第一张表
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        df = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    elif ext in (".xlsx", ".xlsm", ".xls"):
        df = pd.read_excel(path, dtype=str, keep_default_na=False)  # 第一张表
    else:
        raise ValueError(f"不支持的文件类型: {ext} -> {path}")

    if key_col not in df.columns:
        raise ValueError(f"输入文件缺少关键列: {key_col} -> {path}")

    # 去除可能的重复pos，保留首条（交/并集按唯一pos进行）
    df = df.drop_duplicates(subset=[key_col])

    # 统一转换为字符串，避免 NaN 干扰（更兼容的写法）
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
    生成并写出并集与交集CSV，返回 (union_path, inter_path)。
    输出列为：pos + cols1加后缀 + cols2加后缀。
    """
    keys1: Set[str] = set(d1.keys())
    keys2: Set[str] = set(d2.keys())
    inter_keys = keys1 & keys2
    union_keys = keys1 | keys2

    # 构建表头
    header_union = [key_col]
    header_union += [f"{c}_{label1}" for c in cols1]
    header_union += [f"{c}_{label2}" for c in cols2]

    # 写并集
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

    # 写交集
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

    # 打印统计信息
    print("====== 合并结果统计 ======")
    print(f"文件1去重pos总数: {len(keys1)}")
    print(f"文件2去重pos总数: {len(keys2)}")
    print(f"文件1唯一pos数(仅在1中): {len(keys1 - keys2)}")
    print(f"文件2唯一pos数(仅在2中): {len(keys2 - keys1)}")
    print(f"交集pos数(两者共有): {len(inter_keys)}")
    print(f"并集pos数: {len(union_keys)}")
    # 验证恒等式：|A∪B| = |A| + |B| − |A∩B|
    lhs = len(union_keys)
    rhs = len(keys1) + len(keys2) - len(inter_keys)
    print(f"并集恒等式校验: {lhs} == {rhs} -> {'PASS' if lhs == rhs else 'FAIL'}")
    print(f"并集输出: {union_path}")
    print(f"交集输出: {inter_path}")

    return union_path, inter_path


def derive_label_from_path(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    # 仅保留简单字符作为后缀，避免列名过长或包含空格
    safe = []
    for ch in stem:
        if ch.isalnum() or ch in ("_", "-"):
            safe.append(ch)
    label = "".join(safe) or "file"
    return label


def main():
    parser = argparse.ArgumentParser(
        description="按pos对两个CSV做并集与交集合并，保留两侧全部列，并加来源后缀。",
    )
    parser.add_argument("--input1", "-i1", required=True, help="CSV 文件1路径")
    parser.add_argument("--input2", "-i2", required=True, help="CSV 文件2路径")
    parser.add_argument("--label1", default=None, help="文件1来源后缀(默认取文件名stem)")
    parser.add_argument("--label2", default=None, help="文件2来源后缀(默认取文件名stem)")
    parser.add_argument("--out", "-o", default=None, help="输出前缀(默认 label1_label2_pos)")
    parser.add_argument("--key", default="pos", help="主键列名(默认: pos)")
    parser.add_argument("--exclude-columns", nargs="*", default=[], help="不写入合并结果的辅助列")

    args = parser.parse_args()

    label1 = args.label1 or derive_label_from_path(args.input1)
    label2 = args.label2 or derive_label_from_path(args.input2)
    out_prefix = args.out or f"{label1}_{label2}_pos"

    d1, cols1 = read_table_as_dict(args.input1, key_col=args.key, exclude_cols=set(args.exclude_columns))
    d2, cols2 = read_table_as_dict(args.input2, key_col=args.key, exclude_cols=set(args.exclude_columns))

    write_union_intersection(d1, d2, cols1, cols2, label1, label2, out_prefix, key_col=args.key)


if __name__ == "__main__":
    main()
