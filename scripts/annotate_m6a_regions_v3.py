#!/usr/bin/env python3
"""
annotate_m6a_regions_v3.py

Usage:
  python annotate_m6a_regions_v3.py --gtf hg38.ncbiRefSeq.gtf.gz FILE1.csv [FILE2.csv ...]
  python annotate_m6a_regions_v3.py --gtf hg38.ncbiRefSeq.gtf.gz --verify FILE.csv

Outputs per CSV:
  <orig>.annotated_v3.csv  (input columns plus genomic_feature)
  tmpdir contains: sites.bed, CDS.bed, UTR.bed, intermediate files

Requirements:
  - Python packages: pandas
  - bedtools in PATH
"""
import argparse, gzip, os, sys, tempfile, subprocess, re
from pathlib import Path
import pandas as pd
from collections import defaultdict

def detect_gtf_chr_style(gtf_path):
    # read first non-comment line and check seqname prefix
    opener = gzip.open if str(gtf_path).endswith('.gz') else open
    with opener(gtf_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            if ln.startswith('#') or ln.strip()=="":
                continue
            seq = ln.split('\t',1)[0]
            return seq.startswith('chr')
    return False

def parse_pos_field_raw(s):
    """Return (chrom, pos_int, strand) or (None,None,None)"""
    if s is None:
        return (None,None,None)
    s = str(s).strip()
    if s == "":
        return (None,None,None)

    # Robust parsing by working from the right.
    # Supports: chr1_12345_+, chr1_12345_-, chr1:12345:+, chr1-12345-+, and contigs with '_' in name.
    strand = None
    m = re.search(r'[_:\-]([+-])$', s)
    if m:
        strand = m.group(1)
        s = s[:-2]

    m = re.match(r'^(.*)[_:\-](\d+)$', s)
    if m:
        chrom = m.group(1)
        pos = int(m.group(2))
        return (chrom, pos, strand)

    # fallback: pattern chr1:12345(+)
    m = re.match(r'([^:]+):(\d+).*?([+-])?', s)
    if m:
        chrom = m.group(1)
        pos = int(m.group(2))
        strand = m.group(3) if m.lastindex>=3 else None
        return (chrom, pos, strand)
    return (None,None,None)

def normalize_chrom(chrom, want_chr):
    if chrom is None:
        return None
    if want_chr:
        return chrom if chrom.startswith('chr') else 'chr' + chrom
    else:
        return chrom[3:] if chrom.startswith('chr') else chrom

def write_sites_bed_from_csv(csv_path, sites_bed_path, want_chr):
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    first_col = df.columns[0]
    with open(sites_bed_path, 'w') as out:
        for i, val in enumerate(df[first_col].tolist()):
            chrom,pos,strand = parse_pos_field_raw(val)
            if chrom is None or pos is None:
                continue
            chrom = normalize_chrom(chrom, want_chr)
            start = pos - 1
            end = pos
            strand_col = strand if strand in ('+','-') else '.'
            # name column we set to original 0-based row index (string)
            out.write(f"{chrom}\t{start}\t{end}\t{i}\t.\t{strand_col}\n")
    # sanity: lines count
    return

def build_cds_and_utr_from_gtf(gtf_path, cds_bed, utr_bed):
    """
    Build CDS.bed and UTR.bed (UTR = exon - CDS per transcript if explicit UTRs not present)
    Writes 6-col bed: chrom, start, end, transcript_id, ., strand
    """
    transcripts_exons = defaultdict(list)
    transcripts_cds = defaultdict(list)
    transcripts_explicit_utr = defaultdict(list)
    has_explicit_utr = False

    opener = gzip.open if str(gtf_path).endswith('.gz') else open
    with opener(gtf_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            if ln.startswith('#') or ln.strip()=="":
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqname, src, feature, start, end, score, strand, phase, attrs = cols[:9]
            try:
                s = int(start) - 1
                e = int(end)
            except:
                continue
            # parse transcript_id from attributes
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            tid = m.group(1) if m else None
            key = (seqname, tid, strand)
            if feature == 'exon':
                transcripts_exons[key].append((s,e))
            elif feature == 'CDS':
                transcripts_cds[key].append((s,e))
            elif feature in ('five_prime_UTR','three_prime_UTR','UTR'):
                has_explicit_utr = True
                transcripts_explicit_utr[key].append((s,e,feature))

    def merge(intervals):
        if not intervals:
            return []
        intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
        out = []
        cur_s, cur_e = intervals[0]
        for s,e in intervals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                out.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((cur_s, cur_e))
        return out

    def subtract(a_list, b_list):
        a = merge(a_list)
        b = merge(b_list)
        res = []
        bi = 0
        for as_, ae in a:
            cur = as_
            while bi < len(b) and b[bi][1] <= cur:
                bi += 1
            tmp = bi
            while tmp < len(b) and b[tmp][0] < ae:
                bs, be = b[tmp]
                if bs > cur:
                    res.append((cur, min(bs, ae)))
                cur = max(cur, be)
                if cur >= ae:
                    break
                tmp += 1
            if cur < ae:
                res.append((cur, ae))
        return res

    # write CDS and UTR
    with open(cds_bed, 'w') as fc, open(utr_bed, 'w') as fu:
        # explicit UTRs if present
        if has_explicit_utr:
            for key, ivs in transcripts_explicit_utr.items():
                chrom, tid, strand = key
                for s,e,feat in ivs:
                    fu.write(f"{chrom}\t{s}\t{e}\t{tid};{feat}\t.\t{strand}\n")
        # write CDS merged
        for key, ivs in transcripts_cds.items():
            chrom, tid, strand = key
            for s,e in merge(ivs):
                fc.write(f"{chrom}\t{s}\t{e}\t{tid}\t.\t{strand}\n")
        # if no explicit UTRs, compute per transcript exon-CDS
        if not has_explicit_utr:
            for key, exons in transcripts_exons.items():
                chrom, tid, strand = key
                exm = merge(exons)
                cdsm = merge(transcripts_cds.get(key, []))
                if not cdsm:
                    continue
                utrs = subtract(exm, cdsm)
                # determine min/max CDS to classify 5' vs 3' (not strictly needed here)
                min_cds = min([s for s,e in cdsm])
                max_cds = max([e for s,e in cdsm])
                for s,e in utrs:
                    # label not strictly required; we include no label or include tid
                    fu.write(f"{chrom}\t{s}\t{e}\t{tid}\t.\t{strand}\n")
    return

def run_bedtools_intersect_sites(sites_bed, cds_bed, utr_bed, tmpdir):
    a_cds = os.path.join(tmpdir, "sites_vs_cds.tsv")
    a_utr = os.path.join(tmpdir, "sites_vs_utr.tsv")
    cmd1 = f"bedtools intersect -s -a {sites_bed} -b {cds_bed} -wa -wb -loj > {a_cds}"
    cmd2 = f"bedtools intersect -s -a {sites_bed} -b {utr_bed} -wa -wb -loj > {a_utr}"
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    return a_cds, a_utr

def parse_intersect_to_map(path):
    d = {}
    with open(path) as fh:
        for ln in fh:
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            site_idx = cols[3]
            # feature attr from feature bed usually at col9 (0-based indexing after -wa -wb)
            feat_attr = cols[9] if len(cols) > 9 else ''
            if site_idx not in d:
                d[site_idx] = []
            if feat_attr and feat_attr != '.':
                d[site_idx].append(feat_attr)
    return d

def annotate_csv_with_maps(csv_path, sites_bed, cds_map, utr_map, out_path):
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    ann = []
    for i in range(len(df)):
        key = str(i)
        if key in cds_map and len(cds_map[key])>0:
            ann.append('CDS')
        elif key in utr_map and len(utr_map[key])>0:
            ann.append('UTR')
        else:
            ann.append('other')
    df['genomic_feature'] = ann
    df.to_csv(out_path, index=False)
    return

def verify_annotation_consistency(csv_path, sites_bed, cds_bed, utr_bed, annotated_csv, tmpdir):
    # produce ground truth sets using bedtools -u (unique)
    in_cds = os.path.join(tmpdir, "sites.in.CDS.bed")
    in_utr = os.path.join(tmpdir, "sites.in.UTR.bed")
    subprocess.check_call(f"bedtools intersect -s -a {sites_bed} -b {cds_bed} -wa -u > {in_cds}", shell=True)
    subprocess.check_call(f"bedtools intersect -s -a {sites_bed} -b {utr_bed} -wa -u > {in_utr}", shell=True)
    # build sets of indices
    cds_idx = set()
    utr_idx = set()
    for f,p in [(in_cds, cds_idx), (in_utr, utr_idx)]:
        with open(f) as fh:
            for ln in fh:
                cols = ln.split('\t')
                if len(cols) >= 4:
                    p.add(cols[3])
    # load annotated csv
    df = pd.read_csv(annotated_csv, dtype=str, keep_default_na=False)
    mismatches = []
    total = len(df)
    match = 0
    for i,row in df.iterrows():
        annotated = row.get('genomic_feature','').strip()
        gt = 'other'
        key = str(i)
        if key in cds_idx:
            gt = 'CDS'
        elif key in utr_idx:
            gt = 'UTR'
        if annotated == gt:
            match += 1
        else:
            pos = row.iloc[0]
            mismatches.append((i, pos, annotated, gt))
    # write mismatches
    mm_path = os.path.join(tmpdir, "mismatches.tsv")
    with open(mm_path, 'w') as out:
        out.write("index\tpos\tannotated\tground_truth\n")
        for a in mismatches:
            out.write("\t".join(map(str,a)) + "\n")
    print(f"VERIFY: total={total} match={match} mismatches={len(mismatches)} -> mismatches saved to {mm_path}")
    return mm_path, total, match, len(mismatches)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gtf', required=True)
    p.add_argument('--output-dir', default=None, help='directory for annotated CSVs (default: alongside input)')
    p.add_argument('--verify', action='store_true', help='run verification after annotation')
    p.add_argument('csvs', nargs='+')
    args = p.parse_args()
    if args.output_dir:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # check bedtools
    try:
        subprocess.check_call(['bedtools','--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("ERROR: bedtools not found in PATH. Install: conda install -c bioconda bedtools")
        sys.exit(1)

    tmpdir = tempfile.mkdtemp(prefix='annot_v3_')
    print("tmpdir:", tmpdir)

    # detect gtf chr style
    gtf_has_chr = detect_gtf_chr_style(args.gtf)
    print("GTF uses 'chr' prefix?" , gtf_has_chr)

    # build cds & utr
    cds_bed = os.path.join(tmpdir, "CDS.bed")
    utr_bed = os.path.join(tmpdir, "UTR.bed")
    print("Building CDS/UTR beds...")
    build_cds_and_utr_from_gtf(args.gtf, cds_bed, utr_bed)
    print("CDS bed:", cds_bed, "UTR bed:", utr_bed)

    for csvf in args.csvs:
        print("Processing", csvf)
        sites_bed = os.path.join(tmpdir, Path(csvf).stem + ".sites.bed")
        write_sites_bed_from_csv(csvf, sites_bed, want_chr=gtf_has_chr)
        # sanity check: line counts
        n_sites = sum(1 for _ in open(sites_bed))
        n_csv = sum(1 for _ in open(csvf)) - 1
        print("rows:", n_csv, "sites.bed lines:", n_sites)
        # intersect
        a_cds, a_utr = run_bedtools_intersect_sites(sites_bed, cds_bed, utr_bed, tmpdir)
        cds_map = parse_intersect_to_map(a_cds)
        utr_map = parse_intersect_to_map(a_utr)
        out_csv = str((Path(args.output_dir) if args.output_dir else Path(csvf).parent) / (Path(csvf).stem + '.annotated_v3.csv'))
        annotate_csv_with_maps(csvf, sites_bed, cds_map, utr_map, out_csv)
        print("Wrote annotated file:", out_csv)

        if args.verify:
            print("Running verification...")
            mm_path, total, match, mismatches = verify_annotation_consistency(csvf, sites_bed, cds_bed, utr_bed, out_csv, tmpdir)
            print(f"Verification summary: total={total}, match={match}, mismatches={mismatches}")
            print("If mismatches > 0, inspect", mm_path)

    print("ALL DONE. tmpdir retained at", tmpdir)

if __name__ == '__main__':
    main()
