#!/usr/bin/env python3
"""Assign coding consequences of RNA-strand A-to-G edits at m6A sites.

Input is a CSV whose first column, pos, contains IDs such as
chr10_11462944_-. For each matching CDS transcript, the script reports
whether changing the transcript's A to G changes the encoded amino acid.

The GTF/input CSV use UCSC chromosome names (for example, chr10), while
the indexed GRCh38 FASTA may use RefSeq accessions (NC_000010.11).
Chromosome names are resolved against the FASTA index before sequence access.
Negative-strand genomic bases are complemented after ordering coordinates in
transcript 5'-to-3' direction; they are not reversed a second time.

Requires Python 3, pandas, pysam, and an indexed reference FASTA.

Example:
  python scripts/predict_AtoG_coding_effects.py \
    --gtf hg38.ncbiRefSeq.gtf.gz \
    --fasta GRCh38.fa \
    --out HEK_HeLa_over80_intersection.AtoG.effects.csv \
    HEK_HeLa_over80_intersection.csv
"""

import argparse
import gzip
import re
import sys
import time
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pysam

# Genetic code table
GCODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def revcomp(seq: str) -> str:
    comp = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))
    return comp[::-1]


def genomic_bases_to_mrna_codon(bases_in_transcript_order: str, strand: str) -> str:
    """Convert genome bases already ordered 5' to 3' along the transcript."""
    if strand == '+':
        return bases_in_transcript_order.upper()
    if strand == '-':
        return bases_in_transcript_order.upper().translate(str.maketrans('ACGT', 'TGCA'))
    raise ValueError(f'Unexpected strand: {strand}')


def load_fasta_contigs(fasta_path: str) -> List[str]:
    """Load contig names from a FASTA index (.fai)."""
    fai = Path(str(fasta_path) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(f"FASTA index not found: {fai}. Run: samtools faidx {fasta_path}")

    contigs: List[str] = []
    with open(fai, 'rt', encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            if not ln.strip():
                continue
            contigs.append(ln.split('\t', 1)[0])
    return contigs


def build_refseq_base_index(contigs: List[str]) -> Dict[str, str]:
    """Index RefSeq contigs by base accession (e.g., NC_000010 -> NC_000010.11)."""
    idx: Dict[str, str] = {}
    for name in contigs:
        if name.startswith(('NC_', 'NT_', 'NW_')) and '.' in name:
            base = name.split('.', 1)[0]
            idx.setdefault(base, name)
    return idx


def resolve_fasta_contig(query_chrom: str,
                         fasta_contigs_set: set,
                         refseq_base_index: Dict[str, str]) -> Optional[str]:
    """Resolve an input chrom name to a contig present in FASTA.

    Supports:
      - exact match
      - add/remove 'chr'
      - GRCh38 UCSC chrN/chrX/chrY/chrM -> RefSeq NC_0000xx.xx / NC_012920.1

    Returns the FASTA contig name, or None if no reasonable mapping exists.
    """
    if query_chrom is None:
        return None
    chrom = str(query_chrom).strip()
    if chrom == '':
        return None

    if chrom in fasta_contigs_set:
        return chrom

    chrom_upper = chrom.upper()
    if chrom_upper.startswith('CHR'):
        no_chr = chrom[3:]
        if no_chr in fasta_contigs_set:
            return no_chr
    else:
        with_chr = 'chr' + chrom
        if with_chr in fasta_contigs_set:
            return with_chr

    # RefSeq mapping (common for GRCh38 FASTA like the one in your workspace)
    c = chrom_upper[3:] if chrom_upper.startswith('CHR') else chrom_upper
    if c == 'MT':
        c = 'M'

    if c.isdigit():
        n = int(c)
        if 1 <= n <= 22:
            return refseq_base_index.get(f"NC_{n:06d}")
        if n == 23:
            return refseq_base_index.get('NC_000023')
        if n == 24:
            return refseq_base_index.get('NC_000024')

    if c == 'X':
        return refseq_base_index.get('NC_000023')
    if c == 'Y':
        return refseq_base_index.get('NC_000024')
    if c in ('M', 'MITO'):
        return refseq_base_index.get('NC_012920')

    return None


# --- GTF parsing: build transcripts -> merged CDS intervals (transcript order) ---

CDSIntervals = List[Tuple[int, int]]
TranscriptKey = Tuple[str, str, str]  # (chrom, tid, strand)
TranscriptEntry = Tuple[str, str, CDSIntervals]  # (tid, strand, ivs)


def parse_gtf_build_cds(gtf_path: str) -> Dict[TranscriptKey, CDSIntervals]:
    # returns dict: key=(chrom,tid,strand) -> list of (s0,e) intervals merged
    transcripts_cds: Dict[TranscriptKey, CDSIntervals] = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith('.gz') else open
    with opener(gtf_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            if ln.startswith('#') or ln.strip() == '':
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqname, src, feature, start, end, score, strand, phase, attrs = cols[:9]
            if feature != 'CDS':
                continue

            m = re.search(r'transcript_id "([^"]+)"', attrs)
            tid = m.group(1) if m else None
            if tid is None:
                m2 = re.search(r'gene_id "([^"]+)"', attrs)
                tid = m2.group(1) if m2 else "no_tid"

            s0 = int(start) - 1
            e = int(end)
            transcripts_cds[(seqname, tid, strand)].append((s0, e))

    merged: Dict[TranscriptKey, CDSIntervals] = {}
    for key, ivs in transcripts_cds.items():
        ivs_sorted = sorted(ivs, key=lambda x: (x[0], x[1]))
        if not ivs_sorted:
            continue
        out: CDSIntervals = []
        cs, ce = ivs_sorted[0]
        for s, e in ivs_sorted[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                out.append((cs, ce))
                cs, ce = s, e
        out.append((cs, ce))
        merged[key] = out
    return merged


def build_cds_bin_index(merged_cds: Dict[TranscriptKey, CDSIntervals], bin_size: int = 1_000_000):
    """Build a simple genomic bin index to avoid scanning all transcripts per site.

    index[chrom][bin_id] -> list of (tid, strand, ivs)

    This dramatically speeds up queries like: find transcripts whose CDS contains chrom:pos.
    """
    index: Dict[str, Dict[int, List[TranscriptEntry]]] = defaultdict(lambda: defaultdict(list))

    for (chrom, tid, strand), ivs in merged_cds.items():
        # Add transcript to every bin overlapped by any CDS interval.
        # Duplicates across intervals/bins are OK; we deduplicate per query.
        for s0, e in ivs:
            if e <= s0:
                continue
            b0 = s0 // bin_size
            b1 = (e - 1) // bin_size
            for b in range(b0, b1 + 1):
                index[chrom][b].append((tid, strand, ivs))

    return index


def find_transcripts_containing_pos(index, chrom: str, pos1: int, bin_size: int = 1_000_000) -> List[TranscriptEntry]:
    """Return list of (tid, strand, ivs) where pos1 is inside any CDS interval."""
    if chrom is None:
        return []
    pos0 = pos1 - 1
    b = pos0 // bin_size
    candidates = index.get(chrom, {}).get(b, [])

    res: List[TranscriptEntry] = []
    seen = set()
    for tid, strand, ivs in candidates:
        key = (tid, strand)
        if key in seen:
            continue
        seen.add(key)
        for s, e in ivs:
            if s <= pos0 < e:
                res.append((tid, strand, ivs))
                break
    return res


# --- compute transcript coding offset (0-based) for a given pos (1-based) ---

def compute_coding_offset(ivs: CDSIntervals, pos1: int, strand: str) -> Optional[int]:
    pos0 = pos1 - 1
    if strand == '+':
        order = sorted(ivs, key=lambda x: x[0])
    else:
        order = sorted(ivs, key=lambda x: x[0], reverse=True)

    cum = 0
    for s, e in order:
        if s <= pos0 < e:
            if strand == '+':
                offset_in_exon = pos0 - s
            else:
                offset_in_exon = (e - 1) - pos0
            return cum + offset_in_exon
        cum += (e - s)
    return None


def offsets_to_genomic_coords(ivs: CDSIntervals, strand: str, start_offset: int, L: int) -> Optional[List[int]]:
    if strand == '+':
        order = sorted(ivs, key=lambda x: x[0])
    else:
        order = sorted(ivs, key=lambda x: x[0], reverse=True)

    coords: List[int] = []
    cum = 0
    needed = L
    cur_offset = start_offset

    for s, e in order:
        seglen = e - s
        if cur_offset >= cum + seglen:
            cum += seglen
            continue

        start_in_exon = cur_offset - cum
        if start_in_exon < 0:
            start_in_exon = 0

        if strand == '+':
            gpos0 = s + start_in_exon
            avail = seglen - start_in_exon
            take = min(avail, needed)
            for i in range(take):
                coords.append(gpos0 + i + 1)
        else:
            gpos0 = (e - 1) - start_in_exon
            avail = seglen - start_in_exon
            take = min(avail, needed)
            for i in range(take):
                coords.append(gpos0 - i + 1)

        needed -= take
        cur_offset += take

        if needed == 0:
            break
        cum += seglen

    if needed != 0:
        return None
    return coords

@lru_cache(maxsize=4)
def _open_fasta(fasta: str):
    return pysam.FastaFile(fasta)


@lru_cache(maxsize=500000)
def fetch_base_indexed(fasta: str, contig: Optional[str], pos1: int) -> Optional[str]:
    """Fetch one 1-based genome base from indexed FASTA using pysam."""
    if contig is None:
        return None
    try:
        seq = _open_fasta(fasta).fetch(contig, pos1 - 1, pos1).upper()
    except (KeyError, ValueError):
        return None
    return seq[0] if seq else None

class ProgressPrinter:
    def __init__(self, total: int, every: int = 2000):
        self.total = max(int(total), 0)
        self.every = max(int(every), 1)
        self.t0 = time.time()
        self.last_t = self.t0

    def update(self, i: int, analyzed: int, mRNA_A: int, nonsyn: int):
        if i == 0 or (i % self.every != 0 and i != self.total):
            return
        now = time.time()
        elapsed = now - self.t0
        rate = i / elapsed if elapsed > 0 else 0.0
        eta = (self.total - i) / rate if rate > 0 else float('inf')
        eta_str = f"{eta:,.0f}s" if eta != float('inf') else "?"
        msg = (
            f"[progress] {i}/{self.total} rows | analyzed={analyzed} | mRNA_A={mRNA_A} | nonsyn={nonsyn} "
            f"| elapsed={elapsed:,.0f}s | {rate:,.1f} rows/s | ETA={eta_str}"
        )
        print(msg, file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gtf', required=True, help='GTF (gtf or gtf.gz)')
    p.add_argument('--fasta', required=True, help='reference fasta (indexed with samtools faidx)')
    p.add_argument('--out', required=True, help='output csv')
    p.add_argument('--bin-size', type=int, default=1_000_000, help='genomic bin size for CDS indexing (default: 1,000,000)')
    p.add_argument('--progress-every', type=int, default=2000, help='print progress every N rows (default: 2000)')
    p.add_argument('--verbose', action='store_true', help='print extra info (chrom mapping warnings, stage timing)')
    p.add_argument('--assume-cds', action='store_true', help="assume input sites are already CDS; do not label rows as not_in_CDS (use no_transcript_match if needed)")
    p.add_argument('--unique-per-site', action='store_true', help="output at most one result row per input site (choose a representative transcript)")
    p.add_argument('csv', help='input csv with pos in first column (chr_pos_strand)')
    args = p.parse_args()

    _open_fasta(args.fasta)

    t_start = time.time()

    # Stage 1: GTF parse
    print("[stage] Parsing GTF and building CDS per transcript (this may take a bit)...", file=sys.stderr)
    t0 = time.time()
    merged_cds = parse_gtf_build_cds(args.gtf)
    t1 = time.time()
    print(f"[stage] Parsed {len(merged_cds)} transcripts with CDS entries. ({t1 - t0:,.1f}s)", file=sys.stderr)

    # Stage 2: build fast query index
    print(f"[stage] Building CDS bin index (bin_size={args.bin_size})...", file=sys.stderr)
    t0 = time.time()
    cds_index = build_cds_bin_index(merged_cds, bin_size=args.bin_size)
    t1 = time.time()
    if args.verbose:
        n_chrom = len(cds_index)
        n_bins = sum(len(v) for v in cds_index.values())
        print(f"[stage] CDS bin index ready: {n_chrom} chroms, {n_bins} bins. ({t1 - t0:,.1f}s)", file=sys.stderr)
    else:
        print(f"[stage] CDS bin index ready. ({t1 - t0:,.1f}s)", file=sys.stderr)

    # Stage 3: FASTA contigs
    print("[stage] Loading FASTA contigs (.fai) for contig-name mapping...", file=sys.stderr)
    t0 = time.time()
    fasta_contigs = load_fasta_contigs(args.fasta)
    fasta_contigs_set = set(fasta_contigs)
    refseq_base_index = build_refseq_base_index(fasta_contigs)
    t1 = time.time()
    print(f"[stage] Loaded {len(fasta_contigs)} FASTA contigs. ({t1 - t0:,.1f}s)", file=sys.stderr)

    # Stage 4: read CSV
    df = pd.read_csv(args.csv, dtype=str, keep_default_na=False)
    pos_col = df.columns[0]

    out_rows = []
    summary = {'total_sites': 0, 'analyzed_sites': 0, 'mRNA_A_sites': 0, 'nonsynonymous_hits': 0, 'no_transcript_match': 0}

    progress = ProgressPrinter(total=len(df), every=args.progress_every)
    warned_unmapped = set()

    print(f"[stage] Processing CSV rows: {len(df)}", file=sys.stderr)
    if args.assume_cds:
        print("[info] assume input is CDS; will not output not_in_CDS (use no_transcript_match when no transcript overlaps)", file=sys.stderr)

    for i, (_, row) in enumerate(df.iterrows(), start=1):
        summary['total_sites'] += 1
        posfield = str(row[pos_col]).strip()
        rna_strand = posfield[-1] if posfield.endswith(('+', '-')) else None

        parts = re.split(r'[_:\-]', posfield)
        if len(parts) < 2 or rna_strand is None:
            out_rows.append({**row.to_dict(), **{
                'analysis_transcript': '', 'analysis_result': 'pos_parse_failed'
            }})
            progress.update(i, summary['analyzed_sites'], summary['mRNA_A_sites'], summary['nonsynonymous_hits'])
            continue

        chrom = parts[0]
        try:
            pos1 = int(re.sub(r'[^\d]', '', parts[1]))
        except Exception:
            out_rows.append({**row.to_dict(), **{
                'analysis_transcript': '', 'analysis_result': 'pos_parse_failed'
            }})
            progress.update(i, summary['analyzed_sites'], summary['mRNA_A_sites'], summary['nonsynonymous_hits'])
            continue

        fasta_contig = resolve_fasta_contig(chrom, fasta_contigs_set, refseq_base_index)
        if fasta_contig is None and args.verbose and chrom not in warned_unmapped:
            warned_unmapped.add(chrom)
            print(f"[warn] Cannot map chrom '{chrom}' to a FASTA contig name; bases will be 'N'.", file=sys.stderr)

        tlist = [hit for hit in find_transcripts_containing_pos(cds_index, chrom, pos1, bin_size=args.bin_size)
                 if hit[1] == rna_strand]
        if not tlist:
            if args.assume_cds:
                summary['no_transcript_match'] += 1
                result = 'no_transcript_match'
            else:
                result = 'not_in_CDS'
            out_rows.append({**row.to_dict(), **{
                'analysis_transcript': '', 'analysis_result': result
            }})
            progress.update(i, summary['analyzed_sites'], summary['mRNA_A_sites'], summary['nonsynonymous_hits'])
            continue

        summary['analyzed_sites'] += 1
        n_transcripts_hit = len(tlist)
        site_recs = []

        for (tid, strand, ivs) in tlist:
            cod_offset = compute_coding_offset(ivs, pos1, strand)
            if cod_offset is None:
                site_recs.append({**row.to_dict(), **{
                    'analysis_transcript': tid,
                    'transcript_strand': strand,
                    'analysis_result': 'offset_fail',
                    'n_transcripts_hit': n_transcripts_hit,
                }})
                continue

            codon_index = cod_offset // 3
            codon_start_offset = codon_index * 3

            coords = offsets_to_genomic_coords(ivs, strand, codon_start_offset, 3)
            if coords is None or len(coords) != 3:
                site_recs.append({**row.to_dict(), **{
                    'analysis_transcript': tid,
                    'transcript_strand': strand,
                    'analysis_result': 'codon_coords_fail',
                    'n_transcripts_hit': n_transcripts_hit,
                }})
                continue

            bases = []
            for p1 in coords:
                b = fetch_base_indexed(args.fasta, fasta_contig, p1)
                bases.append(b if b is not None else 'N')

            genomic_codon = ''.join(bases).upper()
            mrna_codon_for_translate = genomic_bases_to_mrna_codon(genomic_codon, strand)

            pos_in_codon = cod_offset - codon_start_offset
            current_base = mrna_codon_for_translate[pos_in_codon]
            is_mrna_A = (current_base.upper() == 'A')

            mutated_codon = list(mrna_codon_for_translate)
            mutated_codon[pos_in_codon] = 'G'
            mutated_codon = ''.join(mutated_codon)

            orig_aa = GCODE.get(mrna_codon_for_translate.upper(), 'X')
            mut_aa = GCODE.get(mutated_codon.upper(), 'X')
            effect = 'synonymous' if orig_aa == mut_aa else 'nonsynonymous'

            genomic_ref_b = fetch_base_indexed(args.fasta, fasta_contig, pos1)
            genomic_ref = genomic_ref_b.upper() if genomic_ref_b else 'N'
            if strand == '+':
                genomic_change = f"{genomic_ref}->G"
            else:
                genomic_change = f"{genomic_ref}->C"

            rec = {
                **row.to_dict(),
                'analysis_transcript': tid,
                'transcript_strand': strand,
                'analysis_result': 'ok',
                'n_transcripts_hit': n_transcripts_hit,
                'codon_genomic_coords_1based': ';'.join(map(str, coords)),
                'genomic_codon': genomic_codon,
                'mrna_codon': mrna_codon_for_translate,
                'pos_in_codon': pos_in_codon,
                'mRNA_base': current_base,
                'is_mRNA_A': is_mrna_A,
                'genomic_change_equiv': genomic_change,
                'orig_AA': orig_aa,
                'mut_AA': mut_aa,
                'effect': effect,
            }
            site_recs.append(rec)

            if is_mrna_A:
                summary['mRNA_A_sites'] += 1
                if effect == 'nonsynonymous':
                    summary['nonsynonymous_hits'] += 1

        if args.unique_per_site:
            # Choose a representative transcript per site.
            # Preference: successful records first, then RefSeq NM_, then lexicographic transcript id.
            def _rank(r):
                tid2 = str(r.get('analysis_transcript', ''))
                ok = (r.get('analysis_result') == 'ok')
                if tid2.startswith('NM_'):
                    cls = 0
                elif tid2.startswith('NR_'):
                    cls = 1
                elif tid2.startswith('XM_'):
                    cls = 2
                else:
                    cls = 3
                return (0 if ok else 1, cls, tid2)

            out_rows.append(sorted(site_recs, key=_rank)[0])
        else:
            out_rows.extend(site_recs)

        progress.update(i, summary['analyzed_sites'], summary['mRNA_A_sites'], summary['nonsynonymous_hits'])

    outdf = pd.DataFrame(out_rows)
    outdf.to_csv(args.out, index=False)

    t_end = time.time()
    print("[done] Wrote", args.out, file=sys.stderr)
    print("[done] SUMMARY:", summary, file=sys.stderr)
    if args.verbose:
        print(f"[done] Total elapsed: {t_end - t_start:,.1f}s", file=sys.stderr)


if __name__ == '__main__':
    main()
