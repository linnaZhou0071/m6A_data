#!/usr/bin/env python3
"""Historical coding-effect implementation, retained only for provenance.

This superseded script predicts the amino-acid consequence of RNA A-to-G
changes at CDS sites supplied in a CSV with a pos column. Its negative-strand
handling is not the corrected analysis. Use scripts/predict_AtoG_coding_effects.py
for current results.

Requires Python 3, pandas, bedtools and samtools on PATH.

Historical invocation:
  python scripts/legacy/predict_AtoG_coding_effects.py \
    --gtf hg38.ncbiRefSeq.gtf.gz \
    --fasta hg38.fa \
    HEK_HeLa_over80_intersection.csv \
    --out intersection.AtoG.effects.csv
"""
import argparse, gzip, re, subprocess, sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Genetic code table
GCODE = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
 'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
 'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def revcomp(seq):
    comp = seq.translate(str.maketrans("ACGTacgt","TGCAtgca"))
    return comp[::-1]

def run_cmd(cmd):
    r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\nstderr: {r.stderr}")
    return r.stdout

# --- GTF parsing: build transcripts -> merged CDS intervals (transcript order) ---
def parse_gtf_build_cds(gtf_path):
    # returns dict: key=(chrom,tid,strand) -> list of (s0,e) intervals merged
    transcripts_cds = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith('.gz') else open
    with opener(gtf_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            if ln.startswith('#') or ln.strip()=='':
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqname, src, feature, start, end, score, strand, phase, attrs = cols[:9]
            if feature != 'CDS':
                continue
            # transcript id
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            tid = m.group(1) if m else None
            if tid is None:
                # fallback to gene_id or skip
                m2 = re.search(r'gene_id "([^"]+)"', attrs)
                tid = m2.group(1) if m2 else "no_tid"
            s0 = int(start) - 1
            e = int(end)
            transcripts_cds[(seqname, tid, strand)].append((s0,e))
    # merge and sort per transcript in transcript order (for '-' we keep intervals but will reverse order later)
    merged = {}
    for key, ivs in transcripts_cds.items():
        ivs_sorted = sorted(ivs, key=lambda x:(x[0], x[1]))
        out=[]
        cs,ce = ivs_sorted[0]
        for s,e in ivs_sorted[1:]:
            if s <= ce:
                ce = max(ce,e)
            else:
                out.append((cs,ce))
                cs,ce = s,e
        out.append((cs,ce))
        merged[key] = out
    return merged

def parse_pos_field(s):
    if pd.isna(s):
        return (None, None, None)
    s = str(s).strip()
    if s == '':
        return (None, None, None)
    
    # Split the position ID into chromosome, position and strand.
    for sep in ['_', ':', '\t', ' ' , '-']:
        if sep in s:
            parts = s.split(sep)
            break
    else:
        parts = [s]
    
    chrom = None
    pos = None
    strand = None
    
    if len(parts) >= 2:
        chrom = parts[0]
        p = re.sub(r'[^\d]', '', parts[1])
        if p.isdigit():
            pos = int(p)
        if len(parts) >= 3:
            s2 = parts[2].strip()
            if s2 in ['+', '-']:
                strand = s2
    else:
        m = re.match(r'([^:]+):(\d+)(?:[:(]?([+-])', s)
        if m:
            chrom = m.group(1)
            pos = int(m.group(2))
            strand = m.group(3) if m.lastindex >= 3 else None
    
    if chrom is None:
        return (None, None, None)
    
    # Remove the UCSC chromosome prefix.
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    
    chrom = chrom.upper()
    return (chrom, pos, strand)

# --- helper to check if pos (1-based) in any CDS interval and list matching transcripts ---
def find_transcripts_containing_pos(merged_cds, chrom, pos1):
    # pos1 is 1-based genomic coordinate
    res = []
    pos0 = pos1 - 1
    for (c, tid, strand), ivs in merged_cds.items():
        if c != chrom:
            continue
        for s,e in ivs:
            if s <= pos0 < e:
                res.append((tid, strand, ivs))
                break
    return res

# --- compute transcript coding offset (0-based) for a given pos (1-based) in a given transcript CDS intervals ---
def compute_coding_offset(ivs, pos1, strand):
    # ivs: merged CDS intervals for transcript in genomic coordinates (s0,e)
    # strand: '+' or '-'
    pos0 = pos1 - 1
    # order intervals in transcript order
    if strand == '+':
        order = sorted(ivs, key=lambda x:x[0])
    else:
        order = sorted(ivs, key=lambda x:x[0], reverse=True)
    cum = 0
    for s,e in order:
        if s <= pos0 < e:
            if strand == '+':
                offset_in_exon = pos0 - s
            else:
                # for minus strand, transcript first base of this exon corresponds to e-1
                offset_in_exon = (e - 1) - pos0
            return cum + offset_in_exon
        else:
            cum += (e - s)
    return None

# --- map a transcript offset (0-based) and length L to list of genomic 1-based coords (may cross exons) ---
def offsets_to_genomic_coords(ivs, strand, start_offset, L):
    # ivs: merged CDS intervals as above; strand determines transcript order
    if strand == '+':
        order = sorted(ivs, key=lambda x:x[0])
    else:
        order = sorted(ivs, key=lambda x:x[0], reverse=True)
    coords = []
    cum = 0
    needed = L
    for s,e in order:
        seglen = e - s
        if start_offset >= cum + seglen:
            cum += seglen
            continue
        # starting inside this exon
        start_in_exon = start_offset - cum
        # produce bases from start_in_exon to end of exon (in transcript order)
        if strand == '+':
            gpos0 = s + start_in_exon  # 0-based
            # how many bases available in this exon from here
            avail = seglen - start_in_exon
            take = min(avail, needed)
            for i in range(take):
                coords.append(gpos0 + i + 1)  # convert to 1-based
            needed -= take
        else:
            # minus strand: transcript order uses bases e-1, e-2, ...
            gpos0 = (e - 1) - start_in_exon
            avail = seglen - start_in_exon
            take = min(avail, needed)
            for i in range(take):
                coords.append(gpos0 - i + 1)  # 1-based
            needed -= take
        if needed == 0:
            break
        cum += seglen
    if needed != 0:
        return None
    return coords

# --- fetch single base from fasta using samtools faidx ---
def fetch_base_samtools(fasta, chrom, pos1):
    # returns uppercase base char
    region = f"{chrom}:{pos1}-{pos1}"
    out = run_cmd(f"samtools faidx {fasta} {region}")
    lines = [l.strip() for l in out.splitlines() if l.strip() and not l.startswith('>')]
    seq = ''.join(lines).upper()
    return seq[0] if seq else None

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gtf', required=True, help='GTF (gtf or gtf.gz)')
    p.add_argument('--fasta', required=True, help='reference fasta (indexed with samtools faidx)')
    p.add_argument('--out', required=True, help='output csv')
    p.add_argument('csv', help='input csv with pos in first column (chr_pos_strand)')
    args = p.parse_args()

    # check tools
    try:
        subprocess.run(['bedtools','--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("ERROR: bedtools not found in PATH", file=sys.stderr)
        sys.exit(1)
    try:
        subprocess.run(['samtools','faidx','-h'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("ERROR: samtools not found in PATH", file=sys.stderr)
        sys.exit(1)

    # parse GTF -> CDS per transcript
    print("Parsing GTF and building CDS per transcript (this may take a bit)...")
    merged_cds = parse_gtf_build_cds(args.gtf)
    print(f"Parsed {len(merged_cds)} transcripts with CDS entries.")

    # read CSV
    df = pd.read_csv(args.csv, dtype=str, keep_default_na=False)
    pos_col = df.columns[0]

    out_rows = []
    summary = {'total_sites':0, 'analyzed_sites':0, 'mRNA_A_sites':0, 'nonsynonymous_hits':0}

    for idx, row in df.iterrows():
        summary['total_sites'] += 1
        posfield = str(row[pos_col]).strip()
        # parse pos
        parts = re.split(r'[_:\-]', posfield)
        if len(parts) < 2:
            # skip
            out_rows.append({**row.to_dict(), **{
                'analysis_transcript':'', 'analysis_result':'pos_parse_failed'
            }})
            continue
        chrom = parts[0]
        pos1 = int(re.sub(r'[^\d]','', parts[1]))
        # find transcripts containing the genomic pos in CDS
        tlist = find_transcripts_containing_pos(merged_cds, chrom, pos1)
        if not tlist:
            out_rows.append({**row.to_dict(), **{
                'analysis_transcript':'', 'analysis_result':'not_in_CDS'
            }})
            continue
        summary['analyzed_sites'] += 1
        # process each transcript
        for (tid, strand, ivs) in tlist:
            # compute coding offset in transcript
            cod_offset = compute_coding_offset(ivs, pos1, strand)
            if cod_offset is None:
                out_rows.append({**row.to_dict(), **{
                    'analysis_transcript':tid, 'analysis_result':'offset_fail'
                }})
                continue
            # determine codon start offset (0-based within CDS)
            codon_index = cod_offset // 3
            codon_start_offset = codon_index * 3
            # map codon_start_offset to three genomic 1-based coords
            coords = offsets_to_genomic_coords(ivs, strand, codon_start_offset, 3)
            if coords is None or len(coords) != 3:
                out_rows.append({**row.to_dict(), **{
                    'analysis_transcript':tid, 'analysis_result':'codon_coords_fail'
                }})
                continue
            # fetch bases for codon (genome)
            bases = []
            for p in coords:
                b = fetch_base_samtools(args.fasta, chrom, p)
                bases.append(b if b is not None else 'N')
            # assemble codon in mRNA orientation: if strand '+' codon = bases as is; if '-' codon = revcomp(bases joined)
            genomic_codon = ''.join(bases).upper()
            if strand == '+':
                mrna_codon = genomic_codon.replace('T','U')  # mRNA with U (we'll use T for translation)
                mrna_codon_for_translate = genomic_codon  # use T-based codon table
            else:
                # base seq is from genomic positive-strand; for minus strand, mRNA codon is reverse complement
                rc = revcomp(genomic_codon)
                mrna_codon_for_translate = rc  # T-based
            # compute mRNA base at site position within codon:
            # need to determine index of pos within the codon (0/1/2)
            # cod_offset - codon_start_offset gives position inside codon (0-based)
            pos_in_codon = cod_offset - codon_start_offset
            # but pos_in_codon is relative to transcript order; for '-' strand the mapping handled above
            # get the current base in mRNA (should be 'A' for m⁶A)
            current_base = mrna_codon_for_translate[pos_in_codon]
            is_mrna_A = (current_base.upper() == 'A')
            # simulate mutation mRNA A->G (if current base is A)
            mutated_codon = list(mrna_codon_for_translate)
            mutated_codon[pos_in_codon] = 'G'
            mutated_codon = ''.join(mutated_codon)
            orig_aa = GCODE.get(mrna_codon_for_translate.upper(), 'X')
            mut_aa = GCODE.get(mutated_codon.upper(), 'X')
            effect = 'synonymous' if orig_aa == mut_aa else 'nonsynonymous'
            # genomic-level change that corresponds to mRNA A->G:
            # if strand == '+': genomic A -> G
            # if strand == '-': genomic T -> C (because complement)
            if strand == '+':
                genomic_ref = fetch_base_samtools(args.fasta, chrom, pos1).upper()
                genomic_change = f"{genomic_ref}->G"
            else:
                genomic_ref = fetch_base_samtools(args.fasta, chrom, pos1).upper()
                # genomic_ref should be complement of mRNA base; mRNA A corresponds to genomic T
                # the genomic nucleotide change equivalent to mRNA A->G is T->C
                genomic_change = f"{genomic_ref}->C"
            # record
            rec = {
                **row.to_dict(),
                'analysis_transcript': tid,
                'transcript_strand': strand,
                'codon_genomic_coords_1based': ';'.join(map(str,coords)),
                'genomic_codon': genomic_codon,
                'mrna_codon': mrna_codon_for_translate,
                'pos_in_codon': pos_in_codon,
                'mRNA_base': current_base,
                'is_mRNA_A': is_mrna_A,
                'genomic_change_equiv': genomic_change,
                'orig_AA': orig_aa,
                'mut_AA': mut_aa,
                'effect': effect
            }
            out_rows.append(rec)
            if is_mrna_A:
                summary['mRNA_A_sites'] += 1
                if effect == 'nonsynonymous':
                    summary['nonsynonymous_hits'] += 1

    # write out
    outdf = pd.DataFrame(out_rows)
    outdf.to_csv(args.out, index=False)
    print("Wrote", args.out)
    print("SUMMARY:", summary)

if __name__ == '__main__':
    main()
