#!/usr/bin/env python3
import gzip
from collections import defaultdict

gtf = "hg38.ncbiRefSeq.gtf.gz"

exons = defaultdict(list)   # (chr, tid, strand) -> [(s,e)]
cds   = defaultdict(list)

def parse_attr(s):
    d={}
    for x in s.strip().split(';'):
        if '"' in x:
            k,v=x.strip().split(' "',1)
            d[k]=v.rstrip('"')
    return d

with gzip.open(gtf,'rt') as f:
    for l in f:
        if l.startswith('#'): continue
        c=l.rstrip().split('\t')
        if len(c)<9: continue
        chr_,_,typ,s,e,_,strand,_,attr=c
        a=parse_attr(attr)
        tid=a.get('transcript_id')
        if not tid: continue
        s,e=int(s)-1,int(e)
        key=(chr_,tid,strand)
        if typ=="exon":
            exons[key].append((s,e))
        elif typ=="CDS":
            cds[key].append((s,e))

def merge(xs):
    xs=sorted(xs)
    out=[]
    for s,e in xs:
        if not out or s>out[-1][1]:
            out.append([s,e])
        else:
            out[-1][1]=max(out[-1][1],e)
    return out

def subtract(a,b):
    out=[]
    for s,e in a:
        cur=s
        for cs,ce in b:
            if ce<=cur: continue
            if cs>=e: break
            if cs>cur: out.append((cur,cs))
            cur=max(cur,ce)
        if cur<e: out.append((cur,e))
    return out

with open("CDS.bed","w") as fc, open("UTR.bed","w") as fu:
    for k,ex in exons.items():
        chr_,tid,strand=k
        exm=merge(ex)
        cdsm=merge(cds.get(k,[]))
        for s,e in cdsm:
            fc.write(f"{chr_}\t{s}\t{e}\t{tid}\t.\t{strand}\n")
        if not cdsm: continue
        utr=subtract(exm,cdsm)
        for s,e in utr:
            fu.write(f"{chr_}\t{s}\t{e}\t{tid}\t.\t{strand}\n")
