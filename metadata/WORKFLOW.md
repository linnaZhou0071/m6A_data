# Analysis workflow, initial outputs and corrected reanalysis

The nine numbered diagram steps show the analysis sequence: curated
measurements → exploratory ≥80% sets → intersection/union → region
annotation → same-strand transcript mapping and RNA-A check → amino-acid
prediction → comparison with MANE Select v1.5 → unique-site exclusion →
PEGG design on the intersection input. The final diagram box records a
later correction to the negative-strand codon calculation. The corrected
coding effects and MANE exclusions were recomputed, but PEGG was not
rerun; previously designed guides were reused for retained sites.

The site ID already contains chromosome, 1-based coordinate and RNA strand.
The set operation first forms the intersection/union. Region annotation
then classifies CDS/UTR/other without predicting amino-acid effects. For
CDS sites, a separate step maps same-strand transcripts, compares complete
site IDs and checks the mapped RNA base. The historical non-A separation
was record-level: another transcript record for the same site could remain.

The reference roles are [listed with versions and hashes](REFERENCE_FILES.md):
steps 4–5 read `hg38.ncbiRefSeq.gtf.gz` for regions and transcript models;
step 6 reads `GRCh38.fa`, a decompressed copy of the NCBI
`GCF_000001405.26_GRCh38_genomic.fna.gz` sequence, for codon bases;
step 7 reads `MANE.GRCh38.v1.5.refseq_genomic.gtf.gz`; and step 9 uses
the compressed NCBI GCF FASTA as PEGG's genome input. The GCF FASTA is
**not** the file that assigns CDS/UTR labels in step 4.

## High-methylation sets to cross-group intersection/union

`data/processed/HEK_over_80%.csv` has 27,623 unique sites and
`data/processed/HeLa_over_80%.csv` has 11,064. Despite their older
filenames, the criterion documented for this release is **≥80%**.
`scripts/merge_csv_by_pos.py` joins these tables by the exact
`pos` string. A site present in both becomes one row of
`HEK_HeLa_over80_intersection.csv` (5,646 sites); a site present in
either becomes one row of `HEK_HeLa_over80_union.csv` (33,041 sites).
The latter equals 27,623 + 11,064 − 5,646.

The review-stage tables are local in `data/processed/` and are
excluded from the new public snapshot pending source-value rights review.
Historical copies are retained locally; the clean public snapshot does
not contain them. To
**rebuild without overwriting** either copy:

```bash
mkdir -p /tmp/m6a-merge-check
python scripts/merge_csv_by_pos.py \
  --input1 'data/processed/HEK_over_80%.csv' \
  --input2 'data/processed/HeLa_over_80%.csv' \
  --label1 HEK --label2 HeLa \
  --exclude-columns chrom posnum strand \
  --out /tmp/m6a-merge-check/HEK_HeLa_over80
cmp /tmp/m6a-merge-check/HEK_HeLa_over80_intersection.csv data/processed/HEK_HeLa_over80_intersection.csv
cmp /tmp/m6a-merge-check/HEK_HeLa_over80_union.csv data/processed/HEK_HeLa_over80_union.csv
```

The explicit column exclusion reproduces the 11 historical columns;
reading CSV fields as strings preserves `TRUE/FALSE` and numbers such
as `100` rather than rewriting them as `True/False` and `100.0`.
Both regenerated files passed byte-for-byte comparison locally.

## From sites to PEGG input

The branch-specific input to `filter_non_mane_nonsynonymous.py` is
the **full** cross-group site table, not only its CDS subset. First,
the CDS effect file assesses site–transcript pairs; then
`filter_mane_nonsynonymous.py` selects effect rows labelled
`nonsynonymous` whose transcript is in the MANE Select GTF.
[MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/) provides a
matched NCBI/Ensembl representative transcript at each locus. Using
release 1.5 defines a consistent primary coding reference for this
screen; it does not certify every alternative isoform as neutral.
`filter_non_mane_nonsynonymous.py` collects unique `pos` values from
those rows and removes those sites from the full intersection/union
table. The precise observed transitions are:

| Branch | Full sites | MANE nonsynonymous effect rows | Distinct transcript IDs | Unique sites removed | Retained sites |
| --- | ---: | ---: | ---: | ---: | ---: |
| Intersection | 5,646 | 1,571 | 1,005 | 1,571 | 4,075 |
| Union | 33,041 | 8,507 | 3,968 | 8,493 | 24,548 |

The unit changes at the mapping step: the original table has **one row
per site**, whereas the CDS effect table has **one row per site–transcript
mapping**. A site may appear in several effect rows, and a transcript may
appear for several sites. Therefore 8,507 is neither 8,507 unique
transcripts nor 8,507 unique sites. The filter implements:

```python
excluded_sites = set(mane_nonsyn_rows['pos'].unique())
filtered_sites = original_sites[~original_sites['pos'].isin(excluded_sites)]
```

For the historical union, the calculation is 33,041 − 8,493 = 24,548;
for the intersection, 5,646 − 1,571 = 4,075. One union site,
`chr5_141009819_+`, accounts for 15 rows; the 14 extra rows are not 14
extra sites. The historical MANE-filtered CSVs were audited locally for
these row, transcript-ID and site counts. The filter keeps all UTR and
`other` sites: the retained intersection is 2,570 UTR, 318 other,
and 1,187 CDS. It also retains sites that lack an effect record.

Re-running the current CSV filter scripts in a temporary directory
reproduced both retained site sets and every value after normalizing
numeric display and boolean case. The historical filtered CSVs are not
byte-identical to these fresh outputs: some fields display as `False/True`
or `100.0` instead of `FALSE/TRUE` or `100`. That is a serialization
difference, not a different site selection.

`scripts/prepare_pegg_input.py` parses each `pos`, writes a
1-based single-nucleotide SNP in PEGG's cBioPortal-style CSV, and
carries the original ID as `Original_Pos_ID`. It does not perform
a new coding-effect test. On RNA `+` sites, the genomic plus-strand
alleles are `A→G` (2,047 intersection rows); on RNA `−` sites,
the equivalent plus-strand genomic alleles are `T→C` (2,028 rows).
The laboratory note about converting everything to the positive strand means representing edits in
the reference genome's plus-strand coordinates. It does **not** turn
every site into genomic `A→G`, and the original RNA strand remains
recoverable from `Original_Pos_ID`.

All 4,075 filtered intersection sites enter PEGG input. In the union,
the converter supports `chr1–chr22/X/Y` but not `chrM`, so it
skips `chrM_16559_+` and `chrM_5795_-`; 24,546 remain. The
remaining union PEGG-input alleles are 13,198 `A→G` and 11,348
`T→C`.

“Does not change codon” is not established by this historical filter.
A coding A→G substitution changes the nucleotide codon; it may or may
not change the encoded amino acid. The actual policy is “exclude
observed MANE Select nonsynonymous calls.” Of retained sites, 28
intersection and 181 union sites still have a nonsynonymous call on
another transcript. Moreover, the existing effect tables contain no
record for 632 intersection and 2,496 union CDS sites, all on the RNA
negative strand. Therefore `Safe/Neutral` in the old notes is a
working label, **not** a validated property of every retained site.
## Why 8,507 effect rows remove 8,493 sites

The union MANE-nonsynonymous CSV has one row per **site–transcript
effect**, while `filter_non_mane_nonsynonymous.py` excludes distinct
`pos` values from the full union site table. Exactly 8,492 positions
have one qualifying row. One position, `chr5_141009819_+`, has **15**
qualifying rows from distinct `NM_...` transcript IDs. Therefore
`8,492 × 1 + 15 = 8,507` effect rows but `8,492 + 1 = 8,493`
unique positions. That site's recorded RNA codon is `GAC` and the
A→G effect is `GAC → GGC` (Asp → Gly) in these records.

Verify directly with the versioned table:

```bash
python - <<'PY'
import csv
from collections import Counter
path = "data/intermediate/HEK_HeLa_over80_union_CDS.AtoG.effects.MANE_nonsynonymous.csv"
with open(path, newline="") as handle:
    rows = list(csv.DictReader(handle))
counts = Counter(row["pos"] for row in rows)
print("effect rows:", len(rows))
print("distinct sites:", len(counts))
print("positions with >1 row:", [(pos, n) for pos, n in counts.items() if n > 1])
print("transcripts for the repeated site:")
print([row["clean_transcript_id"] for row in rows if row["pos"] == "chr5_141009819_+"])
PY
```

Expected counts: `8507`, `8493`, and one repeated position with
`15` rows. The transcript IDs can be inspected in the source CSV;
their multiplicity is not 15 independent m⁶A measurements.

## Corrected reanalysis (2026-10-02)

The historical coding-effect script reversed an already transcript-ordered
negative-strand codon a second time. For a negative-strand site, genomic
coordinates returned by `offsets_to_genomic_coords` are already in
transcript 5′→3′ order, so their genome bases must be **complemented**
without another reversal. It also used to include opposite-strand CDS
transcripts. The current [coding-effect script](../scripts/predict_AtoG_coding_effects.py)
corrects both issues and uses indexed FASTA access. Five focused strand/
exon-boundary tests pass.

Yuling Zhou recalls using the hg38 transcript mapping to set aside
rows whose inferred transcript RNA base was not A and saving them as
`potential_SNP.xlsx`. The locally retained workbook has 635 intersection
and 2,522 union unique positions; 635 and 2,521 respectively are on the
negative RNA strand. It includes 632/2,496 CDS sites absent from the old
effect CSVs. The other 3/26 workbook sites also occur in the old effect
CSVs, demonstrating that the old separation applied to site–transcript
records rather than automatically removing each whole site. This was an
intentional QC separation based on the old non-A flags, but those flags
were largely caused by the old
negative-strand codon-orientation error, not by observed genotypes.
The workbook is **not** evidence that any site is a biological SNP.
The historical workbook and effect files have not been overwritten;
genotyping the target cell line remains future work. See the
[site–transcript audit](POTENTIAL_SNP_AUDIT.md).
Specifically, `chr10_11462944_-` appears in that workbook with genomic
bases `TGA` at transcript-ordered coordinates
`11462944;11462943;11462942`: the historical RNA codon is `TCA`
and `is_mRNA_A=False`, whereas the corrected codon is `ACT` with
RNA base `A`. One positive-strand workbook record,
`chr6_20402333_+`, came from an opposite-strand transcript and is
excluded by the new same-strand match. These are computational
orientation/annotation errors, not cell-genotype evidence.

| Branch | Corrected CDS sites / effect rows | MANE-nonsynonymous sites excluded | Retained sites | PEGG input | Retained sites with a historical guide |
| --- | ---: | ---: | ---: | ---: | ---: |
| Intersection | 2,758 / 21,160 | 1,832 | 3,814 | 3,814 | 3,564 |
| Union | 13,279 / 100,065 | 9,697 | 23,344 | 23,342 | Not designed |

Relative to the historical retained tables, the correction newly
excludes **261 intersection** and **1,204 union** sites; it restores
no previously excluded site. The exact IDs are flagged by
`newly_excluded_after_strand_correction` in the corrected
site-fate ledgers.

All corrected CDS sites have effect rows; every row has
`analysis_result=ok` and `is_mRNA_A=True`. Among retained CDS sites,
40 intersection and 262 union sites have no matching MANE effect row;
33 and 215 respectively have a nonsynonymous effect on a non-MANE
transcript. These are still retained under the **chosen MANE-only**
policy, and the per-site status is explicit in
`data/summary/site_fate_*_corrected.csv`.

From the repository root, with the excluded reference files present:

```bash
python -m unittest discover -s tests -v
for branch in intersection union; do
  stem="HEK_HeLa_over80_${branch}"
  python scripts/predict_AtoG_coding_effects.py \
    --gtf hg38.ncbiRefSeq.gtf.gz --fasta GRCh38.fa --assume-cds \
    --out "data/reanalysis/${stem}_CDS.AtoG.effects.corrected.csv" \
    "data/intermediate/${stem}_CDS.v3fixed.csv"
  python scripts/filter_non_mane_nonsynonymous.py \
    --mane-gtf MANE.GRCh38.v1.5.refseq_genomic.gtf.gz \
    --effects-csv "data/reanalysis/${stem}_CDS.AtoG.effects.corrected.csv" \
    --original-csv "data/processed/${stem}.csv" \
    --output-csv "data/reanalysis/${stem}_filtered.corrected.csv"
  python scripts/prepare_pegg_input.py \
    --input-file "data/reanalysis/${stem}_filtered.corrected.csv" \
    --output-file "data/reanalysis/pegg_input_${branch}.corrected.csv"
done
python scripts/verify_reference_alleles.py --input-dir data/reanalysis --suffix .corrected --all-union-sites
python scripts/audit_pegg_coverage.py
python scripts/build_corrected_release.py --full-guides
```

The final command subsets the **historical** full PEGG outputs to the
corrected retained sites; it does not execute PEGG. The original outputs
had 136,868 default-like and 264,104 custom guide rows across 3,796
sites. Their corrected-filter subsets have 129,128 and 248,876 rows
across 3,564 sites; full subsets are locally staged under the
Git-ignored `data/reanalysis/full_guides/` for possible future review. Rank-1
extracts and site-fate ledgers are in `data/summary/`.

The local coverage audit (data/reanalysis/historical_pegg_coverage_audit.csv)
implements the installed PEGG 2.1.0 NGG window and RTT/RHA eligibility
rules and asserts that predicted guide presence matches both historical
outputs for **all 4,075** original intersection inputs. For the 279
historical no-guide sites, the default-like parameter set yields 258
with no eligible NGG PAM and 21 with no RTT satisfying minimum RHA;
the custom set yields 150 and 129 respectively. For the **250** retained
corrected inputs without guides, those splits are 232/18 and 131/119.
This is a parameter-specific computational explanation, not a failed
experimental screen.

## Distribution-plot provenance

The locally updated Prism project `m6A_plotting.pzfx` contains 186,249
raw HEK values and 45,944 raw HeLa values. Its saved HEK ten-bin table
now contains 1,437, 37,302, 31,521, 24,126, 19,449, 17,196, 14,307,
13,288, 12,963 and 14,660 sites: total **186,249**, matching the
[union-derived count table](../figures/m6a_distribution_counts.csv) in
every bin. The revised researcher-authored
[HEK](../figures/m6a_distribution_hek_prism.jpg) and
[HeLa](../figures/m6a_distribution_hela_prism.jpg) JPEGs are featured in
the README. The HEK JPEG title is “HEK”, correctly describing the
mixed HEK293/HEK293T analysis group; an internal Prism table-column
title still reads `HEK293T` and should not be interpreted as the sample
label for the combined distribution.

Run `python scripts/plot_m6a_distribution.py` to rebuild the diagnostic
[SVG](../figures/m6a_distribution_corrected.svg), PNG and
[bin-count CSV](../figures/m6a_distribution_counts.csv) from the
versioned union CSVs. Bins are [0,10), ..., [80,90), [90,100], with
100% in the final bin. The current union-based HEK counts sum to 186,249
and HeLa to 45,944. Histogramming the Prism HEK raw values directly previously placed one
boundary value in a different 40–50%/50–60% bin: a union mean stored as
49.99999999999999 displays as 50 in Prism. The corrected saved ten-bin
table follows the union-CSV bin convention and matches every published
count. This is a floating-point bin-boundary issue, not a site-count
difference.

The source-to-union step is now reconstructed by
`python scripts/rebuild_source_unions.py`: using the exact cited
workbooks/TXT, it rebuilds both group unions and ≥80% tables in memory
and compares **every site, value and flag** against the versioned CSVs.
All four tables passed (186,249, 27,623, 45,944 and 11,064 rows).
Use `--output-dir /some/new/directory` only if normalized output copies
are desired; it refuses to overwrite existing files. Binary-float
serialization and UTF-8 BOMs may differ, but values match within
1e-9 percentage points. The exact historical PEGG environment is
still not reconstructed.
