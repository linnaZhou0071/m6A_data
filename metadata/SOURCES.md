# Source provenance and redistribution boundary

| Local historical label | Source and exact selection | Primary citation |
| --- | --- | --- |
| `file1_s23` (HEK) | `SAC-seq_data.xlsx`, sheets 2 “HEK293 polyA” and 3 “HEK293 ribo-” | [Hu et al., m6A-SAC-seq](https://doi.org/10.1038/s41587-022-01243-z) |
| `file1_sheet1` (HeLa) | Same workbook, sheet 1 “HeLa polyA” | [Hu et al.](https://doi.org/10.1038/s41587-022-01243-z) |
| `file2` (HEK) | `GLORI_HEK293T_mRNA.xlsx`, first worksheet | [Liu et al., GLORI](https://doi.org/10.1038/s41587-022-01487-9) |
| `file3` (HeLa) | GEO GSE211303 file `GSE211303_hela.polya.wt.ftom.ftop.rep1.deep.hits.txt` | [Xiao et al., eTAM-seq](https://doi.org/10.1038/s41587-022-01587-6); [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211303) |

The eTAM `ftom.ftop` TXT was identified by exact position/methylation
comparison with all 40,096 nonempty `m6a_file3` entries in
`data/processed/HeLa_union.csv`. Re-run
`python scripts/verify_etam_source.py` while the source TXT is available;
it exits nonzero on a mismatch. The `ftom.ivt` TXT is not the matching
source.

The SAC-seq sheet label is HEK293, whereas GLORI is HEK293T. GLORI
[Supplementary Fig. 6a](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-022-01487-9/MediaObjects/41587_2022_1487_MOESM1_ESM.pdf)
compares GLORI with m6A-SAC-seq in HEK293T. We therefore allow
cross-method comparison while reporting the source workbook label
unchanged. This is a comparability limitation, not proof of identical
sample provenance.

The exact NCBI GCF, UCSC hg38, and MANE release 1.5 files, local SHA-256
values, and the requested Morales et al. MANE citation are recorded in
[reference file provenance](REFERENCE_FILES.md).

The downloaded supplementary workbooks, GEO TXT files, article PDFs,
genome FASTA, and annotation GTF are **not** relicensed under this
repository's MIT/CC BY 4.0 notices. The derived site tables may reproduce
third-party values; confirm each source's redistribution conditions
before pushing the dataset publicly. Until then, attribution and
local exclusion are the conservative default.

## Row-level source audit (2026-10-02)

Run `python scripts/verify_sac_glori_sources.py` with the excluded
workbooks at the repository root and
`python scripts/verify_etam_source.py` with the cited GEO TXT.
`python scripts/rebuild_source_unions.py` then reconstructs both
union and ≥80% tables in memory and compares every position, value
and flag with the versioned CSVs; all four tables pass.
The first verifier found **zero missing, extra or unequal site values**
at a 1e-9 percentage-point tolerance (needed only for historical
binary-float CSV tails):

| Source selection | Source rows | Unique union positions | Transformation |
| --- | ---: | ---: | --- |
| SAC HEK293 polyA + ribo- | 12,234 + 18,952 | 28,891 | Average two sheet values at 2,295 shared positions |
| SAC HeLa polyA | 10,892 | 10,892 | Direct percentage |
| GLORI HEK293T Sheet1 | 170,240 | 170,240 | `100 × mean(m6A_level_rep1, m6A_level_rep2)` |
| eTAM HeLa selected TXT | 40,096 | 40,096 | Direct methylation value; see its dedicated verifier |

Every group `m6a_mean` in the current union tables equals the
mean of its nonempty method columns to 1e-9 percentage points:
HEK 186,249 rows and HeLa 45,944 rows, no mismatches.
Among ≥80% sites, HEK has **26,162 one-method** and **1,461
two-method** values; HeLa has **10,840 one-method** and **224
two-method** values. No uncertainty or replicate-variance filter
was imposed. These tallies are about *available measurements*, not
biological replication.

## Redistribution review: what may be published?

This is a conservative inventory, **not legal permission**. A public
download link is not itself a CC BY license for every underlying record.

| Material | Current or future treatment | Reason and remaining decision |
| --- | --- | --- |
| Our scripts and tests; new QC annotations, workflow SVG, plot and prose | Code under [MIT](../LICENSE); original documentation, figures and aggregate counts under [CC BY 4.0](../LICENSE-CC-BY-4.0.md) | Original contributions only; no third-party image copied |
| Paper PDFs and source supplementary XLSX | **Do not upload copies** | The [SAC article](https://www.nature.com/articles/s41587-022-01243-z) and [GLORI article](https://www.nature.com/articles/s41587-022-01487-9) identify exclusive Springer Nature/licensor article rights; their supplement pages do not show a blanket CC BY grant. Link to publisher files and cite instead. [Permissions route](https://www.nature.com/reprints/permission-requests.html). |
| eTAM source TXT and reference FASTA/GTF | **Do not upload copies**; link to [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211303), NCBI/UCSC and MANE release | [NCBI policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) imposes no NCBI restriction on molecular-data distribution but cannot transfer or certify a submitter's rights; reference files have their own provenance. |
| Union/high/intersection/filtered CSVs with source methylation values | **Conditional; hold public upload pending owner/institution review or permission** | These tables copy or numerically transform many source values. Attribution is necessary but does not alone settle republication of bulk supplement-derived data. If not cleared, publish column definitions, code, hashes and source-download instructions; deposit only a permitted subset or no source-derived values. |
| Site-fate, corrected effect and guide tables | **Conditional; review source-derived columns and genome-sequence content** | New analysis is original, but some rows copy source coordinates/measurements or reference sequence. A narrowly selected guide table may be easier to justify than wholesale supplement mirrors, but no blanket permission is inferred. |
| Historical Prism/JPEG and `potential_SNP.xlsx` | Keep local until content/history review; do not feature as current outputs | Prism HEK histogram is stale; workbook's historical `SNP` label is unsupported by genotype evidence. |

For any future row-level data deposit, the owner should confirm the
institutional rights interpretation and exact file list before publication.
MIT and CC BY 4.0 apply only to rights the owner controls; the citations
above must accompany any permitted derivative dataset.
