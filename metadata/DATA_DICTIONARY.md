# Data dictionary

All CSVs are UTF-8 and have a header. `pos` is a site ID of the form
`chr<chromosome>_<1-based genomic position>_<RNA strand>`; the
genomic coordinate is not a transcript coordinate. A site is one unique
`pos`, while an effect row is one site–transcript assessment and a
PEGG result row is one candidate guide.

| Field | Meaning |
| --- | --- |
| `m6a_file1_s23`, `m6a_file1_sheet1` | SAC-seq methylation percentage: HEK is the mean across available sheets 2/3 at each site; HeLa is sheet 1 |
| `m6a_file2` | GLORI percentage = 100 × mean of the two workbook replicate fractions |
| `m6a_file3` | eTAM-seq `methylation` percentage from the `ftom.ftop` TXT |
| `in_file*` | Whether a method contributed a value for that site |
| `m6a_mean`, `m6a_mean_HEK`, `m6a_mean_HeLa` | Mean of available method percentages within a group |
| `genomic_feature` | CDS, UTR, or other, from strand-aware overlap with the local GTF |
| `analysis_transcript`, `effect` | Transcript assessed and predicted synonymous/nonsynonymous effect |
| `Original_Pos_ID` | Source `pos` carried into PEGG |
| `Reference_Allele`, `Tumor_Seq_Allele2` | Genomic A→G on RNA + strand, or T→C on RNA − strand |
| `pegRNA_rank` | PEGG's rank within a site; rank 1 is the review extract |
| `PEGG2_Score`, `RF_Score` | PEGG model scores, not experimentally measured editing efficiency |

Historical `site_fate_*.csv` ledgers record old effect/filter status. Prefer
`site_fate_*_corrected.csv`: one row per original site, with
`coding_review_status` (`non_CDS_by_annotation`,
`MANE_synonymous`, `CDS_no_MANE_effect`, or
`MANE_nonsynonymous_excluded`), `has_mane_effect`,
`nonsynonymous_in_any_transcript`, corrected retention and PEGG-input
flags. `retained_after_historical_filter` and
`newly_excluded_after_strand_correction` expose the exact change from
the original pipeline (261 intersection, 1,204 union newly excluded).
`corrected_design_status` is `excluded_before_design`,
`unsupported_chrM`, `not_attempted` (union),
`historical_guide_reused`, or
`no_guide_under_historical_custom_parameters`. The intersection
ledger adds historical guide counts and parameter-specific coverage
reasons for **both** runs. A site with zero guides may still have been
in PEGG input.

`pegg_intersection_rank1_*_corrected_filter.csv` has one historical
PEGG rank-1 guide per **corrected retained site with a design** (3,564
sites per run). It does not contain all candidate guides or a validated
best guide. The default/custom labels denote historical runs; only the
custom parameters presently appear in the active PEGG script. The
locally staged `data/reanalysis/full_guides/` files contain the full
corrected-filter subsets; they are ignored by Git and are not published.

`data/reanalysis/historical_pegg_coverage_audit.csv` has one row per
original PEGG input and historical run. `coverage_reason` is
`candidate_eligible`, `no_eligible_NGG_PAM`, or
`no_RTT_meets_min_RHA`. These reasons mean *no design under the
inspected PEGG eligibility rules*, not that editing is biologically
impossible. `retained_in_corrected_filter` identifies the corrected
subset.

## Optional future Zenodo guide-only export (not published)

The two locally generated files under data/zenodo_preview/ contain 19
columns, one historical candidate guide per row, and no source
methylation percentages or long genome-context strings. They remain
local for a possible future rights-reviewed deposit; no public data
archive or dataset DOI is currently planned.
Original_Pos_ID encodes the selected GRCh38 genomic site and RNA
strand, so this is still a source-derived coordinate selection.
The inherited PEGG header Tumor_Seq_Allele2 names the edited DNA
allele; it does not mean the source sample is a tumour.
