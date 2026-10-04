# Why the historical potential_SNP workbook is not SNP evidence

Yuling Zhou recalls mapping the sites to hg38 RefSeq transcripts, separating
records whose inferred transcript RNA base was not A, excluding those records
from the main coding-effect table, and saving them as
`local/legacy/potential_SNP.xlsx`. This explains the historical screening
**decision**. The exact export command was not retained, and the workbook is
not a cell-line genotype call set.

| Sheet | Rows | Distinct site IDs | Negative-strand site IDs | Workbook rows labelled `is_mRNA_A=False` | Site IDs absent from the old effect CSV | Site IDs present in corrected effects with RNA A |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Intersection | 4,760 | 635 | 635 | 4,760 | 632 | 635 |
| Union | 19,372 | 2,522 | 2,521 | 19,372 | 2,496 | 2,522 |

The workbook has multiple transcript rows per site; row counts are not
counts of independent m⁶A sites. Three intersection sites and 26 union
sites also occur in the old effect CSVs, so the historical non-A
separation did **not** uniformly remove entire sites. It separated
particular site–transcript records; another record could remain for
the same site. In the corrected analysis, 4,758/4,760
intersection and 19,369/19,372 union **site–transcript pairs** still match;
all those matched corrected rows have RNA base A. The two unmatched pairs
at `chr2_112738913_-` use plus-strand transcripts for a minus-strand site.
The additional unmatched union pair at `chr6_20402333_+` uses a
minus-strand transcript for a plus-strand site. These three old pairs are
excluded by the corrected same-strand rule.

A concrete negative-strand example is `chr10_11462944_-`. In transcript
5′→3′ coordinate order, its genomic bases are `TGA`. The old script
[reverse-complemented](../scripts/legacy/predict_AtoG_coding_effects_v2_historical.py)
those already ordered bases and reported RNA codon `TCA`, RNA base
`T`, and `is_mRNA_A=False`. The corrected
[script](../scripts/predict_AtoG_coding_effects.py) only complements
them, yielding RNA codon `ACT`, RNA base `A`, and
`is_mRNA_A=True` for the same site and transcript. Reversing the
genomic bases a second time moved the target to the wrong codon
position. The test
[test_reverse_codon_is_complemented_not_reversed_twice](../tests/test_coding_effects.py)
also encodes this example.

Run `python scripts/audit_potential_snp.py` from the repository root to
recheck the aggregate counts above. It requires the workbook, old effect
CSVs and corrected effect CSVs retained locally; these are not in the
public GitHub snapshot. For the example, inspect matching `pos` and
`analysis_transcript` rows in the workbook and corrected effect CSV.

The workbook accurately records which rows the **old calculation** called
non-A; the researcher then set those rows aside. However, the matched
site–transcript comparisons and the codon example show that the old
negative-strand calculation produced erroneous non-A calls. Calling the
workbook `potential_SNP` was a provisional interpretation, not evidence
of an observed variant. Real variants in the intended cell line remain
possible and would require cell-line genotyping or an independent variant
call set.
