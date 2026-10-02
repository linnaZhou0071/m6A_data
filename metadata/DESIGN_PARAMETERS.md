# PEGG design parameters and selection status

This audit distinguishes what the **recorded custom script** specifies from the locally installed PEGG 2.1.0 `prime.run` signature. It does not establish the version or complete settings of either historical result file. The local `pegg_env` currently cannot import PEGG because its NumPy/cyvcf2 binary versions conflict, so a fresh end-to-end run has not been performed.

| Parameter or step | Recorded custom call | PEGG 2.1.0 default / behavior | Interpretation |
| --- | --- | --- | --- |
| Input format | `cBioPortal` explicitly | Required input choice | Genomic-plus variant rows |
| Genome | Loaded GRCh38 `chrom_dict` explicitly | Context lookup | Match genome build and reference alleles |
| `RTT_lengths` | `[20,22,25,27,30,32,35]` explicitly | `[5,10,15,25,30]` | **Not a default** |
| `min_RHA_size` | `6` explicitly | `1` | PEGG discards generated guides below the selected RHA minimum |
| `sensor` | `False` explicitly | `True` | Sensor sequence omitted |
| `PAM` | Not passed | `NGG` | PEGG searches eligible PAMs during generation; no separate post-design PAM-quality threshold was applied |
| `PBS_lengths` | Not passed | `[8,10,13,15]` | Inherits default if using this installed version |
| `rankby` | Not passed | `PEGG2_Score` | PEGG also calculates `RF_Score`; rank 1 is a model ranking, not validation |
| `pegRNAs_per_mut` | Not passed | `All` | No top-N restriction within PEGG |
| `RE_sites` | Not passed | `None` | No restriction-site exclusion list |
| `polyT_threshold` | Not passed | `4` | PEGG annotates `contains_polyT_terminator`; the recorded call does not remove flagged rows |
| Post-PEGG filters | None recorded | Not part of this call | No efficiency cutoff, off-target filter, oligo-feasibility gate, or library-balance filter |

The output headers include `PAM`, `RTT_length`, `PBS_length`, `RHA_size`, `PEGG2_Score`, `RF_Score`, and `contains_polyT_terminator`. Thus “PAM and efficiency not filtered” should be read narrowly: PAM eligibility is already a design constraint, scores are calculated, but no **additional** release filter or experimental validation has been applied. See [PEGG's official quickstart](https://pegg.readthedocs.io/en/latest/quickstart.html) and [the recorded local call](../scripts/run_pegg_design.py).

The current full intersection files have 136,868 and 264,104 guide rows, both covering 3,796 of 4,075 input sites.
The smaller file actually contains RTT lengths 5/10/15/25/30 and minimum RHA 1; the larger contains RTT lengths 20/22/25/27/30/32/35 and minimum RHA 6. Both contain PBS lengths 8/10/13/15. These distributions are consistent with the script/default contrast, but still do not establish the exact historical software version or all passed arguments. The larger file is named `custom`, but filenames and row counts alone do not prove its exact generation parameters. The smaller historical run's complete arguments are unknown. This candidate snapshot is an explicitly labelled subset of those historical outputs; PEGG will not be rerun for this release. A future mature-workflow version could restore a compatible environment and compare new designs against the originals.

## Corrected-filter subset and no-guide audit

The corrected MANE filter retains 3,814 of the original 4,075
intersection PEGG inputs. All 3,814 are exact site/allele subsets of
the historical input; this reanalysis **does not rerun PEGG**.
Subsetting the historical full outputs gives 129,128 default-like
and 248,876 custom guide rows, each covering 3,564 retained sites.
The 250 retained sites without guides are listed in
local data/reanalysis/historical_pegg_coverage_audit.csv (not part of the new public GitHub snapshot).
For the default-like RTT/RHA pattern, 232 have no eligible NGG PAM
and 18 have eligible PAM(s) but no RTT meeting minimum RHA; for the
custom pattern, the split is 131 and 119. The audit implements the
installed PEGG 2.1.0 eligibility bounds and reproduces observed
presence/absence for every one of the 4,075 original sites. It is
**not** proof of exact parameter provenance, sequence-level
reproducibility, or biological editing efficiency.

The user's design inspiration is [Pierce et al., *Nature*,
2025](https://www.nature.com/articles/s41586-025-09732-2).
That paper explores **6–10 nt of RTT homology** in one sup-tRNA
screen and **total RTT lengths 21–36 nt** in a later optimization;
those are distinct measures in a different target system. Our
recorded custom call used total RTT lengths 20–35 from the discrete
list above and only a **minimum** RHA of 6. It did **not** enforce
RHA ≤10; the corrected-filter subset of the historical custom output
contains **182,384/248,876** guide rows with RHA >10 (observed maximum 34).
We cite the paper as rationale, not as evidence that our settings
duplicate its protocol.
