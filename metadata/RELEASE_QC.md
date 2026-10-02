# QC and publication boundaries

**Current application scope (2026-10-02):** GitHub v0.3.0 is a public code/documentation snapshot under a new single-root public history; the corrected row-level dataset remains local and unpublished. No Zenodo deposit is planned for this application. The release gates below are an audit inventory for a possible future public dataset or mature workflow, not conditions for calling computational candidates validated reagents.

Status key: **PASS** = independently checked on the current local snapshot; **PARTIAL** = some evidence exists but the acceptance criterion is not met; **OPEN** = not yet verified; **FAIL** = a known problem blocks the named claim. Ownership: **Agent** = code/data audit and documentation I can do locally; **Owner** = Yuling Zhou or other rights/scientific decision-maker; **Both** = agent prepares evidence, owner signs off.

## What has actually been checked

- The **historical** stored [QC report](qc_report.json) says the ≥80% high tables have 27,623 HEK and 11,064 HeLa distinct sites; exact-`pos` intersection/union are 5,646/33,041; MANE-subtraction, PEGG-input subsets, rank-1 cardinality, and guide-ID uniqueness pass its listed Boolean checks. These are **set/count checks**, not proof that every source call or coding-effect assignment is biologically correct. The report currently does not cause `build_release.py` to exit nonzero if a Boolean is false.
- [Corrected QC report](corrected_qc_report.json) and local site-fate ledgers now cover all 2,758/13,279 CDS sites and 21,160/100,065 effect rows with `analysis_result=ok` and `is_mRNA_A=True`; corrected MANE-only filtering retains 3,814/23,344 intersection/union sites, newly excluding 261/1,204 sites relative to history. Thirteen focused unit tests pass, including exact-80 threshold inclusion and Zenodo export checks. The corrected input REF check passes all 3,814/23,342 nuclear variants, and an RNA-strand-aware check passes all 33,041 high-methylation union sites, including the two mitochondrial positions.
- [Source-to-union rebuilder](../scripts/rebuild_source_unions.py) recreated both union and ≥80% tables from the cited SAC/GLORI/eTAM inputs and matched all positions, values and flags (186,249/27,623 HEK and 45,944/11,064 HeLa rows) without overwriting them.
- [SAC/GLORI source verifier](../scripts/verify_sac_glori_sources.py) matched all 28,891 SAC HEK, 10,892 SAC HeLa and 170,240 GLORI positions/values to their union columns (1e-9 percentage-point tolerance for float serialization); all 186,249 HEK and 45,944 HeLa group means also recompute. Of the ≥80% sites, 26,162 HEK and 10,840 HeLa values have only one supporting method.
- The local PEGG coverage audit explains every one of the 279 historical no-guide sites under each parameter set and exactly matches observed presence/absence. Corrected-filter subsets of historical output contain 129,128/248,876 guide rows covering 3,564 retained intersection sites. The original full outputs remain untouched; the lean Zenodo preview has the same row/site counts without methylation values.
- [Corrected distribution counts](../figures/m6a_distribution_counts.csv) sum to the current union totals; the historical Prism HEK manual histogram sums to 186,803 and must not be used as the current 186,249-site plot.
- [eTAM source verification](../scripts/verify_etam_source.py) was rerun read-only on 2026-10-02: 40,096 source rows = 40,096 unique HeLa `m6a_file3` records; zero missing, extra, or unequal methylation values. This verifies the selected TXT-to-union mapping, **not** the entire raw-to-union pipeline.
- [Reference-allele verification](../scripts/verify_reference_alleles.py) was rerun read-only on 2026-10-01: zero GRCh38 REF mismatches across 4,075 intersection and 24,546 union PEGG inputs. This does **not** prove the intended cell line carries that allele or that editing is feasible.
- The two union `chrM` skips persist in the corrected input. Historically, 8,507 union MANE-nonsynonymous effect rows reduced to 8,493 distinct removed sites; these are **not** the corrected counts (9,697 excluded sites). The observed historical RTT/PBS/minimum-RHA distributions are recorded in [design parameters](DESIGN_PARAMETERS.md).
- The workflow and corrected distribution SVGs parse as XML; README links need a final local and public-clone audit after the last edit. Original source downloads and full PEGG results are excluded from the proposed new Git content. This does **not** audit the pre-existing Git history or establish distribution rights.
- The older 27-entry SHA-256 manifest passed on 2026-10-01 but predates this reanalysis. The new local inventory script separates code/docs, rights-held derived tables, large Zenodo candidates and historical-only files; it must be regenerated at the release freeze.
- A basic search of committed text found no obvious common credential patterns; binary files, complete history, GitHub settings, and collaborator access have **not** been cleared by a full security/rights review.

## Level A — gates before any public candidate-dataset release

### A1. Define the claim and uncertainty labels — PARTIAL; Both

**Why:** “m6A site”, “protein-preserving edit”, “PEGG candidate”, and “validated screening guide” are not interchangeable. Overclaiming is a data-quality defect even if every CSV parses.

**Do / pass:** Owner approved the scope statement: published cross-method measurements in specified samples; exploratory ≥80% mean rule; computational A→G/T→C candidates for future PE screening; no functional or editing validation. Final public wording still needs owner review. In every table/figure/README, label unknown CDS effects, transcript disagreements, absent designs, and whether counts mean sites, transcript-effect rows, or guide rows. Do not call the corrected 3,814 or historical 4,075 “safe/neutral” or the rank-1 extracts an experiment-ready library. Agent can revise labels and add machine-readable uncertainty columns after policy approval.

### A2. Source rights and split MIT/CC BY 4.0 boundary — OPEN for future data deposits; Owner with Agent inventory

**Why:** Derived CSVs reproduce source methylation values; MIT for original code and CC BY 4.0 for original documentation/figures cannot grant rights in SAC-seq/GLORI/eTAM-seq supplements or annotations. The corrected row-level dataset is not part of the public snapshot.

**Do / pass:** For each workbook, GEO TXT, image, derived column and reference, record URL, accession/version, citation, original terms, whether redistribution/derivative-value publication is allowed, and decision. Keep PDFs, source workbooks/TXT and reference FASTA/GTF outside the public Git tree unless independently permitted. If a source restriction is unclear, publish metadata/code and a download/rebuild recipe rather than assuming permission for source-derived values. Owner (or institution) gives final rights approval. See [SOURCES.md](SOURCES.md), [Creative Commons FAQ](https://creativecommons.org/faq/) and the scoped [LICENSE](../LICENSE).

### A3. Public history replacement — owner-directed; verify remote after push; Both

**Why:** The original public Git history contained two preliminary cross-group CSVs and `potential_SNP.xlsx`. The workbook is owner-created, but its SNP label is unsupported by genotype evidence. A clean single-root public commit removes these files from the reachable `main` history; it cannot remove clones, forks, cached views or copies fetched earlier. The previous commit hashes also change.

**Do / pass:** The owner requested the history rewrite after the v0.2.0 review. Before any force update, confirm the exact new tree, that no current branch or tag points to old commits, and that the remote has not advanced unexpectedly. Afterward verify one root commit and zero tracked row-level data files, and avoid merging an old clone back into `main`. GitHub explains the [limits and side effects of history removal](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/removing-sensitive-data-from-a-repository).


### A4. Source identity, sample comparability and raw lineage — PARTIAL; Both

**Why:** The HEK group combines a SAC worksheet labelled HEK293 with GLORI labelled HEK293T; biological replicates, protocol and coverage may differ. A merged coordinate alone does not make measurements directly comparable.

**Do / pass:** Agent has made the source-to-column matrix and checked SAC, GLORI and eTAM values row-for-row. The source notes document sheet/replicate handling, coordinate origin and strand. The owner accepts the wording that HEK293 and HEK293T are grouped for analysis but not asserted to be identical samples. A reviewer must be able to trace each released measurement to a source record or an explicit exclusion reason; the fully executable importer is the stronger B1 gate.

### A5. Percentage, missingness and ≥80% rule — PARTIAL; Both

**Why:** Means over available methods can privilege sites measured once; 80 versus >80 changes membership; missingness/coverage can dominate cross-method comparisons.

**Do / pass:** Source values are range-checked by the source verifiers, and every group mean recomputes from the available method columns. Presence flags, unique IDs and coordinate fields now pass automated checks for every HEK/HeLa union row. Report per-site number of contributing methods, missingness, support/coverage and replicate evidence where available; stratify threshold counts by method support. The high sets pass the stored ≥80% check (38 HEK sites equal exactly 80). One-method support dominates the high sets: 26,162/27,623 HEK and 10,840/11,064 HeLa. Owner decides whether a minimum-support or replicate-confidence rule is needed; otherwise disclose it as a limitation.

### A6. Coordinate, strand and reference validation — PARTIAL; Agent; Owner for target genotype

**Why:** An off-by-one position, wrong genome build, or RNA/genomic strand inversion creates a guide for the wrong base even if set arithmetic is correct.

**Do / pass:** Keep an explicit 1-based `chr_pos_RNAstrand` convention and a contig alias map; test representative positive/negative-strand and chromosome-edge sites against GRCh38. The historical REF check passed 28,621 inputs; the corrected check passed 27,156 inputs with zero mismatches. Every corrected CDS effect row is RNA-A and the input converter emits genomic A→G on RNA-positive and T→C on RNA-negative sites. All 33,041 high-methylation union sites also have the expected GRCh38 genomic A/T reference base. Motif/flanking-sequence and target-cell-genotype checks remain open. Before wet-lab use, confirm the intended cell-line genotype and structural variants; reference agreement alone is insufficient.

### A7. File manifest, schema and integrity — PARTIAL; Agent

**Why:** Counts may stay unchanged while column definitions, missing-value encodings, units or a file's bytes change. A checksum command that ignores missing files can produce a misleading partial pass.

**Do / pass:** Declare required public files and externally hosted files separately, with versions, row counts, unique-key rules, coordinate system, columns/units, NA values and SHA-256 hashes. Assert no duplicate/empty keys, no out-of-range values, and no silent skipped rows. Pick one canonical CSV serialization for `TRUE/FALSE`, numeric precision and line endings; regenerate and compare. The [data dictionary](DATA_DICTIONARY.md) and local-only file inventory are current review aids. The older local checksum manifest is historical; regenerate the local inventory after the final edit and create an exact public/Zenodo manifest only after rights and file-list approval. A public-clone check must fail if a **required** file is absent.

### A8. Accessibility and attribution of released outputs — PASS for labelled code snapshot; OPEN for future data release

**Why:** The public code snapshot and aggregate QC counts are inspectable, but it does not contain the corrected site-level or full PEGG CSVs. A clone alone cannot reconstruct historical guides without excluded source/reference files and full historical outputs.

**Do / pass:** Keep source citations, sample labels, field definitions and limitations public; label the corrected dataset unpublished and computational candidates unvalidated. If a future row-level deposit is chosen, review actual files and rights, provide a dictionary and hashes, preview the draft, and only then add its version-specific DOI. No DOI is required for the current application.


## Level B — required to call this a mature, reproducible workflow

### B1. Raw-to-high-table reconstruction — PASS on selected local sources; clean-clone test remains B6; Agent, then Owner review

**Why:** The original workbook/TXT-to-union, mean and threshold stage was previously unscripted. It is now reconstructed from the exact local source files; a new researcher still needs to obtain those excluded downloads.

**Done / remain:** The non-overwriting importer checks duplicate source rows, reconstructs group means and ≥80% membership, and compares every key/value/flag to all four versioned tables; an exactly-80 fixture passes. The selected local source versions are documented. A clean-clone/download test and broader synthetic source fixtures remain B6 work; old files remain intact.

### B2. Annotation and CDS-effect completeness — PASS for CDS coverage; PARTIAL for protein-neutral interpretation; Both

**Why:** The historical effect tables missed 632 intersection and 2,496 union negative-strand CDS sites because codons were reversed twice. Corrected effects now cover every annotated CDS site. Coverage alone does not make a retained edit protein-neutral.

**Done / remain:** The fixed script has plus/minus and exon-boundary tests; separately named effects contain no failed rows and all target RNA bases are A. The owner chose MANE-only exclusion. The site-fate ledgers explicitly mark retained CDS sites with no MANE effect (40 intersection, 262 union). Do not call those sites, or the entire retained set, universally protein-preserving.

### B3. Transcript and MANE exclusion policy — PARTIAL (MANE-only chosen and implemented); Both

**Why:** The historical filter removes only recorded MANE Select nonsynonymous calls. After correction it leaves 33 retained intersection and 215 retained union sites with a nonsynonymous call on another transcript. A coding triplet itself changes even when the amino acid does not.

**Done / remain:** The owner chose MANE-only exclusion; corrected filters and site-fate categories implement it, including no-MANE-effect and non-MANE disagreement flags. Confirm whether that is the final public scientific policy and do not describe retained edits as “not changing codon.”

### B4. PEGG-input attrition and no-guide explanations — PASS for historical coverage; Agent

**Why:** The corrected filter retains 23,344 union sites but the converter accepts 23,342; 250/3,814 corrected intersection inputs have no historical guide. The historical counts were 24,548→24,546 and 279/4,075. Without reasons, a reader cannot tell unsupported formats from PAM/RTT/RHA design failure.

**Do / pass:** The two omitted `chrM` IDs are recorded in corrected QC. Machine-readable status/reason exists for each corrected site and for each historical PEGG input under both parameter sets; predicted eligibility matches observed guide coverage. Union sites were **not** run through PEGG and are labelled not attempted. Distinguish “no candidate under these parameters” from “failed conversion” and from “not attempted.”

### B5. Historical design provenance — OPEN for exact reproduction; not a gate for the labelled candidate snapshot

**Why:** The custom script explicitly sets RTT lengths, min RHA and sensor, while the smaller run's command/version is not recorded. Output-column distributions support but do not prove those settings. The installed `pegg_env` currently fails to import PEGG due to NumPy/cyvcf2 binary incompatibility.

**Do / pass:** This release explicitly reuses historical guides and does **not** rerun PEGG, per owner decision. Record the observed parameter distributions and known custom call; label the other run `default-like`, not proven default. Exact software/version/argument provenance and sequence-level reproduction remain unverified and belong to a future mature-workflow version. Do not overwrite the historical full outputs. See [parameter audit](DESIGN_PARAMETERS.md).

### B6. Automated tests, hard failures and clean-clone run — OPEN; Agent

**Why:** The historical `build_release.py` writes Boolean checks but does not currently fail when one is false. The corrected builder does fail on its checked invariants, and 13 versioned small-fixture tests pass in the local Python 3.9 environment and a clean export of the tracked v0.2.0 tree. CI, a fresh environment install and an end-to-end run with separately downloaded data remain absent. Local full-data success may depend on files intentionally omitted from Git.

**Do / pass:** Plus/minus strand, exon boundaries, duplicate plotting keys, exactly-80 inclusion, transcript ID cleaning and `chrM` skipping now have tests. Still add missing-value, MANE multiple-row subtraction and rank-1 extraction fixtures. Make failed invariants produce a nonzero exit. From a fresh clone/environment, run documented steps with explicit downloaded inputs or a small fixture dataset; verify the file manifest and links. Mark full-data steps as requiring external downloads, not silently “reproducible” from the clone.

### B7. Figures, schema stability and release freeze — PARTIAL (plot fixed; freeze pending); Agent, Owner approval

**Why:** A figure without its input table/script cannot be regenerated; reusing a filename for changed content breaks citations. Different filtered CSV serializations already exist.

**Do / pass:** The corrected distribution figure has a regeneration script and count CSV with explicit bins/denominators; historical Prism/JPEG outputs remain marked legacy. Document or remove any other original plots before release. Review the data dictionary and decide which large intermediate tables are essential to review versus archival. Freeze a versioned schema, regenerate checksums after **all** content changes, run a final diff and independent spot checks, then tag a release and archive its exact bytes. No tag/DOI has yet been made.

## Level C — additional gates before calling any guide set screening-ready

### C1. Guide/editor compatibility and sequence QC — OPEN; Both

**Why:** PEGG's NGG PAM search and RHA minimum are generation constraints, not a complete quality filter. `PEGG2_Score`/`RF_Score` are predictions, not measured editing efficiencies.

**Do / pass:** Owner specifies editor/pegRNA scaffold, PAM, delivery cell line, desired edit, acceptable bystanders and control strategy. Agent can screen candidate PAM status/disruption, PBS/RTT/RHA and nick geometry, polyT/restriction sites, sequence repeats, synthesis constraints, predicted efficiency and off-target candidates against the intended genome. Choose thresholds **before** selecting winners and retain failure reasons; validate genotype and editing in the target cells. No arbitrary threshold is applied to the current rank-1 extracts.

### C2. Library composition, controls and biological validation — OPEN; Owner-led

**Why:** One ranked guide per site does not provide adequate redundancy, control distribution, or evidence that a phenotype is mediated by m6A rather than altered RNA sequence/protein/splicing.

**Do / pass:** Plan multiple independent designs per target when feasible, non-targeting and positive/negative controls, library representation and sequencing QC, and an assay-specific power/replicate plan. Experimentally measure edit rates, m6A changes, protein/splicing and phenotype with orthogonal confirmation. Motif-disruption alternatives for nonsynonymous A→G sites need their own coding/transcript/motif analysis. These are **not** prerequisites for honestly releasing a *computational candidate dataset*, but they are prerequisites for claiming a validated screening library.

## Current decision and ownership

1. **Current application output:** public MIT-licensed code and CC BY 4.0-licensed original documentation/figures; the corrected row-level dataset stays local and unpublished. The selected guides are computational candidates, not validated screening reagents.
2. **Public-history rewrite:** the owner requested a clean single-root public branch to remove preliminary CSVs and the historical workbook from reachable `main` history. This does not guarantee deletion of external clones or caches. Verify the remote after the force update and do not reintroduce old commits.
3. **Future claims:** a public row-level dataset requires file-level rights review; a mature end-to-end workflow requires the remaining Level B work; a screening-ready library requires Level C and experiments.
4. **Sign-off:** any force update to the public branch requires an exact file/tree review and the owner's explicit final approval. This document alone authorizes no push, Zenodo publication or DOI assignment.
