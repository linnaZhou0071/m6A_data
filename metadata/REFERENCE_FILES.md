# Reference files used locally

These are provenance records for local inputs, not redistributed copies.
SHA-256 values identify the exact local bytes. Download dates were not
recorded at acquisition and must not be inferred from file timestamps.

| Local file | Workflow use | Source directory / version | SHA-256 |
| --- | --- | --- | --- |
| `GCF_000001405.26_GRCh38_genomic.fna.gz` | Step 9: PEGG genome input; compressed source for step 6 | [NCBI assembly GCF_000001405.26, GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/) | `1d77c0aea04dcdcb9a4fb17ddd811a18288ceeeb95f61f3574c71627ba671d4e` |
| `GRCh38.fa` | Step 6: fetch reference bases to reconstruct RNA codons | Byte-for-byte decompressed copy of the preceding NCBI FASTA; verified by streaming decompression into a byte comparison | `cfdd7bbdfddd25b8cf32a54787a360cc7ede3ee01ef7ab5c629a970ed6ef2cc4` |
| `hg38.ncbiRefSeq.gtf.gz` | Steps 4–5: CDS/UTR annotation and transcript mapping; also used with FASTA in step 6 | [UCSC hg38 genes download](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/); server listing dates this file 2022-10-28 | `856919cfc5854079e70dd016048045092fd79b782aa8da9dbbd1c51a9046d8a4` |
| `MANE.GRCh38.v1.5.refseq_genomic.gtf.gz` | Step 7: identify MANE Select transcript IDs for the exclusion rule | [NCBI MANE human release 1.5](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/); RefSeq genomic GTF | `59e5ff1c794a70de777e553ab67062b316578d41287cd898ee1da2827527969a` |

The NCBI FASTA uses RefSeq contig accessions such as `NC_000001.11`,
while the UCSC GTF and site IDs use `chr1`-style names. The effect
script maps these names before fetching codons. The PEGG input script
uses chromosome numbers without `chr`, as required by its selected
cBioPortal input format. No coordinate liftover is implied.

Cite the MANE resource as:
[Morales, Pujar, Loveland et al., Nature 604, 310–315 (2022),
DOI 10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8).
That article describes the MANE project; release 1.5 is identified
separately by the exact NCBI directory and filename above.

The reference FASTA, annotation GTFs and their indexes are third-party
files, kept locally and excluded from this repository's license grants.
