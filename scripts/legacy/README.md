# Superseded coding-effect implementations

Use [the current coding-effect script](../predict_AtoG_coding_effects.py) for corrected analyses; files here are retained only to document historical calculations.

- [predict_AtoG_coding_effects.py](predict_AtoG_coding_effects.py): original bedtools/samtools implementation, superseded because its negative-strand coding effects were unreliable.
- [predict_AtoG_coding_effects_v2_historical.py](predict_AtoG_coding_effects_v2_historical.py): added RefSeq FASTA-contig resolution but still retained the negative-strand codon-orientation error.
