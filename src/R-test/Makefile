# Makefile for automatic and manual tests.
#
# First, download testdata in ../../testdata: make all
#
# Run test for individual steps using the targets here
#
DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

dada2filter: dada2filter.00

dada2taxonomy: dada2taxonomy.00

dada2filter.00:
	../R/dada2filter.R --verbose --indir ../../testdata/reads/ --fwdmark _R1_ --revmark _R2_ --trimleft 13,13 --trunclen 150,150
	@ls -l filtered/ | grep -v 'total' | sed 's/^[^ ]* [0-9] [^ ]* [^ ]* \( *[0-9]*\) ... [0-9]\+ ..:.. \(.*\)/\1 \2/' > $@.out
	@$(CHECK)

dada2taxonomy.00:
	../R/dada2taxonomy.R --verbose --seed 42 \
	  --rdp_fasta ../../testdata/silva.fna.gz \
	  --idtaxa_rdata ../../testdata/SILVA_SSU_r132_March2018.RData \
	  --species_fasta ../../testdata/silva_species.fna.gz ../../testdata/fifty_16S_sequences.fna
	@gunzip -c fifty_16S_sequences.idtaxa.silva.tsv.gz fifty_16S_sequences.rdp.silva.tsv.gz > $@.out
	@$(CHECK)
