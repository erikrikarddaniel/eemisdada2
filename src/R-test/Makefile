# Makefile for automatic and manual tests.
#
# First, download testdata in ../../testdata: 'make all'
#
# Run test for individual steps using the targets here
#
DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: dada2filter dada2errmodels dada2cleanNmerge dada2bimeras dada2idseq dada2taxonomy

dada2filter: dada2filter.00

dada2errmodels: dada2errmodels.00

dada2cleanNmerge: dada2cleanNmerge.00

dada2bimeras: dada2bimeras.00

dada2idseq: dada2idseq.00 dada2idseq.01

dada2taxonomy: dada2taxonomy.00

dada2filter.00:
	../R/dada2filter.R --verbose --indir ../../testdata/reads/ --fwdmark _R1_ --revmark _R2_ --trimleft 13,13 --trunclen 150,150
	@ls -l filtered/ | grep -v 'total' | sed 's/^[^ ]* [0-9] [^ ]* [^ ]* \( *[0-9]*\) ... [0-9]\+ ..:.. \(.*\)/\1 \2/' > $@.out
	@$(CHECK)

dada2errmodels.00:
	@if [ ! -e filtered ] && [ $$(ls filtered/*.fastq.gz 2>/dev/null | wc -l) -eq 0 ]; then echo "*** You must run 'make dada2filter' (or try 'make all') first to create input data ***"; exit 1; fi
	../R/dada2errmodels.R --verbose --seed=42 --filterdir=filtered --fwdmark=_R1_ --revmark=_R2_ --nsamples=10 
	@ls -l *.errorates.rds | sed 's/^[^ ]* [0-9] [^ ]* [^ ]* \( *[0-9]*\) ... [0-9]\+ ..:.. \(.*\)/\1 \2/' > $@.out
	@$(CHECK)

dada2cleanNmerge.00:
	@if [ ! -e seq.dada2errmodels.fwd.errorates.rds ] || [ ! -e seq.dada2errmodels.rev.errorates.rds ]; then echo "*** You must run 'make dada2errmodels' (or try 'make all') to create input data ***"; exit 1; fi
	../R/dada2cleanNmerge.R --verbose --filterdir=filtered --fwdmark=_R1_ --revmark=_R2_ --fwderrmodel=seq.dada2errmodels.fwd.errorates.rds --reverrmodel=seq.dada2errmodels.rev.errorates.rds
	@ls -l dada2.cleaned.merged.tsv.gz | sed 's/^[^ ]* [0-9] [^ ]* [^ ]* \( *[0-9]*\) ... [0-9]\+ ..:.. \(.*\)/\1 \2/' > $@.out
	@$(CHECK)

dada2bimeras.00:
	@if [ ! -e dada2.cleaned.merged.rds ]; then echo "*** You must run 'make dada2cleanNmerge' (or try 'make all') to create input data ***"; exit 1; fi
	../R/dada2bimeras.R --verbose --method=consensus --minab=4 --seqtabfile=dada2.cleaned.merged.rds
	@ls -l dada2.cleaned.merged.bimera.tsv.gz | sed 's/^[^ ]* [0-9] [^ ]* [^ ]* \( *[0-9]*\) ... [0-9]\+ ..:.. \(.*\)/\1 \2/' > $@.out
	@$(CHECK)

dada2idseq.00:
	@rm -f $@.out
	@../R/dada2idseq.R --verbose --prefix=S_ --fnafile=$@.out --outtable=$@.out.tsv.gz $@.tsv.gz
	@gunzip -c $@.out.tsv.gz >> $@.out
	@$(CHECK)

dada2idseq.01:
	@cp $@.fna $@.out
	@../R/dada2idseq.R --verbose --prefix=S_ --fnafile=$@.out --outtable=$@.out.tsv.gz $@.tsv.gz
	@gunzip -c $@.out.tsv.gz >> $@.out
	@$(CHECK)

dada2taxonomy.00:
	@../R/dada2taxonomy.R --verbose --seed 42 \
	  --rdp_fasta ../../testdata/silva.fna.gz \
	  --idtaxa_rdata ../../testdata/SILVA_SSU_r132_March2018.RData \
	  --species_fasta ../../testdata/silva_species.fna.gz ../../testdata/fifty_16S_sequences.fna
	@gunzip -c fifty_16S_sequences.idtaxa.silva.tsv.gz fifty_16S_sequences.rdp.silva.tsv.gz > $@.out
	@$(CHECK)
