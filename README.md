# eemisdada2
Packaging of dada2 for EEMiS, and other, users.

[DADA2](http://benjjneb.github.io/dada2/) is an error correction algorithm for
Illumina amplicon reads.  This repository offers three R scripts that lets you
run the three main steps in the processing separately, saving intermediate
results: 

1. `dada2filter`: quality and length trimming
2. `dada2errmodels`: calculation of error models from the data
3. `dada2correction`: sequence correction, bimera filtering and table creation

In my [`biomakefiles`
repository](https://github.com/erikrikarddaniel/biomakefiles), there's a
makefile-based workflow implementation with documentation.
