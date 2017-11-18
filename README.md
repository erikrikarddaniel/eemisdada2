# eemisdada2
Packaging of dada2 for EEMiS, and other, users.

[DADA2](http://benjjneb.github.io/dada2/) is an error correction algorithm for
Illumina amplicon reads.  This repository offers four R scripts that lets you
run the four main steps (cleaning and merging is performed by the same script)
in the processing separately, saving intermediate results: 

1. `dada2filter`: quality and length trimming
2. `dada2errmodels`: calculation of error models from the data
3. `dada2cleanNmerge`: sequence correction and merge
4. `dada2bimeras`: filter bimeras, create final tables

DADA2 can be installed from the 
[Bioconductor](http://bioconductor.org/packages/release/bioc/html/dada2.html)
site.

DADA2 requires headers for the `libcurl4-gnutls` library. In Debian
they can be installed with this command:

```bash
$ sudo apt install libcurl4-gnutls-dev
```

In my [`biomakefiles` repository](https://github.com/erikrikarddaniel/biomakefiles),
there's a makefile-based workflow implementation with documentation.
