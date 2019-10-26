#!/usr/bin/env Rscript

# dada2cleanNmerge
#
# Runs the DADA2 sequence cleaning algorithm on all fastq.gz files in a
# specified directory using error models read from specified rds files. The
# files are presumed to be filtered and cut using the dada2filter program.
#
# After cleaning, pairs are merged.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))

VERSION = sprintf("dada2cleanNmerge version 1.1.1, DADA2 version: %s", packageDescription('dada2')$Version)

# Get arguments
option_list = list(
  make_option(
    c("--concatenate"), action="store_true", default=FALSE, 
    help="Concatenate sequences instead of try to merge (use when overlap is too short), default %default. Setting this overrides --maxmismatch and --minoverlap."
  ),
  make_option(
    c('--filterdir'), type='character', default='filtered',
    help='Directory for quality truncated reads, default "filtered". Will be created if it does not exist.'
  ),
  make_option(
    c('--fwderrmodel'), type='character', 
    help='R RDS file with error model for forward reads, as output by the dada2errmodels script (NNN.fwd.errorates.rds), no default.'
  ),
  make_option(
    c('--fwdmark'), type='character', default='.r1',
    help='Name pattern to identify forward reads in file names, default ".r1".'
  ),
  make_option(
    c('--maxmismatch'), type='integer', default=0,
    help="Maximum number of non matching bases in merge, default %default."
  ),
  make_option(
    c('--minoverlap'), type='integer', default=20,
    help="Minimum number of overlapping bases in merge, default %default."
  ),
  make_option(
    c('--prefix'), type='character', default='dada2.cleaned.merged',
    help='Prefix for output files, default %default.'
  ),
  make_option(
    c('--reverrmodel'), type='character', 
    help='R RDS file with error model for reverse reads, as output by the dada2errmodels script (NNN.rev.errorates.rds), no default.'
  ),
  make_option(
    c('--revmark'), type='character', default='.r2',
    help='Name pattern to identify reverse reads in file names, default ".r2".'
  ),
  make_option(
    c("--seed"), default = '0', 
    help = "Set seed for random number generation, default don't set."
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("--version"), action="store_true", default=FALSE, 
    help="Show script version plus installed DADA2 version."
  )
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$version ) {
  write(VERSION, stdout())
  q('no')
}

# Load packages here so they don't slow down --help and --version
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(tidyr))

logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

# I can't get type = 'integer' to work above...
opt$seed <- as.integer(opt$seed)
if ( opt$seed > 0 ) set.seed(opt$seed)

filtered = Sys.glob(sprintf("%s/*.fastq.gz", opt$filterdir))

fwdfiltered = filtered[grep(opt$fwdmark, filtered, fixed=T)]
revfiltered = filtered[grep(opt$revmark, filtered, fixed=T)]

samples = sub(sprintf("%s.*", opt$fwdmark), "", basename(fwdfiltered))

names(fwdfiltered) = samples
names(revfiltered) = samples

logmsg(sprintf("Reading error models from %s (forward) and %s (reverse).", opt$fwderrmodel, opt$reverrmodel))
fwderr = readRDS(opt$fwderrmodel)
reverr = readRDS(opt$reverrmodel)

logmsg("Applying the error models")
mergers = vector('list', length(samples))
names(mergers) = samples
for ( s in samples ) {
  logmsg(sprintf("--> Cleaning %s <--", s))
  logmsg("	forward")
  fwdderep = derepFastq(fwdfiltered[[s]])
  fwddada2 = dada(fwdderep, err=fwderr, multithread=T)
  logmsg("	reverse")
  revderep = derepFastq(revfiltered[[s]])
  revdada2 = dada(revderep, err=reverr, multithread=T)
  if ( opt$concatenate ) {
    logmsg("	* concatenating *")
    merge <- mergePairs(fwddada2, fwdderep, revdada2, revderep, justConcatenate = TRUE)
  } else {
    logmsg("	merging")
    merge = mergePairs(
      fwddada2, fwdderep, revdada2, revderep, 
      minOverlap=opt$minoverlap, maxMismatch=opt$maxmismatch
    )
  }
  mergers[[s]] = merge
}

logmsg("Constructing sequence table and removing chimeras")
seqtab = makeSequenceTable(mergers)
saveRDS(seqtab, sprintf("%s.rds", opt$prefix))

# Turn into a standard, long, data frame
seqtab.df = data.frame(t(seqtab)) %>% tibble::rownames_to_column('seq')
seqtab.dfl = seqtab.df %>% gather(sample, count, 2:length(colnames(seqtab.df))) %>% filter(count>0)
write_tsv(seqtab.dfl, sprintf("%s.tsv.gz", opt$prefix))

seqtab.sum = seqtab.dfl %>% summarise(count=sum(count))
logmsg(sprintf("%d sequences remaining, totalling %d observations", length(seqtab.df$seq), seqtab.sum$count[1]))

logmsg("Done")
