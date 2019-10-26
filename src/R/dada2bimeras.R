#!/usr/bin/env Rscript

# dada2bimeras
#
# Runs the DADA2 bimera filter.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))

VERSION = sprintf("dada2bimeras version 1.0.1, DADA2 version: %s", packageDescription('dada2')$Version)

# Get arguments
option_list = list(
  make_option(
    c('--method'), type='character', default='consensus',
    help='Method used for bimera detection, default "%default". See R doc for removeBimeraDenovo.'
  ),
  make_option(
    c('--minab'), type='integer', default=8,
    help='Minimum parent abundance, default %default. See R doc for isBimeraDenovo.'
  ),
  make_option(
    c("--oneoff"), action="store_true", default=TRUE,
    help='Allow one off bimeras, default %default. See R doc for isBimeraDenovo.'
  ),
  make_option(
    c("--nooneoff"), action="store_false", dest='oneoff',
    help='Do not allow one off bimeras. See R doc for isBimeraDenovo.'
  ),
  make_option(
    c('--overab'), type='integer', default=1,
    help='Parent overabundance multiplier, default %default. See R doc for isBimeraDenovo.'
  ),
  make_option(
    c('--prefix'), type='character', default='dada2.cleaned.merged.bimera',
    help='Prefix for output files, default %default.'
  ),
  make_option(
    c("--seed"), default = '0', 
    help = "Set seed for random number generation, default don't set."
  ),
  make_option(
    c('--seqtabfile'), type='character', default='dada2.cleaned.merged.rds',
    help='RDS file with seqtab object from dada2cleanNmerge, default %default.'
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

# Read sequence table object from cleanNmerge step

logmsg(sprintf("Reading seqtab from %s", opt$seqtabfile))
seqtab = readRDS(opt$seqtabfile)

logmsg(sprintf("Removing bimeras, overab: %d, minab: %d, oneoff: %s", opt$overab, opt$minab, opt$oneoff))
seqtab = removeBimeraDenovo(
  seqtab, multithread=TRUE,
  minFoldParentOverAbundance=opt$overab,
  minParentAbundance = opt$minab, 
  allowOneOff=opt$oneoff,
  verbose = opt$verbose
)

logmsg(sprintf("Saving to rds file %s.rds", opt$prefix))
saveRDS(seqtab, sprintf("%s.rds", opt$prefix))

# Turn into a standard, long, data frame
seqtab.df = data.frame(t(seqtab)) %>% tibble::rownames_to_column('seq')
seqtab.dfl = seqtab.df %>% gather(sample, count, 2:length(colnames(seqtab.df))) %>% filter(count>0)

logmsg(sprintf("Writing tsv file %s.tsv.gz", opt$prefix))
write_tsv(seqtab.dfl, sprintf("%s.tsv.gz", opt$prefix))

seqtab.sum = seqtab.dfl %>% summarise(count=sum(count))
logmsg(sprintf("%d sequences remaining, totalling %d observations", length(seqtab.df$seq), seqtab.sum$count[1]))

logmsg("Done")
