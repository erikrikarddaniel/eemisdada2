#!/usr/bin/env Rscript

# run_dada2
#
# This is the second script in a series of three that implements a DADA2-based
# pipeline for correction of amplicon reads.
#
# Runs the DADA2 sequence error estimation algorithm on all fastq.gz files in a
# specified directory. The files are presumed to be filtered and cut using the
# dada2filter program.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))

VERSION = sprintf("dada2errmodels version 1.0, DADA2 version: %s", packageDescription('dada2')$Version)

# Get arguments
option_list = list(
  make_option(
    c('--filterdir'), type='character', default='filtered',
    help='Directory for quality truncated reads, default "filtered". Will be created if it does not exist.'
  ),
  make_option(
    c('--fwd'), action='store_true', default=TRUE,
    help='Calculate matrix for forward reads, default %default.'
  ),
  make_option(
    c('--nofwd'), action='store_false', dest='fwd',
    help='Do not calculate matrix for forward reads.'
  ),
  make_option(
    c('--fwdmark'), type='character', default='.r2',
    help='Name pattern to identify forward reads in file names, default "%default".'
  ),
  make_option(
    c('--maxconsist'), type='integer', default=10,
    help='Maximum number of iterations to achieve self consistency, default 10.'
  ),
  make_option(
    c('--nsamples'), type='integer', default=0,
    help='Number of samples to select for learning algorithm. Should be high enough to make sure at least 1 million reads are included. If 0, an eigth of the number of samples will be used, assuming you have at least 10 million reads in total, fairly evenly distributed over samples. Default 0.'
  ),
  make_option(
    c('--prefix'), type='character', default='seq',
    help='String to prefix output files with, default "seq".'
  ),
  make_option(
    c('--rev'), action='store_true', default=TRUE,
    help='Calculate matrix for reverse reads, default %default.'
  ),
  make_option(
    c('--norev'), action='store_false', dest='rev',
    help='Do not calculate matrix for reverse reads.'
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
  write(VERSION, stderr())
  q('no')
}

# Load DADA2 here, takes a long while...
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

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

logmsg(
  sprintf(
    "Starting estimation of error models (forward: %d, reverse: %d) from files in %s. Output files will be prefixed with '%s'.", 
    opt$fwd, opt$rev, opt$filterdir, opt$prefix
  )
)

filtered = Sys.glob(sprintf("%s/*.fastq.gz", opt$filterdir))

fwdfiltered = filtered[grep(opt$fwdmark, filtered, fixed=T)]
revfiltered = filtered[grep(opt$revmark, filtered, fixed=T)]

samples = sub(sprintf("%s.*", opt$fwdmark), "", basename(fwdfiltered))

names(fwdfiltered) = samples
names(revfiltered) = samples

# Learning the error models
if ( opt$nsamples == 0 ) {
  nsamples = length(samples) %/% 8 + 1
} else {
  nsamples = opt$nsamples
}

if ( opt$fwd ) {
  logmsg(sprintf("Learning forward error rates, using %d samples.", nsamples))

  fwdlearn = dada(
    derepFastq(sample(fwdfiltered, nsamples)),
    err=NULL, selfConsist=TRUE, multithread=TRUE,
    MAX_CONSIST = opt$maxconsist
  )
  saveRDS(fwdlearn, sprintf('%s.dada2errmodels.fwd.rds', opt$prefix))
  fwderr = fwdlearn[[1]]$err_out
  saveRDS(fwderr, sprintf('%s.dada2errmodels.fwd.errorates.rds', opt$prefix))
}

if ( opt$rev ) {
  logmsg(sprintf("Learning reverse error rates, using %d samples.", nsamples))

  revlearn = dada(
    derepFastq(sample(revfiltered, nsamples)),
    err=NULL, selfConsist=TRUE, multithread=TRUE,
    MAX_CONSIST = opt$maxconsist
  )
  saveRDS(revlearn, sprintf('%s.dada2errmodels.rev.rds', opt$prefix))
  reverr = revlearn[[1]]$err_out
  saveRDS(reverr, sprintf('%s.dada2errmodels.rev.errorates.rds', opt$prefix))
}

logmsg("Done")
