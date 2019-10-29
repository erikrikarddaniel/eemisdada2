#!/usr/bin/env Rscript

# dada2filter
#
# Filters sequences using the fastqPairedFilter function in the dada2 package.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))

VERSION = "1.0.1"

# Get arguments
option_list = list(
  make_option(
    c('--indir'), type='character', default='.',
    help='Input directory, default current directory.'
  ),
  make_option(
    c('--filterdir'), type='character', default='filtered',
    help='Directory for quality truncated reads, default "filtered". Will be created if it does not exist.'
  ),
  make_option(
    c('--fwdmark'), type='character', default='.r1',
    help='Name pattern to identify forward reads in file names, default ".r1".'
  ),
  make_option(
    c('--maxee'), type='character', default=NA,
    help='Maximum number of expected errors in sequences, comma-separated pair for forward and reverse, default Inf.'
  ),
  make_option(
    c('--minq'), type='integer', default=0,
    help='minQ option in fastqPairedFilter() call. Default 0.'
  ),
  make_option(
    c('--revmark'), type='character', default='.r2',
    help='Name pattern to identify reverse reads in file names, default ".r2".'
  ),
  make_option(
    c("--rmphix"), action="store_true", default=FALSE, 
    help="Remove phiX sequences. Currently fails (recursion in function call?)"
  ),
  make_option(
    c('--trimleft'), type='character', default="0,0",
    help='Comma separated pair of values for left truncation of forward and reverse reads respectively. Default "0,0".'
  ),
  make_option(
    c('--trunclen'), type='character', default=NA,
    help='Comma separated pair of values for truncation position for forward and reverse reads respectively. No default, mandatory'
  ),
  make_option(
    c('--truncq'), type='integer', default=2,
    help='truncQ option in fastqPairedFilter() call. Default 2.'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print version of script and DADA2 library"
  )
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$version ) {
  write(sprintf("dada2filter version %s, DADA2 version %s", VERSION, packageVersion('dada2')), stdout())
  q('no', 0)
}

logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

infiles = Sys.glob(sprintf("%s/*.fastq.gz", opt$indir))
if ( length(infiles) == 0 ) {
  print("No input files, exiting")
  q('no', 1)
}
logmsg(sprintf("You have %d files in %s that are going to be cleaned.", length(infiles), opt$indir))

fwdinfiles = infiles[grep(opt$fwdmark, infiles, fixed=T)]
revinfiles = infiles[grep(opt$revmark, infiles, fixed=T)]

if ( length(fwdinfiles) != length(revinfiles) ) {
  logmsg(
    sprintf(
      "You do not have the same number (%d) of forward as reverse (%d) files. Exiting.", 
      length(fwdinfiles), length(revinfiles)
    ), 
    llevel="WARNING"
  )
  q('no', 1)
}

# The DADA2 library takes long to load, do here to avoid --help and --version takes so long...
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

samples = sub(sprintf("%s.*", opt$fwdmark), "", basename(fwdinfiles))

# Cutting sequences
fwdfiltered = file.path(opt$filterdir, paste0(samples, opt$fwdmark, ".fastq.gz"))
revfiltered = file.path(opt$filterdir, paste0(samples, opt$revmark, ".fastq.gz"))

if ( ! file_test("-d", opt$filterdir) ) dir.create(opt$filterdir)

if ( ! is.na(opt$trunclen) ) {
  trunclen = as.integer(unlist(strsplit(opt$trunclen, split=',')))
} else {
  stop("You need to set --trunclen. See help (-h).")
}
if ( ! is.na(opt$trimleft) ) {
  trimleft = as.integer(unlist(strsplit(opt$trimleft, split=',')))
} else {
  stop("You need to set --trimleft. See help (-h).")
}
if ( ! is.na(opt$maxee) ) {
  maxee = as.integer(unlist(strsplit(opt$maxee, split=',')))
} else {
  maxee = c(Inf, Inf)
}

logmsg(
  sprintf(
    "Running quality filtering, truncLen: %s, trimLeft: %s, maxN: %d, maxEE: %s, truncQ: %d, minQ: %d, rm.phix: %s", 
    paste(trunclen, collapse=','), paste(trimleft, collapse=','), 0, paste(maxee, collapse=','), opt$truncq, opt$minq, opt$rmphix
  )
)
for ( i in seq_along(fwdinfiles) ) {
  logmsg(sprintf("--> Filtering %s <--", samples[i]))
  fastqPairedFilter(
    c(fwdinfiles[i], revinfiles[i]),
    c(fwdfiltered[i], revfiltered[i]),
    truncLen = trunclen, trimLeft = trimleft,
    maxN = 0, maxEE = maxee, truncQ = opt$truncq, minQ = opt$minq,
    rm.phix = opt$rmphix, compress = T, verbose = opt$verbose
  )
}

logmsg(sprintf("Done, output in %s", opt$filterdir))
