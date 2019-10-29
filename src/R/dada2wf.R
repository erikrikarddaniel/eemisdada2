#!/usr/bin/env Rscript

# dada2wf.R
#
# Runs DADA2 denoising, from trimming all the way through to bimera filtering.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

SCRIPT_VERSION = "0.9"

# Get arguments
option_list = list(
  make_option(
    c('--bimerafilename_prefix'), type='character', default='dada2.cleaned.merged.bimera',
    help='Prefix for output files from the bimera check, default %default.'
  ),
  make_option(
    c("--concatenate"), action="store_true", default=FALSE, 
    help="Concatenate sequences instead of try to merge (use when overlap is too short), default %default. Setting this overrides --maxmismatch and --minoverlap."
  ),
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
    c('--fwdmark'), type='character', default='.r1',
    help='Name pattern to identify forward reads in file names, default ".r1".'
  ),
  make_option(
    c('--maxconsist'), type='integer', default=10,
    help='Maximum number of iterations to achieve self consistency, default 10.'
  ),
  make_option(
    c('--maxee'), type='character', default=NA,
    help='Maximum number of expected errors in sequences, comma-separated pair for forward and reverse, default Inf.'
  ),
  make_option(
    c('--maxmismatch'), type='integer', default=0,
    help="Maximum number of non matching bases in merge, default %default."
  ),
  make_option(
    c('--mergefile_prefix'), type='character', default='dada2.cleaned.merged',
    help='Prefix for denoised and merge sequence output files, default %default.'
  ),
  make_option(
    c('--method'), type='character', default='consensus',
    help='Method used for bimera detection, default "%default". See R doc for removeBimeraDenovo.'
  ),
  make_option(
    c('--minab'), type='integer', default=8,
    help='Minimum parent abundance, default %default. See R doc for isBimeraDenovo.'
  ),
  make_option(
    c('--minoverlap'), type='integer', default=20,
    help="Minimum number of overlapping bases in merge, default %default."
  ),
  make_option(
    c('--minq'), type='integer', default=0,
    help='minQ option in fastqPairedFilter() call. Default 0.'
  ),
  make_option(
    c('--nsamples'), type='integer', default=0,
    help='Number of samples to select for learning algorithm. Should be high enough to make sure at least 1 million reads are included. If 0, an eigth of the number of samples will be used, assuming you have at least 10 million reads in total, fairly evenly distributed over samples. Default 0.'
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
    c('--overab'), type='integer', default=2,
    help='Parent overabundance multiplier, default %default. See R doc for isBimeraDenovo.'
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
    c("--rmphix"), action="store_true", default=FALSE, 
    help="Remove phiX sequences. Currently fails (recursion in function call?)"
  ),
  make_option(
    c("--seed"), default = '0', 
    help = "Set seed for random number generation, default don't set."
  ),
  make_option(
    c('--seqfile_prefix'), type='character', default='seq',
    help='String to prefix output files with, default "seq".'
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
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] \n\n\tRuns the DADA2 denoising workflow.",
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
  quit('no')
}

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

# I can't get type = 'integer' to work above...
opt$options$seed <- as.integer(opt$options$seed)
if ( opt$options$seed > 0 ) set.seed(opt$options$seed)

infiles = Sys.glob(sprintf("%s/*.fastq.gz", opt$options$indir))
if ( length(infiles) == 0 ) {
  logmsg("No files found, exiting", 'WARNING')
  q('no', 0)
}
logmsg(sprintf("You have %d files in %s that are going to be cleaned.", length(infiles), opt$options$indir))

fwdinfiles = infiles[grep(opt$options$fwdmark, infiles, fixed=T)]
revinfiles = infiles[grep(opt$options$revmark, infiles, fixed=T)]

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

# The DADA2 library takes long to load, do here to avoid --help and --version take so long...
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

samples = sub(sprintf("%s.*", opt$options$fwdmark), "", basename(fwdinfiles))

# Cutting sequences
fwdfiltered = file.path(opt$options$filterdir, paste0(samples, opt$options$fwdmark, ".fastq.gz"))
revfiltered = file.path(opt$options$filterdir, paste0(samples, opt$options$revmark, ".fastq.gz"))

if ( ! file_test("-d", opt$options$filterdir) ) dir.create(opt$options$filterdir)

if ( ! is.na(opt$options$trunclen) ) {
  trunclen = as.integer(unlist(strsplit(opt$options$trunclen, split=',')))
} else {
  stop("You need to set --trunclen. See help (-h).")
}
if ( ! is.na(opt$options$trimleft) ) {
  trimleft = as.integer(unlist(strsplit(opt$options$trimleft, split=',')))
} else {
  stop("You need to set --trimleft. See help (-h).")
}
if ( ! is.na(opt$options$maxee) ) {
  maxee = as.integer(unlist(strsplit(opt$options$maxee, split=',')))
} else {
  maxee = c(Inf, Inf)
}

logmsg(
  sprintf(
    "Running quality filtering, truncLen: %s, trimLeft: %s, maxN: %d, maxEE: %s, truncQ: %d, minQ: %d, rm.phix: %s", 
    paste(trunclen, collapse=','), paste(trimleft, collapse=','), 0, paste(maxee, collapse=','), opt$options$truncq, opt$options$minq, opt$options$rmphix
  )
)
for ( i in seq_along(fwdinfiles) ) {
  logmsg(sprintf("--> Filtering %s <--", samples[i]))
  fastqPairedFilter(
    c(fwdinfiles[i], revinfiles[i]),
    c(fwdfiltered[i], revfiltered[i]),
    truncLen = trunclen, trimLeft = trimleft,
    maxN = 0, maxEE = maxee, truncQ = opt$options$truncq, minQ = opt$options$minq,
    rm.phix = opt$options$rmphix, compress = T, verbose = opt$options$verbose
  )
}

logmsg(
  sprintf(
    "Starting estimation of error models (forward: %d, reverse: %d) from files in %s. Output files will be prefixed with '%s'.", 
    opt$options$fwd, opt$options$rev, opt$options$filterdir, opt$options$seqfile_prefix
  )
)

filtered = Sys.glob(sprintf("%s/*.fastq.gz", opt$options$filterdir))

fwdfiltered = filtered[grep(opt$options$fwdmark, filtered, fixed=T)]
revfiltered = filtered[grep(opt$options$revmark, filtered, fixed=T)]

samples = sub(sprintf("%s.*", opt$options$fwdmark), "", basename(fwdfiltered))

names(fwdfiltered) = samples
names(revfiltered) = samples

# Learning the error models
if ( opt$options$nsamples == 0 ) {
  nsamples = length(samples) %/% 8 + 1
} else if ( opt$options$nsamples > length(samples) ) {
  logmsg(
    sprintf(
      "You asked for %d samples to estimate error rates from, but only %d present, using %d",
      opt$options$nsamples, length(samples), length(samples)
    )
  )
  nsamples = length(samples)
} else {
  nsamples = opt$options$nsamples
}

if ( opt$options$fwd ) {
  logmsg(sprintf("Learning forward error rates, using %d samples.", nsamples))

  fwdlearn = dada(
    derepFastq(sample(fwdfiltered, nsamples)),
    err=NULL, selfConsist=TRUE, multithread=TRUE,
    MAX_CONSIST = opt$options$maxconsist
  )
  fwderr = fwdlearn[[1]]$err_out
  saveRDS(fwderr, sprintf('%s.dada2errmodels.fwd.errorates.rds', opt$options$seqfile_prefix))
}

if ( opt$options$rev ) {
  logmsg(sprintf("Learning reverse error rates, using %d samples.", nsamples))

  revlearn = dada(
    derepFastq(sample(revfiltered, nsamples)),
    err=NULL, selfConsist=TRUE, multithread=TRUE,
    MAX_CONSIST = opt$options$maxconsist
  )
  reverr = revlearn[[1]]$err_out
  saveRDS(reverr, sprintf('%s.dada2errmodels.rev.errorates.rds', opt$options$seqfile_prefix))
}

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
  if ( opt$options$concatenate ) {
    logmsg("	* concatenating *")
    merge <- mergePairs(fwddada2, fwdderep, revdada2, revderep, justConcatenate = TRUE)
  } else {
    logmsg("	merging")
    merge = mergePairs(
      fwddada2, fwdderep, revdada2, revderep, 
      minOverlap=opt$options$minoverlap, maxMismatch=opt$options$maxmismatch
    )
  }
  mergers[[s]] = merge
}

logmsg("Constructing sequence table and removing chimeras")
seqtab = makeSequenceTable(mergers)

# Turn into a standard, long, data frame
seqtab.df = data.frame(t(seqtab)) %>% tibble::rownames_to_column('seq')
seqtab.dfl = seqtab.df %>% gather(sample, count, 2:length(colnames(seqtab.df))) %>% filter(count>0)
write_tsv(seqtab.dfl, sprintf("%s.tsv.gz", opt$options$mergefile_prefix))

seqtab.sum = seqtab.dfl %>% summarise(count=sum(count))
logmsg(sprintf("%d sequences remaining, totalling %d observations", length(seqtab.df$seq), seqtab.sum$count[1]))

logmsg(sprintf("Removing bimeras, overab: %d, minab: %d, oneoff: %s", opt$options$overab, opt$options$minab, opt$options$oneoff))
seqtab = removeBimeraDenovo(
  seqtab, multithread=TRUE,
  minFoldParentOverAbundance=opt$options$overab,
  minParentAbundance = opt$options$minab, 
  allowOneOff=opt$options$oneoff,
  verbose = opt$options$verbose
)

# Turn into a standard, long, data frame
seqtab.df = data.frame(t(seqtab)) %>% tibble::rownames_to_column('seq')
seqtab.dfl = seqtab.df %>% gather(sample, count, 2:length(colnames(seqtab.df))) %>% filter(count>0)

logmsg(sprintf("Writing tsv file %s.tsv.gz", opt$options$bimerafilename_prefix))
write_tsv(seqtab.dfl, sprintf("%s.tsv.gz", opt$options$bimerafilename_prefix))

seqtab.sum = seqtab.dfl %>% summarise(count=sum(count))
logmsg(sprintf("%d sequences remaining, totalling %d observations", length(seqtab.df$seq), seqtab.sum$count[1]))

logmsg("Done")
