#!/usr/bin/env Rscript

# dada2fasta.R
#
# Reads a table with sequences, such as those outtable by dada2cleanNmerge.R and
# dada2bimeras.R, optionally combines with an fasta file with unique sequences,
# and outtables a new fasta file. Each sequence will be given a unique name.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

SCRIPT_VERSION = "0.2.1"

# Get arguments
# For testing: opt <- list(args = 'dada2idseq.00.tsv.gz', options = list(fnafile = 'dada2idseq.01.fna', idlen = 9, prefix = 'S_'))
option_list = list(
  make_option(
    c('--fnafile'), type='character', default = 'sequences.fna.gz',
    help='Name of fasta file, existing or not, that will be read and written 
		to; default %default. Sequence names has to start with the 
    		same prefix as specified by the prefix option here.'
  ),
  make_option(
    c('--idlen'), type = 'integer', default = 9,
    help = 'Number of characters to use for numerical index in sequence names,
    		default %default.'
  ),
  make_option(
    c('--outtable'), type='character', default = 'asvtable.tsv.gz',
    help='Name of outtable tsv file, default "%default".'
  ),
  make_option(
    c('--prefix'), type='character', default = 'seq_',
    help='Prefix to use for sequence names, default "%default".'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages."
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] seqtable.tsv[.gz]

    Reads seqtable.tsv[.gz] (long format: seq, sample, count) and, if 
    existing, a fasta file with old sequences, and outputs a fasta file with
    unique sequences plus a new table with sequence names instead of
    sequences.",
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

# Read the fasta file we might be called with
if ( file.exists(opt$options$fnafile) ) {
  logmsg(sprintf("Reading fasta file %s", opt$options$fnafile))
  seqs <- data.frame(seq = readDNAStringSet(opt$options$fnafile)) %>%
    tibble::rownames_to_column('seqname') %>%
    mutate(seqnum = sub(opt$options$prefix, '', seqname) %>% as.integer())
  max_seqnum = ifelse(
    nrow(seqs) > 0, 
    seqs %>% filter(seqnum == max(seqnum)) %>% pull(seqnum),
    0
  )
} else {
  seqs <- tibble(seqname = character(), seq = character(), seqnum = integer())
  max_seqnum = 0
}
#logmsg(sprintf("Read %s, max_seqnum: %d, colnames: %s", opt$options$fnafile, max_seqnum, paste(colnames(seqs), collapse = ', ')), 'DEBUG')

# Read the ASV table, assumed to be in long format
logmsg(sprintf("Read ASV table %s", opt$args[1]))
seqtab <- read_tsv(
  opt$args[1],
  col_types = cols(.default = col_character(), count = col_integer())
)

# Join unique sequences from the two tables, create new seqname for those
# missing from input fna
seqname_format = sprintf("%%s%%0%dd", opt$options$idlen)
seqs <- seqs %>% select(-seqnum) %>%
  union(
    seqtab %>% distinct(seq) %>%
      anti_join(seqs, by = 'seq') %>%
      mutate(seqname = sprintf(seqname_format, opt$options$prefix, max_seqnum + rank(seq)))
  )

logmsg(sprintf("Writing %d sequences to %s fasta file", nrow(seqs), opt$options$fnafile))
seqs %>% arrange(seqname) %>%
  transmute(d = sprintf(">%s\n%s", seqname, seq)) %>%
  write.table(opt$options$fnafile, col.names = FALSE, row.names = FALSE, quote = FALSE)

logmsg(sprintf("Writing ASV table to %s", opt$options$outtable))
seqtab %>% inner_join(seqs, by = 'seq') %>%
  select(sample, seqname, count) %>%
  write_tsv(opt$options$outtable)

logmsg("Done")
