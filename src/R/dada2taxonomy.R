#!/usr/bin/env Rscript

# dada2taxonomy.R
#
# Determines taxonomy for amplicon sequences using methods implemented in the DADA2 package.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

SCRIPT_VERSION = "0.1"

# opt variable for testing: opt <- list(options = list(idtaxa_rdata = '../../testdata/SILVA_SSU_r132_March2018.RData', rdp_fasta = '../../testdata/silva.fna.gz', species_fasta = '../../testdata/silva_species.fna.gz', verbose = TRUE), args = c('../../testdata/fifty_16S_sequences.fna'))
# Get arguments
option_list = list(
  make_option(
    c("--idtaxa_rdata"), action="store", default = '',
    help = "Name of RData file with classified sequences for the IDTAXA method."
  ),
  make_option(
    c("--rdp_fasta"), action="store", default = '',
    help = "Name of fasta file with classified sequences for the naive bayesian method (RDP)."
  ),
  make_option(
    c("--species_fasta"), action="store", default = '',
    help = "Name of fasta file with classified sequences for species assignment."
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
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] sequences0.fna ... sequencesn.fna
    
	Taxonomy for sequences in sequences0.fna to sequencesn.fna will be determined using
        methods available, see option list.
        
        Output files are tab separated and uses the name of the input file as basename
        followed by the method name.
        
        Data for the RDP and the species classifier can be downloaded from 
        https://benjjneb.github.io/dada2/training.html and for the IDTAXA classifier 
        from http://www2.decipher.codes/Downloads.html.", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
  quit('no')
}

# This is here because it takes a long time to load, so we want to avoid
# loading it if the user uses the --help or --version flags.
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))

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

for ( f in opt$args ) {
  bf = sub('(.+)(\\.fasta)|(\\.fna)|(\\.fa)|(\\.RData)(.gz)?', '\\1', basename(f))
  logmsg(sprintf("Handling %s, writing output to files starting with %s", f, bf))
  seqs <- readDNAStringSet(f)
  seqs_df <- seqs %>% data.frame() %>% tibble::rownames_to_column('s')
  colnames(seqs_df) <- c('seqname', 'sequence')
  if ( stringr::str_length(opt$options$rdp_fasta) > 0 ) {
    db = case_when(
      grepl('gtdb',  opt$options$rdp_fasta, ignore.case = TRUE) ~ 'gtdb.',
      grepl('silva', opt$options$rdp_fasta, ignore.case = TRUE) ~ 'silva.',
      grepl('unite', opt$options$rdp_fasta, ignore.case = TRUE) ~ 'unite.',
      TRUE ~ ''
    )
    logmsg(sprintf("\tAssigning with RDP classifier using %s", opt$options$rdp_fasta, multithread = TRUE))
    rdp <- assignTaxonomy(seqs, opt$options$rdp_fasta, multithread = TRUE)

    if ( stringr::str_length(opt$options$species_fasta) > 0 ) {
      logmsg(sprintf("\tAssigning species using %s", opt$options$species_fasta))
      rdp <- addSpecies(rdp, opt$options$species_fasta)
    }

    colnames(rdp) <- colnames(rdp) %>% tolower()

    of <- sprintf("%s.rdp.%stsv.gz", bf, db)
    logmsg(sprintf("\tWriting table to %s", of))
    rdp %>% data.frame() %>% tibble::rownames_to_column('sequence') %>%
      inner_join(seqs_df, by = 'sequence') %>%
      select(seqname, 2:ncol(.), sequence) %>%
      write_tsv(of)
  }

  if ( stringr::str_length(opt$options$idtaxa_rdata) > 0 ) {
    db = case_when(
      grepl('gtdb',  opt$options$idtaxa_rdata, ignore.case = TRUE) ~ 'gtdb.',
      grepl('silva', opt$options$idtaxa_rdata, ignore.case = TRUE) ~ 'silva.',
      grepl('unite', opt$options$idtaxa_rdata, ignore.case = TRUE) ~ 'unite.',
      TRUE ~ ''
    )
    logmsg(sprintf("\tAssigning with IDTAXA using %s", opt$options$idtaxa_rdata))
    load(opt$options$idtaxa_rdata)
    ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
    idt <- IdTaxa(seqs, trainingSet, processors = NULL) %>%
      sapply(
        function(x) { 
          m <- match(ranks, x$rank)
          taxa <- x$taxon[m]
          taxa[startsWith(taxa, "unclassified_")] <- NA
          taxa 
        }
      ) %>% t() %>%
      data.frame()
    colnames(idt) <- ranks
    rownames(idt) <- seqs

    if ( stringr::str_length(opt$options$species_fasta) > 0 ) {
      logmsg(sprintf("\tAssigning species using %s", opt$options$species_fasta))
      idt <- addSpecies(idt, opt$options$species_fasta)
    }

    of <- sprintf("%s.idtaxa.%stsv.gz", bf, db)
    logmsg(sprintf("\tWriting table to %s", of))
    idt <- idt %>% data.frame() %>% tibble::rownames_to_column('sequence') 
    if ( stringr::str_length(opt$options$species_fasta) > 0 ) idt <- rename(idt, dada2_species = Species)
    idt %>%
      inner_join(seqs_df, by = 'sequence') %>%
      select(seqname, 2:ncol(.), sequence) %>%
      write_tsv(of)
  }
}

logmsg("Done")
