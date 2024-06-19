#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='all-agg.csv',
              help="Output file name [default %default]" ),
  make_option("--debug", type="logical", default=FALSE, action="store_true",
              help="Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to aggregate the counts from ends to genes for a DeTCT experiment',
  sep = "\n"
)
if (interactive()) {
  cmd_args <- c('/data/scratch/bty114/detct/grcz11/everything/filter-strict/all-test.csv.gz', 
                '/data/scratch/bty114/detct/grcz11/everything/filter-strict/sample-subset-dnmt8.tsv')
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'aggregate_counts_to_genes_detct.R',
    usage = "Usage: %prog [options] all_file samples_file",
    description = desc),
  args = cmd_args,
  positional_arguments = 2
)

# load packages
packages <- c('tidyverse', 'magrittr', 'rnaseqtools', 'DESeq2')
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# unpack options
output_file <- cmd_line_args$options[['output_file']]

Sys.setenv("VROOM_CONNECTION_SIZE" = 2^20)

# load samples file
# do this before loading the detct data in case there's a problem
# saves loading all the detct data first
samples <- load_rnaseq_samples(cmd_line_args$args[2])

# load RNAseq data (all file)
data <- load_detct_data(data_file = cmd_line_args$args[1])

# get raw counts and add gene id column
# remove rows not associated with a gene
# remove rows associated with more than one gene
counts <- get_counts(data, samples) %>% 
  mutate(., GeneID = data$GeneID) %>% 
  filter(., GeneID != "-" & !grepl(",", GeneID)) %>% 
  relocate(., GeneID) %>% 
  mutate(., GeneID = factor(GeneID))

# aggregate to gene level
counts_agg <- counts |> 
  group_by(GeneID) |> 
  summarise(across(.cols = all_of(setdiff(colnames(counts), c("GeneID"))), sum))

# make DESeq2 object and calculate normalised counts
dds <- DESeqDataSetFromMatrix(
  countData = column_to_rownames(counts_agg, var="GeneID"),
  colData = column_to_rownames(samples, var="sample"),
  design = formula(~ 1)
)
dds <- estimateSizeFactors(dds)

# get normalised counts
counts <- as.data.frame(counts(dds, normalized=FALSE)) %>% 
  set_colnames(., sub("$", " count", colnames(.))) %>% 
  rownames_to_column(., var = "GeneID")

norm_counts <- as.data.frame(counts(dds, normalized=TRUE)) %>% 
  set_colnames(., sub("$", " normalised count", colnames(.))) %>% 
  rownames_to_column(., var = "GeneID")

# out gene info, counts and normalised counts
if (grepl("\\.gz$", output_file)) {
  fh <- gzfile(output_file, open = "wb")
} else {
  fh <- open(output_file, open = "w")
}
if (grepl("\\.csv", output_file)) {
  write_func <- readr::write_csv
} else {
  write_func <- readr::write_tsv
}
select(data, GeneID, Name, Description) %>% 
  filter(., GeneID != "-") %>% 
  unique() %>% 
  inner_join(., counts, by = join_by(GeneID)) %>% 
  inner_join(., norm_counts, by = join_by(GeneID)) %>% 
  write_func(., file = fh)
close(fh)
