#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--transcripts_file", action="store", default=NULL,
              help="File of transcript lengths [default %default]"),
  make_option("--transcripts_db", action="store", default=NULL,
              help="TxDb file [default %default]"),
  make_option("--gtf_file", action="store", default="../reference/grcz11/Danio_rerio.GRCz11.111.gtf",
              help="GTF annotation file [default %default]"),
  make_option("--fpkm", action="store_true", default=FALSE,
              help="Create and output fpkm [default %default]"),
  make_option("--tpm", action="store_true", default=TRUE,
              help="Create and output fpkm [default %default]"),
  make_option("--output_base", action="store", default="all",
              help="Base filename of output files [default %default]"),
  make_option("--output_format", action="store", default="tsv",
              help="Format (tsv/csv) to output [default %default]"),
  make_option(c("-d", "--debug"), action="store_true", default=FALSE,
              help="Print extra output [default %default]")
)

desc <- paste('Script to convert DESeq2 counts to FPKM and/or TPM', 
              'Expects an all.tsv file from the DESeq R script', sep = "\n")
if (interactive()) {
  cmd_args <- c('samples.tsv', 'all.tsv')
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list,
    description = desc,
    usage = "Usage: %prog [options] samples_file counts_file" ),
  args = cmd_args,
  positional_arguments = 2
)

packages <- c('DESeq2', 'rnaseqtools', 'dplyr', 'readr', 'tibble')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# check format option
if (!(cmd_line_args$options[['output_format']] %in% c('tsv', 'csv'))) {
  stop("--output_format option must be one of 'tsv' or 'csv'")
}

samples_file <- cmd_line_args$args[1]
samples <- readr::read_tsv(samples_file, col_names = c("sample", "condition"),
                           show_col_types = FALSE)

in_file <- cmd_line_args$args[2]
rnaseq_data <- load_rnaseq_data(in_file)

rnaseq_counts <- get_counts(rnaseq_data, samples)

# make a DESeq2 object
design <- formula(~ condition)
dds <- DESeqDataSetFromMatrix(
  rnaseq_counts,
  samples,
  design=design
)
rowData(dds) <- get_gene_metadata(rnaseq_data)

# check for the existence of Transcript info
if (!is.null(cmd_line_args$options$transcripts_file) &
    file.exists(cmd_line_args$options$transcripts_file)){
  # load transcript info
  TxInfo <- read.table(file=cmd_line_args$options$transcripts_file, sep="\t", row.names = 1)
  TxLengths <- TxInfo[[1]]
  names(TxLengths) <- rownames(TxInfo)
} else{
  # load Genomic Features. Not required if there is a pre-existing transcripts file
  library('GenomicFeatures')
  # check for the existence of TxDb file
  if (!is.null(cmd_line_args$options$transcripts_db) & 
      file.exists(cmd_line_args$options$transcripts_db)) {
    txdb <- loadDb(cmd_line_args$options$transcripts_db)
  } else {
    # make TxDb from GFF file
    txdb <- makeTxDbFromGFF(
      cmd_line_args$options$gtf_file, 
      format="gtf",
      dataSource = "Ensembl",
      organism = "Danio rerio",
      taxonomyId = 7955
    )
    # save to db file
    if (!is.null(cmd_line_args$options$transcripts_db)) {
      saveDb(txdb, cmd_line_args$options$transcripts_db)
    }
  }
  # calculate lengths of transcripts
  exonsByGene <- exonsBy(txdb, by="gene")
  TxLengths <- sapply(exonsByGene, 
                      function(gene){ reduce(gene) |> width() |> sum() }
  )
  if (!is.null(cmd_line_args$options$transcripts_file)) {
    write.table(TxLengths, file=cmd_line_args$options$transcripts_file, 
                row.names=TRUE, sep="\t", quote=FALSE, col.names = FALSE)
  }
}

# remove genes from deseq object that are not in TxLengths
dds <- dds[ rowData(dds)$GeneID %in% names(TxLengths), ]
# order tx lengths by deseq object
TxLengths <- TxLengths[ rowData(dds)$GeneID ]

# add tx lengths to deseq object
assays(dds, withDimnames = FALSE)[["avgTxLength"]] <- 
  matrix(rep(TxLengths, ncol(dds)), ncol = ncol(dds))

rnaseq_fpkm <- fpkm(dds)
rownames(rnaseq_fpkm) <- rowData(dds)$GeneID

# calculate tpm from fpkm
fpkm_to_tpm <- function(fpkm) {
  total_fpkm_by_sample <- colSums(fpkm)
  do.call(cbind, lapply(1:ncol(fpkm), function(i) {
    exp( log(fpkm[,i]) - log(total_fpkm_by_sample[i]) + log(1e6) )
  }))
}
tpm <- fpkm_to_tpm(rnaseq_fpkm)
colnames(tpm) <- colnames(rnaseq_fpkm)

# join to metadata
pval_cols <- c("pval", "adjp", "log2fc")
metadata <- rowData(dds) |> as_tibble() |> 
  dplyr::select(-any_of(pval_cols))

rnaseq_fpkm <- inner_join(
  metadata,
  as_tibble(rnaseq_fpkm, rownames = "GeneID"),
  by = join_by(GeneID == GeneID)
)

tpm <- inner_join(
  metadata,
  as_tibble(tpm, rownames = "GeneID"),
  by = join_by(GeneID == GeneID)
)

# output fpkm and tpm
if (cmd_line_args$options$output_format == "tsv") {
  write_func <- readr::write_tsv
  fpkm_out_file <- paste0(cmd_line_args$options$output_base, '-fpkm.tsv')
  tpm_out_file <- paste0(cmd_line_args$options$output_base, '-tpm.tsv')
  counts_out_file <- paste0(cmd_line_args$options$output_base, '-counts.tsv')
} else {
  write_func <- readr::write_csv
  fpkm_out_file <- paste0(cmd_line_args$options$output_base, '-fpkm.csv')
  tpm_out_file <- paste0(cmd_line_args$options$output_base, '-tpm.csv')
  counts_out_file <- paste0(cmd_line_args$options$output_base, '-counts.csv')
}
if (cmd_line_args$options$fpkm) {
  write_func(rnaseq_fpkm, file = fpkm_out_file)
}
if (cmd_line_args$options$tpm) {
  write_func(tpm, file = tpm_out_file)
}

# output counts as well
rnaseq_data |> dplyr::select(-any_of(pval_cols)) |> 
  write_func(file = counts_out_file)

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2024 Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
