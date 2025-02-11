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
              help="Create and output FPKM [default %default]"),
  make_option("--tpm", action="store_true", default=TRUE,
              help="Create and output TPM [default %default]"),
  make_option("--filter_threshold", action="store", default=NULL, type = "integer",
              help="Number of correlations above 0.9 to filter based on [default %default]"),
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

packages <- c('DESeq2', 'rnaseqtools', 'dplyr', 'readr', 'tibble', 'purrr')
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
                      function(gene){ GenomicRanges::reduce(gene) |> GenomicRanges::width() |> sum() }
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

# function to calculate the correlation coefficient for a subset of genes
# compared to all genes and return the number > 0.9
calculate_cor_values_vs_all_genes <- function(
    tpm_df, tpm_no_all_zeros, samples, label, write_func) {
  a <- tpm_df |>
    dplyr::select(starts_with("zmp")) |>
    t()
  colnames(a) <- tpm_df$GeneID
  b <- tpm_no_all_zeros |>
    dplyr::select(starts_with("zmp")) |>
    t()
  colnames(b) <- tpm_no_all_zeros$GeneID
  # reorder to put the genes in "a" first
  b <- b[ , c(colnames(a), setdiff(colnames(b), colnames(a))) ]
  # calculate cor vs all other genes
  c <- cor(a, b, method = "spearman")
  # set cells where gene1 == gene2 to NA
  num_rows <- nrow(tpm_df)
  c[ seq(1,num_rows^2,num_rows + 1) ] <- NA
  # get just combinations > 0.9
  indices <- which(c > 0.9)
  col_i <- as.integer((indices - 1) / nrow(c)) + 1
  row_i <- ((indices - 1) %% nrow(c)) + 1
  # remove duplicates by insisting col_i > row_i
  to_keep <- col_i > row_i
  indices <- indices[ to_keep ]
  col_i <- col_i[ to_keep ]
  row_i <- row_i[ to_keep ]
  # output the suspect ones
  cor_by_num_zeros_out_file <- paste0(cmd_line_args$options$output_base,
                                      "-cor-by-num-zeros-", label, ".",
                                      cmd_line_args$options$output_format)
  tibble(
    Gene1 = rownames(c)[ row_i ],
    Gene2 = colnames(c)[ col_i ],
    cor = c[ indices ],
  ) |>
    arrange(cor) |>
    write_func(file = cor_by_num_zeros_out_file)
  
  # return count of correlations over 0.9
  return(length(indices))
}

# function to check where to pick the zeros threshold
filter_by_zeros <- function(tpm_df, samples, cor_count_threshold = 50, write_func) {
  t_tpm <- tpm_df |>
    dplyr::select(-c(GeneID, Name)) |>
    t()
  tpm_no_all_zeros <- tpm_df |>
    dplyr::mutate(num_zeros = apply(t_tpm, 2, \(x) (x == 0) |> sum())) |>
    dplyr::filter(num_zeros != nrow(samples))
  tpm_by_zeros <- split(tpm_no_all_zeros, tpm_no_all_zeros$num_zeros)
  num_genes <-   tpm_by_zeros |>
    purrr::imap(\(x, idx) { tibble(
      num_zeros = idx,
      num_genes = nrow(x),
      "cor_gt_0.9" = NA_integer_,
    )}) |> 
    purrr::list_rbind() |> 
    dplyr::mutate(cumul_genes = cumsum(num_genes))
  
  rev_num_zeros <- names(tpm_by_zeros) |> as.integer() |>
    sort(decreasing = TRUE) |> as.character()
  for (num_zeros in rev_num_zeros) {
    cor_gt_0.9 <- calculate_cor_values_vs_all_genes(
      tpm_by_zeros[[num_zeros]], tpm_no_all_zeros, samples, num_zeros, write_func
    )
    num_genes[ num_genes$num_zeros == num_zeros, "cor_gt_0.9"] <- cor_gt_0.9
    if (cor_gt_0.9 < cor_count_threshold) {
      break
    }
  }
  
  above_threshold <- num_genes$num_zeros[
    !is.na(num_genes$`cor_gt_0.9`) & num_genes$`cor_gt_0.9` >= cor_count_threshold
  ]
  return(dplyr::filter(tpm_no_all_zeros, !(num_zeros %in% above_threshold)))
}

# output fpkm and tpm
if (cmd_line_args$options$output_format == "tsv") {
  write_func <- readr::write_tsv
} else {
  write_func <- readr::write_csv
}
if (cmd_line_args$options$fpkm) {
  fpkm_out_file <- paste0(cmd_line_args$options$output_base, "-fpkm.",
                          cmd_line_args$options$output_format)
  write_func(rnaseq_fpkm, file = fpkm_out_file)
}
if (cmd_line_args$options$tpm) {
  tpm_out_file <- paste0(cmd_line_args$options$output_base, "-tpm.",
                         cmd_line_args$options$output_format)
  write_func(tpm, file = tpm_out_file)
}

# filter tpms if filter_threshold set
if (!is.null(cmd_line_args$options$filter_threshold)) {
  cor_by_num_zeros_out_file <- paste0(cmd_line_args$options$output_base, 
                                  "-cor-by-num-zeros.",
                                  cmd_line_args$options$output_format)
  tpm_filtered <- filter_by_zeros(
      tpm,
      samples, 
      cor_count_threshold = cmd_line_args$options$filter_threshold,
      write_func
    )
  tpm_filtered_out_file <- paste0(cmd_line_args$options$output_base, 
                                  "-tpm-filtered-by-zeros.",
                                  cmd_line_args$options$output_format)
  
  write_func(tpm_filtered, file = tpm_filtered_out_file)
}

# output counts as well
counts_out_file <- paste0(cmd_line_args$options$output_base, "-counts.",
                          cmd_line_args$options$output_format)
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
