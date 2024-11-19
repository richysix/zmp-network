#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to create some test data for the nextflow pipeline", 
              sep = "\n")

if (interactive()) {
  cmd_args <- c() # add test files here
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc,
    usage = "Usage: %prog [options]"
  ),
  args = cmd_args,
  positional_arguments = 0
)

packages <- c('tidyverse', 'faux', 'rprojroot')
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

root_path <- rprojroot::find_root(rprojroot::is_rstudio_project)
data_dir <- file.path(root_path, "data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

cors <- matrix(c(1, 0.8, 0.75, 0.78, 
                 0.8, 1, 0.78, 0.77, 
                 0.75, 0.78, 1, 0.8, 
                 0.78, 0.77, 0.8, 1), ncol = 4)

samples_per_expt <- 12
expts <- c("test-1", "test-2")
samples <- tibble(
  expt = rep(expts, each = samples_per_expt),
  sample = paste0('sample-', seq_len(length(expts) * samples_per_expt)),
  condition = rep(rep(c("sib", "mut"), each = samples_per_expt/2), length(expts))
) 
samples |> write_tsv(file = file.path(data_dir, 'test-samples.tsv'))

set.seed(665)
expr <- rnorm_multi(n = samples_per_expt * length(expts), vars = ncol(cors), 
                    mu = c(5, 10, 15, 20), r = cors) |> 
  t() |> 
  magrittr::set_colnames(paste0(samples$sample, " normalised count")) |> 
  magrittr::set_rownames(paste0('ENSDARG0000000000', seq_len(ncol(cors)))) |> 
  as_tibble(expr, rownames = "GeneID") |> 
  mutate(across(starts_with("sample"), trunc, .names = "{.col} int"),
         `Gene name` = paste0('gene-', seq_len(ncol(cors))) ) |> 
  rename_with(\(x) sub("normalised count int", "count", x)) |> 
  relocate(matches("normalised"), .after = `Gene name`) |> 
  relocate(`Gene name`, .after = GeneID)

write_csv(expr, file = file.path(data_dir, "test-counts.csv"))
write_csv(expr, file = file.path(data_dir, "test-counts.csv.gz"))

# create fake transcripts to make fpkm normalisation produce sensible numbers
tibble(
  "GeneID" = expr$GeneID,
  "Length" = 1e9
) |> write_tsv(file = file.path(data_dir, "test-transcript-lengths.tsv"),
               col_names = FALSE)

# create fake annotation file
set.seed(67)
tibble(
  "GeneID" = expr$GeneID,
  "Chr" = sample(1:25, ncol(cors)),
  "Start" = seq(100, 400, 100),
  "End" = c(200, 500, 600, 350),
  "Strand" = rep("1", ncol(cors)),
  "Biotype" = rep("protein_coding", ncol(cors)),
  "Name" = paste("gene", seq_len(ncol(cors)), sep = "-"),
  "Description" = paste("gene", seq_len(ncol(cors)), sep = "-")
) |> write_tsv(file = file.path(data_dir, "test-annotation.txt"),
               col_names = FALSE)

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
