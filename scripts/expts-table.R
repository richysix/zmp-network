#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to produce tables of the expt/sample info", sep = "\n")

if (interactive()) {
  cmd_args <- c(
    samples = file.path("data", "expt-sample-condition-initial.tsv"),
    expts = file.path("data", "expt-info.txt")) # add test files here
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc,
    usage = "Usage: %prog [options] samples-file expt-info-file"
  ),
  args = cmd_args,
  positional_arguments = 2
)

packages <- c('tidyverse', 'gt')
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

# load all samples and count samples per expt
bin_breaks <- c(0, 12, 24, 48, 96)
bin_labels <- paste(
  bin_breaks[seq_len(length(bin_breaks) - 1)] + 1,
  bin_breaks[seq_len(length(bin_breaks) - 1) + 1],
  sep = "-"
)
sample_counts <- read_tsv(cmd_line_args$args[1], show_col_types = FALSE) |> 
  count(expt) |> 
  mutate(
    expt = fct_reorder(expt, n),
    n_bin = cut(n, breaks = bin_breaks, labels = bin_labels)
  ) |> 
  arrange(n)

expt_info <- read_tsv(cmd_line_args$args[2], show_col_types = FALSE) |> 
  dplyr::select(expt, gene_name, allele_id)

samples_table <- sample_counts |> 
  inner_join(expt_info, by = join_by(expt == expt)) |> 
  group_by(expt, n) |> 
  summarise(
    "Gene(s)" = paste(gene_name, collapse = "/"),
    "Allele(s)" = paste(allele_id, collapse = "/"),
    .groups = "drop"
  ) |> 
  rename(Expt = expt) |> 
  arrange(n) |> 
  gt() |> 
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(`Gene(s)`, `Allele(s)`),
    )
  )

samples_table |> 
  gtsave(filename = file.path("plots", "sample-n-table.pdf"))

samples_table |> 
  gtsave(filename = file.path("plots", "sample-n-table.png"), 
         vwidth = 1200,
         zoom = 2)

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
