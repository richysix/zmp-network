#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("", sep = "\n")

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

set.seed(665)
expr <- rnorm_multi(n = 12, vars = 4, mu = c(5, 10, 15, 20), r = cors) |> 
  t() |> 
  magrittr::set_colnames(paste0('sample-', 1:12, " normalised count")) |> 
  magrittr::set_rownames(paste0('ENSDARG0000000000', 1:4)) |> 
  as_tibble(expr, rownames = "GeneID") |> 
  mutate(across(starts_with("sample"), trunc, .names = "{.col} int"),
         `Gene name` = paste0('gene-', 1:4)) |> 
  rename_with(\(x) sub("normalised count int", "count", x)) |> 
  relocate(matches("normalised"), .after = `Gene name`) |> 
  relocate(`Gene name`, .after = GeneID)

write_csv(expr, file = file.path(data_dir, "test-counts.csv"))
write_csv(expr, file = file.path(data_dir, "test-counts.csv.gz"))

tibble(
  expt = rep("test", 12),
  sample = paste0('sample-', 1:12),
  condition = rep(c("sib", "mut"), each = 6)
) |> 
  write_tsv(file = file.path(data_dir, 'test-samples.tsv'))

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
