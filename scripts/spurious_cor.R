#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to investigate count data", sep = "\n")

if (interactive()) {
  cmd_args <- c("zmp_ph204", "zmp_ph192") # add test files here
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc,
    usage = "Usage: %prog [options] input_file"
  ),
  args = cmd_args,
  positional_arguments = c(1, Inf)
)

packages <- c('tidyverse', 'rnaseqtools')
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

# Load data
# samples <- rnaseqtools::load_rnaseq_samples(cmd_line_args$args[1])
# count_data <- rnaseqtools::load_rnaseq_data(cmd_line_args$args[2])
# counts <- rnaseqtools::get_counts(count_data)
# norm_data <- rnaseqtools::normalise_counts(count_data, samples)
# 
# norm_counts <- rnaseqtools::get_counts(norm_data, normalised = TRUE)
# t_norm <- norm_counts |> t()
# all_zero <- apply(t_norm, 2, sum) == 0
# t_norm <- t_norm[,!all_zero]
# 
# norm_data <- norm_data |> 
#   filter(!all_zero) |> 
#   mutate(
#     t_norm_var = apply(t_norm, 2, var),
#     var_quantile = cut(
#       t_norm_var, 
#       c(0, quantile(t_norm_var, seq(0.1,1,0.1))),
#       labels = paste(seq(0,90,10), paste0(seq(10,100,10), "%"), sep = "-")
#     ),
#     num_zeros = apply(t_norm, 2, \(x) (x == 0) |> sum()),
#     mean_norm_counts = apply(t_norm, 2, mean),
#     expr_quantile = cut(
#       mean_norm_counts,
#       c(0, quantile(mean_norm_counts, seq(0.1,1,0.1))),
#       labels = paste(seq(0,90,10), paste0(seq(10,100,10), "%"), sep = "-")
#     )
#   ) |> 
#   relocate(t_norm_var:expr_quantile, .after = Name)
# 
# # Look at num zeros vs variance by percentile
# norm_data |> ggplot(aes(x = num_zeros)) +
#   geom_histogram(binwidth = 1) +
#   facet_wrap(vars(var_quantile), ncol = 5) +
#   labs(
#     title = "Genes with the lowest variances have the most numbers of zero values",
#     subtitle = "Histogram of zero values facetted by variance percentile"
#   )
# # Look at num zeros vs expression percentile
# norm_data |> ggplot(aes(x = num_zeros)) +
#   geom_histogram(binwidth = 1) +
#   facet_wrap(vars(expr_quantile), ncol = 5) +
#   labs(
#     title = "Genes with the lowest mean expression levels have the most numbers of zero values",
#     subtitle = "Histogram of zero values facetted by mean expression percentile"
#   )
# 
# # Look at num_zeros vs mean counts
# cutoff <- filter(norm_data, num_zeros >= 12) |> 
#   pull(mean_norm_counts) |> 
#   max()
# norm_data |> ggplot(aes(x = mean_norm_counts, y = num_zeros)) +
#   geom_vline(xintercept = cutoff, colour = "firebrick3") +
#   geom_point() +
#   lims(x = c(0,20))

# load TPMs
load_tpms <- function(expt) {
  expt_dir <- file.path("nf", "results", expt)
  tpms_file <- file.path(expt_dir, paste(expt, "all-tpm.tsv", sep = "-"))
  read_tsv(tpms_file)
}

load_samples <- function(expt) {
  expt_dir <- file.path("nf", "results", expt)
  sample_file <- file.path(expt_dir, "samples.tsv")
  read_tsv(sample_file)
}

filter_and_transpose_tpms <- function(tpms, expt) {
  t_tpm <- tpms |> 
    dplyr::select(-c(GeneID, Name)) |> 
    t()
  # remove where all values are zero
  all_zero <- apply(t_tpm, 2, sum) == 0
  cat(glue::glue("Expt {expt}: {sum(all_zero)} genes are zero across all samples"), "\n")
  t_tpm <- t_tpm[,!all_zero]
  filtered_tpms <-   tpms |> 
    filter(!all_zero)
  return(list(filtered_tpms, t_tpm))
}

add_metadata <- function(tpms, t_tpm) {
  tpms |> mutate(
      tpms_var = apply(t_tpm, 2, var),
      var_quantile = cut(
        tpms_var, 
        c(0, quantile(tpms_var, seq(0.1,1,0.1))),
        labels = paste(seq(0,90,10), paste0(seq(10,100,10), "%"), sep = "-")
      ),
      num_zeros = apply(t_tpm, 2, \(x) (x == 0) |> sum()),
      mean_tpm = apply(t_tpm, 2, mean),
      expr_quantile = cut(
        mean_tpm,
        c(0, quantile(mean_tpm, seq(0.1,1,0.1))),
        labels = paste(seq(0,90,10), paste0(seq(10,100,10), "%"), sep = "-")
      )
    ) |> 
    relocate(tpms_var:expr_quantile, .after = Name)
}

num_zero_plots <- function(tpms) {
  plot_list <- list(
    "num_zeros_by_var_q" = NULL,
    "num_zeros_by_expr_q" = NULL,
    "num_zeros_by_mean_expr" = NULL,
  )
  # Look at num zeros vs variance by percentile
  plot_list[["num_zeros_by_var_q"]] <- tpms |> ggplot(aes(x = num_zeros)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(vars(var_quantile), ncol = 5) +
    labs(
      title = "Genes with the lowest variances have the most numbers of zero values",
      subtitle = "Histogram of zero TPM values facetted by variance percentile"
    )
  # Look at num zeros vs expression percentile
  plot_list[["num_zeros_by_expr_q"]] <- tpms |> ggplot(aes(x = num_zeros)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(vars(expr_quantile), ncol = 5) +
    labs(
      title = "Genes with the lowest mean expression levels have the most numbers of zero values",
      subtitle = "Histogram of zero TPM values facetted by mean expression percentile"
    )
  # Look at num_zeros vs mean tpm
  cutoff <- filter(tpms, num_zeros >= 12) |> 
    pull(mean_tpm) |> 
    max()
  plot_list[["num_zeros_by_mean_expr"]] <- tpms |> ggplot(aes(x = mean_tpm, y = num_zeros)) +
    geom_vline(xintercept = cutoff, colour = "firebrick3") +
    geom_point() +
    lims(x = c(0,20))
  
  return(plot_list)
}

plot_cor_hist <- function(tpms) {
  cor_mat <- tpms |>
    dplyr::select(starts_with("zmp")) |>
    t() |> cor(method = "spearman")
  rownames(cor_mat) <- tpms$GeneID
  colnames(cor_mat) <- tpms$GeneID
  cor_long <- as_tibble(cor_mat, rownames = "Gene1") |>
    pivot_longer(cols = -Gene1, names_to = "Gene2", values_to = "cor")
  cor_long |> filter(Gene1 != Gene2) |>
    ggplot(aes(x = cor)) +
    geom_histogram(binwidth = 0.1, boundary = 0) +
    scale_y_log10() +
    lims(x = c(-1,1))
}

effect_of_zeros <- function(expt) {
  tpms <- load_tpms(expt)
  samples <- load_samples(expt)
  tpm_list <- filter_and_transpose_tpms(tpms, expt)
  tpms_no_all_zeros <- add_metadata(tpm_list[[1]], tpm_list[[2]])
  # filter by num_zeros
  num_zeros_cutoff <- ceiling(nrow(samples) * 0.25)
  # no more than 1/4 of samples an be zero
  filtered_by_zeros <- tpms_no_all_zeros |> 
    filter(num_zeros < num_zeros_cutoff)
  cat(glue::glue("Expt {expt}: {nrow(filtered_by_zeros)} genes have zeros in fewer then 1/4 of samples"), "\n")
  
  plots <- num_zero_plots(tpms_no_all_zeros)
  expt_dir <- file.path("nf", "results", expt)
  pdf(file.path(expt_dir, "num_zero_plots.pdf"))
  walk(plots, print)
  mean_tpm_dist <- filtered_by_zeros |> 
    ggplot(aes(x = mean_tpm)) +
    geom_histogram() +
    scale_x_log10() +
    labs(title = "Mean TPM distribution after filtering by zeros")
  print(mean_tpm_dist)
  # sample genes and plot hist of cor values
  sub_sample <- filtered_by_zeros |> slice_sample(prop = 0.2)
  cor_hist <- plot_cor_hist(sub_sample) +
    labs(title = "Genes filtered by zero don't have an excess of extreme correlation values")
  print(cor_hist)
  # sample genes with lots of zeros and plot hist of cor values
  sub_sample <- tpms_no_all_zeros |>
    filter(num_zeros > num_zeros_cutoff) |>
    slice_sample(prop = 0.4)
  cor_hist <- plot_cor_hist(sub_sample) +
    labs(title = "Genes with more zeros are skewed towards high correlation values")
  print(cor_hist)
  dev.off()
}

walk(cmd_line_args$args, effect_of_zeros)

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
