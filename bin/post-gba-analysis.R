#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(
    "--samples_file", type = "character",
    default = "expt-sample-condition-test.tsv",
    help = "Name of the overall samples file [default %default]"
  ),
  make_option(
    "--min_cluster_size", type = "integer", default = 100,
    help = paste0("Lower bound of the range to calculate the max ecdf over ",
                  "[default %default]")
  ),
  make_option(
    "--max_cluster_size", type = "integer", default = 1000,
    help = paste0("Upper bound of the range to calculate the max ecdf over ",
                  "[default %default]")
  ),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to create overview plots for GBA analysis", sep = "\n")

cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc,
    usage = "Usage: %prog [options]"),
  positional_arguments = 0
)

packages <- c("tidyverse", "gt")
for (package in packages) {
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# helper function for converting inflation numbers for axes
divide_by_10 <- function(x) { (as.integer(x)/10) |> as.character() }

# find all expt names
expts <- read_tsv("expts.txt", show_col_types = FALSE, col_names = c("expt"))
# load all samples and count samples per expt
bin_breaks <- c(0, 12, 24, 48, 96)
bin_labels <- paste(
  bin_breaks[seq_len(length(bin_breaks) - 1)] + 1,
  bin_breaks[seq_len(length(bin_breaks) - 1) + 1],
  sep = "-"
)
sample_counts <- read_tsv(cmd_line_args$options$samples_file, show_col_types = FALSE) |>
  filter(expt %in% expts$expt) |>
  count(expt) |>
  mutate(
    expt = fct_reorder(expt, n),
    sample_n_bin = cut(n, breaks = bin_breaks, labels = bin_labels)
  )
colour_palette <- biovisr::cbf_palette(sample_counts$expt,
                                       named = TRUE)

## AUC VALUES
# function to summarise AUC values
summarise_auc <- function(filename) {
  expt <- str_remove(filename, "-all-tpm.*$")
  file_info <- basename(filename) |>
    str_remove("^.*all-tpm-") |>
    str_split_1("\\.")
  read_tsv(filename, show_col_types = FALSE) |>
    summarise(
      mean_auc = mean(auc),
      mean_null_auc = mean(degree_null_auc)
    ) |> mutate(
      diff = mean_auc - mean_null_auc,
      expt = expt,
      params = file_info[1],
      inflation = file_info[3] |> str_remove("I"),
      domain = file_info[4]
    )
}

auc_files <- list.files(pattern = "*\\.auc\\.tsv", recursive = TRUE)

auc_data <- auc_files |>
  map(summarise_auc) |>
  list_rbind() |> 
  mutate(
    threshold_method = case_when(
      grepl("-k", params) ~ "knn",
      TRUE ~ "cor"
    ), 
    threshold = case_when(
      threshold_method == "knn" ~ str_remove(params, "^t[0-9]*-k"),
      threshold_method == "cor" ~ str_remove(params, "^t"),
    ),
    threshold = as.integer(threshold),
    expt = factor(expt, levels = levels(sample_counts$expt))
  ) |> 
  arrange(desc(threshold_method), threshold) |> 
  mutate(params = factor(params, levels = unique(params)))

dodge <- position_dodge(width = 0.75)
auc_line_points <- auc_data |> 
  pivot_longer(cols = c(mean_auc, mean_null_auc)) |> 
  ggplot(mapping = aes(
    x = params,
    y = value,
    colour = inflation,
    shape = name,
    group = interaction(params, inflation)
  )) +
  geom_point(position = dodge) +
  geom_line(position = dodge, show.legend = FALSE) +
  facet_grid(rows = vars(domain), cols = vars(expt),
             labeller = labeller(domain = \(x) toupper(x))) +
  scale_colour_manual(values = biovisr::cbf_palette(factor(auc_data$inflation)),
                  labels = divide_by_10) +
  scale_shape_discrete(name = "Value", labels = c("Mean AUC", "Mean null AUC"))
  
auc_diff <- ggplot(data = auc_data) +
  geom_col(aes(x = params, y = diff, fill = inflation),
           position = position_dodge()) +
  facet_grid(rows = vars(domain), cols = vars(expt),
             labeller = labeller(domain = \(x) toupper(x))) +
  scale_fill_manual(values = biovisr::cbf_palette(factor(auc_data$inflation)),
                    labels = divide_by_10) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

miscr::open_graphics_device(
  filename = file.path(plot_dir, "auc-diff.pdf"),
  width = 12,
  height = 4.8
)
print(auc_line_points)
print(auc_diff)
dev.off()

# CLUSTER SIZES
# function to load cluster size data
load_cluster_sizes <- function(filename) {
  # parse filename for info
  expt <- str_remove(filename, "-all-tpm.*$")
  file_info <- basename(filename) |> 
    str_remove("^.*all-tpm-") |> 
    str_split_1("\\.")
  # read in data and add metadata
  read_tsv(filename, show_col_types = FALSE) |> 
    mutate(
      expt = expt,
      params = file_info[1],
      inflation = file_info[3] |> str_remove("I")
    )
}

# find all cluster size files
cluster_sizes_files <- list.files(pattern = "*\\.cl-sizes\\.tsv")

# read in data
cl_size_data <- map(cluster_sizes_files, load_cluster_sizes) |>
  list_rbind() |>
  mutate(
    threshold_method = case_when(
      grepl("-k", params) ~ "knn",
      TRUE ~ "cor"
    ), 
    threshold = case_when(
      threshold_method == "knn" ~ str_remove(params, "^t[0-9]*-k"),
      threshold_method == "cor" ~ str_remove(params, "^t"),
    ),
    threshold = as.integer(threshold),
    threshold = case_when(
      threshold_method == "cor" ~ threshold / 100,
      threshold_method == "knn" ~ threshold,
    ),
    expt = factor(expt, levels = levels(sample_counts$expt))
  ) |> 
  arrange(desc(threshold_method), threshold) |> 
  mutate(params = factor(params, levels = unique(params)))

cl_size_data_by_method <- split(cl_size_data, cl_size_data$threshold_method)
cl_size_summary <- function(cl_size_df){
  cl_size_df |>
  mutate(threshold = factor(
    threshold, levels = unique(threshold) |> as.numeric() |> sort())
  ) |>
  group_by(expt, params, inflation, threshold) |>
  summarise(
    num_nodes = sum(cluster_size),
    total_clusters = n(),
    singleton_nodes = sum(cluster_size == 1),
    non_singleton_nodes = sum(cluster_size[ cluster_size > 1]),
    .groups = "drop") |> 
  mutate(non_singleton_clusters = total_clusters - singleton_nodes) 
}
cl_size_summary_by_method <- map(cl_size_data_by_method, cl_size_summary)

num_clusters_plot <- function(cl_size_summary_df) {
  cl_size_summary_df |> 
    ggplot(aes(x = threshold, y = non_singleton_clusters, fill = factor(inflation))) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = biovisr::cbf_palette(factor(cl_size_summary_df$inflation)),
                      name = "Inflation", 
                      labels = divide_by_10 ) +
    facet_wrap(vars(expt)) +
    labs(y = "Number of clusters (excluding singletons)") +
    theme_minimal()
}

num_non_singleton_nodes_plot <- function(cl_size_summary_df) {
  cl_size_summary_df |> 
    ggplot(aes(x = threshold, y = non_singleton_nodes, fill = factor(inflation))) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = biovisr::cbf_palette(factor(cl_size_summary_df$inflation)),
                      name = "Inflation", 
                      labels = divide_by_10 ) +
    facet_wrap(vars(expt)) +
    labs(y = "Number of nodes in non-singleton clusters") +
    theme_minimal()
}

cluster_size_histogram <- function(dat) {
  total_nodes <- dat |> group_by(threshold, inflation) |>
    summarise(sum = sum(cluster_size), .groups = "drop") |> 
    pull(sum) |> max()
  ten_pc_nodes <- ceiling(total_nodes/10)
  dat |> mutate(threshold = factor(
    threshold, levels = unique(threshold) |> sort())
  ) |> 
    ggplot(aes(x = cluster_size)) +
    geom_histogram(binwidth = 1, boundary = 1) +
    geom_hline(yintercept = ten_pc_nodes, colour = "firebrick3") +
    facet_grid(rows = vars(inflation), cols = vars(threshold),
               labeller = labeller(inflation = divide_by_10)) +
    scale_x_log10(limits = c(0.1,NA)) +
    scale_y_log10() +
    theme_minimal()
}

facetted_ecdf <- function(dat) {
  dat |>
    arrange(cluster_size) |>
    group_by(expt, threshold, inflation, cluster_size) |>
    summarise(
      total_nodes = sum(cluster_size), 
      .groups = "drop_last"
    ) |>
    mutate(
      frac = total_nodes / sum(total_nodes),
      ecdf = cumsum(frac)) |>
    ggplot(aes(x = cluster_size, y = ecdf, colour = inflation)) +
    geom_step() +
    facet_grid(rows = vars(expt), cols = vars(threshold)) +
    scale_x_log10() +
    scale_y_continuous(breaks = c(seq(0,1,0.2)), limits = c(0,1), 
                       labels = scales::label_percent()) +
    scale_colour_manual(values = biovisr::cbf_palette(factor(cl_size_data$inflation)),
                        labels = divide_by_10,
                        guide = guide_legend(reverse = TRUE)) +
    theme_minimal()
}

ecdf_plot <- function(expt, threshold, dat) {
  dat <- dat |> 
    dplyr::filter(expt == {{expt}}, threshold == {{threshold}})
  dat |>
    arrange(cluster_size) |>
    group_by(inflation, cluster_size) |>
    summarise(
      total_nodes = sum(cluster_size), 
      .groups = "drop_last"
    ) |>
    mutate(
      frac = total_nodes / sum(total_nodes),
      ecdf = cumsum(frac)) |>
    ggplot(aes(x = cluster_size, y = ecdf, colour = inflation)) +
    geom_vline(xintercept = cmd_line_args$options$min_cluster_size,
               colour = "firebrick4", alpha = 0.75, linetype = "dashed") +
    geom_vline(xintercept = cmd_line_args$options$max_cluster_size,
               colour = "firebrick4", alpha = 0.75, linetype = "dashed") +
    geom_step() +
    scale_x_log10() +
    scale_y_continuous(breaks = c(seq(0,1,0.2)), limits = c(0,1), 
                       labels = scales::label_percent()) +
    scale_colour_manual(values = biovisr::cbf_palette(factor(cl_size_data$inflation)),
                        labels = divide_by_10,
                        guide = guide_legend(reverse = TRUE)) +
    labs(title = glue::glue("{expt}, threshold = {threshold}")) +
    theme_minimal()
}

output_plots <- function(dat, dat_summary, method) {
  pdf(file = file.path(plot_dir, paste0(method, "-clustering-summary-plots.pdf")))
  num_clusters_plot(dat_summary) |> print()
  num_non_singleton_nodes_plot(dat_summary) |> print()
  dat %>%
    split(.$expt) |> 
    map(cluster_size_histogram) |> 
    imap(function(x, expt){
      x + labs(title = expt)
    }) |> 
    map(print)
  # facetted ecdf plot
  facetted_ecdf(dat) |> print()
  # individual plots
  exp_by_threshold <- expand_grid(expt = levels(dat$expt), threshold = unique(dat$threshold))
  map2(exp_by_threshold$expt, exp_by_threshold$threshold, ecdf_plot, dat) |> 
    map(print)
  dev.off()
}

map(c("knn", "cor"), 
    function(x){ 
      output_plots(cl_size_data_by_method[[x]], cl_size_summary_by_method[[x]], x)
    }) |> invisible()

calc_diff <- function(dat_df, min_size, max_size) {
  lower <- dat_df$ecdf[ dat_df$cluster_size <= min_size ] |> tail(1)
  upper <- dat_df$ecdf[ dat_df$cluster_size <= max_size ] |> tail(1)
  dat_df |> dplyr::select(expt:inflation) |> 
    slice_head(n = 1) |> 
    bind_cols(tibble(ecdf_diff = upper - lower))
}

calc_ecdf_diff <- function(dat_df, method, 
                           min_cl_size = 50, max_cl_size = 1000) {
  ecdf_diff <- dat_df |> 
    arrange(cluster_size) |>
    group_by(expt, threshold, inflation, cluster_size) |>
    summarise(
      total_nodes = sum(cluster_size), 
      .groups = "drop_last"
    ) |>
    mutate(
      frac = total_nodes / sum(total_nodes),
      ecdf = cumsum(frac)
    ) %>%
    split(list(.$expt, .$threshold, .$inflation)) |> 
    map(\(x) {
      calc_diff(
        x, 
        min_size = min_cl_size,
        max_size = max_cl_size
      )
    }) |> 
    list_rbind() |> 
    mutate(method = method) |> 
    ungroup(threshold, inflation) |> 
    arrange(expt, desc(ecdf_diff))
  
  write_tsv(ecdf_diff, file = file.path(paste(method, "cluster_ecdf_diff.tsv", sep = "-")))
  
  return(ecdf_diff)
}

ecdf_diff <- imap(cl_size_data_by_method,
    function(x, name){ 
      calc_ecdf_diff(
        x, 
        name,
        cmd_line_args$options$min_cluster_size,
        cmd_line_args$options$max_cluster_size
      )
    })

top_ecdf_diff_by_expt_method <- 
  map(ecdf_diff, \(x){ slice_head(x, n = 1) }) |> 
  list_rbind() |> 
  arrange(expt, desc(ecdf_diff))

table_title <- glue::glue(
  "% of nodes in clusters larger than {cmd_line_args$options$min_cluster_size}",
  "and smaller than {cmd_line_args$options$max_cluster_size}",
  .sep = " "
)
filename <- glue::glue(
  "top_ecdf_diff_table", 
  "{cmd_line_args$options$min_cluster_size}",
  "{cmd_line_args$options$max_cluster_size}.html",
  .sep = "-")
top_ecdf_diff_by_expt_method |> 
  gt() |>
  tab_header(title = table_title) |> 
  fmt_percent(columns = c(ecdf_diff), decimals = 1) |>
  tab_style(
    style = cell_fill(color =  biovisr::cbf_palette()[5]),
    locations = cells_body(rows = method == "knn")
  ) |>
  tab_style(
    style = cell_fill(color = biovisr::cbf_palette()[6]),
    locations = cells_body(rows = method == "cor")
  ) |> gtsave(filename = filename)

top_ecdf_diff_by_expt_method |> 
  write_tsv(file = "top_ecdf_diff.tsv")

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
