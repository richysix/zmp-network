#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(
    "--samples_file", type = "character",
    default = "expt-sample-condition-initial.tsv", 
    help = "Name of the overall samples file [default %default]" 
  ),
  make_option("--debug", type = "logical", default = FALSE, action = "store_true",
              help = "Turns on debugging statements [default %default]" )
)

desc <- paste(
  '\nScript to create plots assessing the effect of filtering on each network',
  sep = "\n"
)
cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'edge-filtering-analysis.R',
    usage = "Usage: %prog [options]",
    description = desc),
  positional_arguments = 0
)

# load packages
packages <- c('tidyverse', 'patchwork', 'rprojroot', 'biovisr', 'miscr')
for (package in packages) {
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

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
  ) |> 
  arrange(n)
colour_palette <- biovisr::cbf_palette(sample_counts$expt,
                                       named = TRUE)

# set levels of expt
expts$expt <- factor(expts$expt, levels = levels(sample_counts$expt))

# function to load the correlation histogram data
load_cor_hist <- function(expt_name, filename) {
  read_tsv(
    filename, show_col_types = FALSE,
    col_names = c("bin_end", "count")
  ) |> 
    mutate(
      bin_start = bin_end - 0.1,
      expt = expt_name
    )
}

# load data
cor_hist_data <- 
  map(expts$expt,
      function(x) {
        filename <- file.path(paste(x, "tpm-orig.cor-hist.txt", sep = "-"))
        load_cor_hist(x, filename)
      }) |>
  list_rbind() |> 
  mutate(expt = fct_relevel(expt, levels(sample_counts$expt))) |> 
  inner_join(sample_counts, by = "expt")

plot_cor_hist <- function(cor_hist_data, bin_labels) {
  max_val <- max(cor_hist_data$count) |> log10() |> ceiling()
  all_cor_col <- cor_hist_data |> 
    mutate(count = count + 1) |> 
    ggplot(aes(x = bin_end, y = count)) +
    geom_col(aes(fill = sample_n_bin), position = position_nudge(x = -0.05)) +
    scale_y_log10(
      breaks = 10^(seq_len(max_val)),
      labels = scales::label_log()
    ) +
    scale_fill_manual(values = biovisr::cbf_palette(bin_labels, named = TRUE)) +
    facet_wrap(vars(expt)) +
    labs(x = "Spearman Correlation") +
    theme_minimal()
  return(all_cor_col)
}

all_cor_col <- plot_cor_hist(cor_hist_data, bin_labels)
miscr::output_plot(
  list(plot = all_cor_col, filename = file.path(plot_dir, "orig-cor-dist.pdf")),
  width = 12, 
  height = 8
)

# plot histogram for filtered data
cor_hist_data <- 
  map(expts$expt,
      function(x) {
        filename <- file.path(paste(x, "tpm-filtered-orig.cor-hist.txt", sep = "-"))
        load_cor_hist(x, filename)
      }) |>
  list_rbind() |> 
  mutate(expt = fct_relevel(expt, levels(sample_counts$expt))) |> 
  inner_join(sample_counts, by = "expt")

all_cor_col <- plot_cor_hist(cor_hist_data, bin_labels)
miscr::output_plot(
  list(plot = all_cor_col, filename = file.path(plot_dir, "filtered-cor-dist.pdf")),
  width = 12, 
  height = 8
)

 # L       Percentage of nodes in the largest component
 # D       Percentage of nodes in components of size at most 3 [-div option]
 # R       Percentage of nodes not in L or D: 100 - L -D
 # S       Percentage of nodes that are singletons
 # E       Fraction of edges retained (input graph has 87883982)
load_stats <- function(expt, file) {
  filename <- file.path(paste(expt, file, sep = "-"))
  read_tsv(filename, show_col_types = FALSE) |> 
    mutate(Expt = expt)
}

node_stats_plots <- function(stats_df, method) {
  plot_list <- list()
  filename_base <- file.path(plot_dir, method)
  singletons <- ggplot(data = stats_df) +
    geom_point(aes(x = Cutoff, y = S, colour = Expt)) +
    scale_colour_manual(values = colour_palette) +
    ylab("Singletons") +
    theme_minimal()
  miscr::output_plot(
    list(plot = singletons,
         filename = paste0(filename_base, '-stats-singletons.pdf'))
  )

  nodes <- ggplot(data = stats_df) +
    geom_point(aes(x = Cutoff, y = L, colour = Expt)) +
    scale_colour_manual(values = colour_palette) +
    ylab("Connected Nodes") +
    theme_minimal()
  miscr::output_plot(
    list(plot = nodes,
         filename = paste0(filename_base, '-stats-nodes.pdf'))
  )
  miscr::add_to_plot_list(
    plot_list,
    nodes,
    filename = paste0(filename_base, '-stats-nodes.pdf')
  )
  
  combined_plot <- nodes + singletons + plot_layout(guides = 'collect')
  miscr::output_plot(
    list(plot = combined_plot, 
         filename = paste0(filename_base, "-stats-nodes-singletons.pdf")),
    width = 12, height = 4.75
  )

  node_degrees <- ggplot(data = stats_df) +
    geom_hline(yintercept = 100, colour = "firebrick3", linetype = "dashed") +
    geom_pointrange(aes(x = Cutoff, y = NDmed, ymin = NDmed - NDiqr/2, 
                        ymax = NDmed + NDiqr/2, colour = Expt)) +
    scale_colour_manual(values = colour_palette) +
    theme_minimal()
  miscr::output_plot(
    list(plot = node_degrees, 
         filename = paste0(filename_base, "-stats-node-degrees.pdf")),
    width = 9.6,
  )
  
  node_degrees_close_up <- ggplot(data = stats_df) +
    geom_hline(yintercept = 100, colour = "firebrick3", linetype = "dashed") +
    geom_point(aes(x = Cutoff, y = NDmed, colour = Expt)) +
    scale_colour_manual(values = colour_palette) +
    lims(y = c(NA, 200)) +
    theme_minimal()
  miscr::output_plot(
    list(plot = node_degrees_close_up, 
         filename = paste0(filename_base, "-stats-node-degrees-close-up.pdf")),
    width = 9.6,
  )
  all_components_plot(stats_df, method)
}

all_components_plot <- function(stats_df, method) {
  stats_long <- stats_df |> 
    dplyr::select(L:S, Cutoff, Expt) |> 
    pivot_longer(-c(Expt, Cutoff), names_to = "category", values_to = "nodes") |> 
    mutate(category = factor(category, levels = c("L", "D", "S", "R")))
  
  categories <- c(
    L = "Largest Component",
    S = "Singletons",
    R = "Remainder",
    D = "Disconnected"
  )
  all_components <- ggplot(stats_long, aes(x = Cutoff, y = nodes, colour = Expt)) +
    geom_point() +
    facet_wrap(vars(category), nrow = 1, scales = "free_y",
               labeller = labeller(category = categories)) +
    scale_colour_manual(values = colour_palette) +
    theme_minimal()
  if (method == "knn") {
    all_components <- all_components +
      scale_x_reverse()
  }
  
  filename_base <- file.path(plot_dir, method)
  miscr::output_plot(
    list(plot = all_components, 
         filename = paste0(filename_base, "-stats-all-components.pdf")),
    width = 12,
    height = 4.75
  )
}

cor_stats <- map(expts$expt, \(x) load_stats(x, "tpm-filtered-t20.vary-cor-stats.tsv")) |>
  list_rbind() |>
  mutate(Expt = fct_relevel(Expt, levels(sample_counts$expt)))
node_stats_plots(cor_stats, "cor")

# KNN stats
knn_stats <- map(expts$expt,
                 \(x) load_stats(x, "tpm-filtered-t20.vary-knn-stats.tsv")) |> 
  list_rbind() |> 
  rename(Cutoff = kNN)
node_stats_plots(knn_stats, "knn")

load_node_degrees <- function(filename) {
  method <- ifelse(grepl("tpm-filtered-t[0-9]+-k[0-9]+", filename), "knn", "cor")
  if (method == "knn") {
    threshold <- str_remove(filename, "^.*tpm-filtered-t[0-9]+-k") |>
      str_remove("\\.stats\\.tsv") |>
      as.integer()
  } else {
    threshold <- str_remove(filename, "^.*tpm-filtered-t") |>
      str_remove("\\.stats\\.tsv") |>
      as.integer()
  }
  expt <- str_remove(filename, "-tpm.*$")
  read_tsv(filename, show_col_types = FALSE) |> 
    mutate(
      Expt = expt,
      Method = method,
      Threshold = threshold
    )
}

stats_files <- list.files(
  pattern = "*\\.stats\\.tsv", recursive = TRUE
)

node_degree_data <- purrr::map(stats_files, load_node_degrees) |>
  list_rbind() |>
  mutate(
    Expt = factor(Expt, levels = levels(sample_counts$expt))
  )

facetted_histograms <- function(degree_df) {
  method <- degree_df$Method[1]
  if (method == "knn") {
    degree_df <- degree_df |> 
      mutate(
        Threshold = factor(
          Threshold,
          levels = unique(degree_df$Threshold) |> 
            sort(decreasing = TRUE))
      )
  }
  median_by_expt_by_threshold <- degree_df |> 
    dplyr::filter(degree > 0) |> 
    group_by(Expt, Threshold) |> 
    summarise(
      median_degree = median(degree),
      .groups = "drop"
    )
  node_degree_distribution <- ggplot(data = degree_df,
                                     aes(x = degree, fill = Expt)) +
    geom_histogram(binwidth = 10, boundary = 1, show.legend = FALSE) +
    facet_grid(rows = vars(Expt), cols = vars(Threshold)) +
    geom_vline(data = median_by_expt_by_threshold,
               aes(xintercept = median_degree),
               colour = "black") +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "firebrick3") +
    scale_fill_manual(values = colour_palette) +
    theme_minimal()

  out_file <- file.path(
    plot_dir, paste(method, "stats-node-degree-distribution.pdf",
                    sep = "-")
  )
  miscr::output_plot(
    list(plot = node_degree_distribution, 
         filename = out_file),
    width = 12,
    height = 4.75
  )
  node_degree_distribution_close_up <- 
    ggplot(data = degree_df, aes(x = degree, fill = Expt)) +
    geom_histogram(binwidth = 10, boundary = 1, show.legend = FALSE) +
    facet_grid(rows = vars(Expt), cols = vars(Threshold)) +
    geom_vline(data = median_by_expt_by_threshold,
               aes(xintercept = median_degree),
               colour = "black") +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "firebrick3") +
    scale_fill_manual(values = colour_palette) +
    scale_x_continuous(limits = c(-10,200)) +
    theme_minimal()
  
  out_file <- file.path(
    plot_dir, paste(method, "stats-node-degree-distribution-close-up.pdf",
                    sep = "-")
  )
  miscr::output_plot(
    list(plot = node_degree_distribution_close_up, 
         filename = out_file),
    width = 12,
    height = 4.75
  )
}

# cor stats
node_degree_data |>
  dplyr::filter(Method == "cor") |>
  facetted_histograms()

# knn
node_degree_data |>
  dplyr::filter(Method == "knn") |>
  facetted_histograms()
