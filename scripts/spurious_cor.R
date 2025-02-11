#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to investigate count data", sep = "\n")

if (interactive()) {
  cmd_args <- c("zmp_ph204", "zmp_ph192", "zmp_ph213", "zmp_ph238", "zmp_ph250",
                "zmp_ph46", "zmp_ph71") # add test files here
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

packages <- c('tidyverse', 'rnaseqtools', 'patchwork', 'ggridges', 'ggtext')
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

# load TPMs
# load TPM values from file for a given expt
load_tpms <- function(expt) {
  expt_dir <- file.path("nf", "results", expt)
  tpms_file <- file.path(expt_dir, paste(expt, "all-tpm.tsv", sep = "-"))
  read_tsv(tpms_file, show_col_types = FALSE)
}

# load sample info from file for an expt
load_samples <- function(expt) {
  expt_dir <- file.path("nf", "results", expt)
  sample_file <- file.path(expt_dir, "samples.tsv")
  read_tsv(sample_file, show_col_types = FALSE)
}

# get only the tpm values, transpose, select all the columns where 
# all values in that column are zero and remove them from the original tibble
# and the transposed matrix
filter_all_zeros_and_transpose_tpms <- function(tpms, expt) {
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

# takes a tibble of tpms and a transposed matrix of tpms
# and adds metadaa to the tibble
# variance and the variance decile
# num of zero values
# mean_tpm and mean_tpm decile
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

# arrange the sample data by the condition column
# based on the genotype values
arrange_by_condition <- function(samples) {
  if (all(unique(samples$condition) %in% c("wt", "het", "hom"))) {
    condition_levels <- c("wt", "het", "hom")
  } else if (all(unique(samples$condition) %in% c("sib", "mut"))) {
    condition_levels <- c("sib", "mut")
  } else if (all(grepl("(wt|het|hom)_tfap2a_(wt|het|hom)_tfap2c", unique(samples$condition)))) {
    condition_levels <- c(
      "wt_tfap2a_wt_tfap2c",
      "wt_tfap2a_het_tfap2c",
      "het_tfap2a_wt_tfap2c",
      "het_tfap2a_het_tfap2c",
      "wt_tfap2a_hom_tfap2c",
      "hom_tfap2a_wt_tfap2c",
      "het_tfap2a_hom_tfap2c",
      "hom_tfap2a_het_tfap2c",
      "hom_tfap2a_hom_tfap2c"
    )
  } else {
    cat("Sample data:\n", readr::format_tsv(samples))
    stop("Couldn't determine the genotype levels")
  }
  
  samples |>
    mutate(condition = factor(condition, levels = condition_levels)) |>
    arrange(condition) |>
    mutate(sample = factor(sample, levels = unique(sample)))
}

# plots of numbers of zeros by variance decile
# number of zeros by expression decile (mean tpm)
# number of zeros by mean tpm
# heatmap of genes with zeros for more than 1/4 samples
num_zero_plots <- function(tpms_df, samples, num_zeros_cutoff) {
  plot_list <- list(
    "num_zeros_by_var_q" = NULL,
    "num_zeros_by_expr_q" = NULL,
    "num_zeros_by_mean_expr" = NULL,
    "too_many_zeros_heatmap" = NULL
  )
  # Look at num zeros vs variance by percentile
  plot_list[["num_zeros_by_var_q"]] <- tpms_df |> ggplot(aes(x = num_zeros)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(vars(var_quantile), ncol = 5) +
    labs(
      title = "Genes with the lowest variances have the most numbers of zero values",
      subtitle = "Histogram of zero TPM values facetted by variance percentile"
    )
  # Look at num zeros vs expression percentile
  plot_list[["num_zeros_by_expr_q"]] <- tpms_df |> ggplot(aes(x = num_zeros)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(vars(expr_quantile), ncol = 5) +
    labs(
      title = "Genes with the lowest mean expression levels have the most numbers of zero values",
      subtitle = "Histogram of zero TPM values facetted by mean expression percentile"
    )
  # Look at num_zeros vs mean tpm
  cutoffs <- filter(tpms_df, num_zeros >= num_zeros_cutoff) |> 
    summarise(
      max_tpm = max(mean_tpm),
      pc_95 = quantile(mean_tpm, 0.95)
    )
  x_max <- ceiling(cutoffs$max_tpm[1]/10)*10
  
  cutoffs <- cutoffs |> tidyr::pivot_longer(
      cols = everything(), names_to = "type", values_to = "value"
    ) |> 
    mutate(
      y_max = nrow(samples),
      text = glue::glue("x = {round(value, digits = 1)}"),
      hjust = c(-1, 1),
      nudge = x_max*0.1*hjust
    )
  plot_list[["num_zeros_by_mean_expr"]] <- tpms_df |> 
    ggplot(aes(x = mean_tpm, y = num_zeros)) +
    geom_point(position = position_jitter(width = 0, height = 0.3)) +
    geom_vline(
      data = cutoffs,
      aes(xintercept = value, linetype = type),
      colour = "firebrick3") +
    geom_text(
      data = cutoffs, 
      aes(x = value, y = y_max, label = text),
      nudge_x = cutoffs$nudge) +
    lims(x = c(0,x_max))
  
  samples <- arrange_by_condition(samples) |> 
    mutate(group = 1)
  sample_metadata_plot <-
    biovisr::df_heatmap(samples, x = 'sample', y = "group",
               fill = "condition", fill_palette = biovisr::cbf_palette(samples$condition),
               colour = "grey50", size = 0.5,
               xaxis_labels = TRUE, yaxis_labels = FALSE,
               na.translate = FALSE
    ) + guides(fill = "none") +
    theme(
      axis.title = element_blank(),
      axis.title.y = element_blank()
    )
  
  expr_heatmap <- 
    filter(tpms_df, num_zeros >= num_zeros_cutoff) |>
    arrange(desc(mean_tpm)) |>
    # order samples by condition
    dplyr::select(all_of(samples$sample)) |>
    mutate(across(everything(), log10), 
           across(everything(), \(x){ ifelse(is.infinite(x), NA, x) } )) |>
    biovisr::matrix_heatmap(xaxis_labels = FALSE, yaxis_labels = FALSE)
  plot_list[["too_many_zeros_heatmap"]] <- 
    expr_heatmap + sample_metadata_plot + plot_layout(nrow = 2, heights = c(20,1))
  return(plot_list)
}

# plot a histogram of correlation values for the supplied tpms data frame
plot_cor_hist <- function(tpms_df) {
  cor_mat <- tpms_df |>
    dplyr::select(starts_with("zmp")) |>
    t() |> cor(method = "spearman")
  rownames(cor_mat) <- tpms_df$GeneID
  colnames(cor_mat) <- tpms_df$GeneID
  cor_long <- as_tibble(cor_mat, rownames = "Gene1") |>
    pivot_longer(cols = -Gene1, names_to = "Gene2", values_to = "cor")
  cor_long |> filter(Gene1 != Gene2) |>
    ggplot(aes(x = cor)) +
    geom_histogram(binwidth = 0.1, boundary = 0) +
    scale_y_log10() +
    lims(x = c(-1,1))
}

calculate_cor_values_vs_all_genes <- function(tpms_df, all_tpms, label, expt) {
  a <- tpms_df |>
    dplyr::select(starts_with("zmp")) |>
    t()
  colnames(a) <- tpms_df$GeneID
  b <- all_tpms |>
    dplyr::select(starts_with("zmp")) |>
    t()
  colnames(b) <- all_tpms$GeneID
  b <- b[ , c(colnames(a), setdiff(colnames(b), colnames(a))) ]
  # calculate cor vs all other genes
  c <- cor(a, b, method = "spearman")
  # rearrange the matrix with the a genes first
  c <- c[ , c(rownames(c), setdiff(colnames(c), rownames(c))) ]
  # set cells where gene1 == gene2 to NA
  num_rows <- nrow(tpms_df)
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
  tibble(
    Gene1 = rownames(c)[ row_i ],
    Gene2 = colnames(c)[ col_i ],
    cor = c[ indices ],
  ) |>
    arrange(cor) |>
    write_tsv(file = file.path(output_dir, paste0(paste("cor-by", label, expt, sep = "-"), ".tsv")))

  # subset matrix to just where coli > rowi
  indices <- seq_along(c)
  col_i <- as.integer((indices - 1) / nrow(c)) + 1
  row_i <- ((indices - 1) %% nrow(c)) + 1
  to_keep <- col_i > row_i
  c <- c[ to_keep ]
  # return counts for histogram
  cut(c, breaks = seq(-1,1,0.1)) |> 
    table()
}

cor_by_mean_tpm <- function(all_tpms, expt) {
  hist_data <- split(all_tpms, all_tpms$expr_quantile) |>
    imap(.f = function(x, q){ calculate_cor_values_vs_all_genes(x, all_tpms, paste0("expr-quantile-", q), expt) } ) |> 
    imap(.f = function(x, q){
      tibble(
        range = names(x),
        count = x,
        midpoint = seq(-0.95,0.95,0.1),
        expr_quantile = q
      )
    }) |> 
    list_rbind()
  
  filter_by_tpm_plot <- hist_data |>
    ggplot(aes(x = midpoint, y = count)) +
    geom_col() +
    facet_wrap(vars(expr_quantile), ncol = 5) +
    scale_y_log10() +
    lims(x = c(-1,1)) +
    labs(title = "The distribution of correlations looks more normal for more highly expressed genes")
  
  filter_by_tpm_table <- filter(hist_data, midpoint == 0.95)
  
  return(
    list(
      plot = filter_by_tpm_plot,
      table = filter_by_tpm_table
    )
  )
}

cor_by_num_zeros <- function(tpms_df, expt, high_cor_threshold = 100) {
  zero_levels <- unique(tpms_df$num_zeros) |> 
    sort() |> as.character()
  num_genes <-   split(tpms_df, tpms_df$num_zeros) |>
    imap(\(x, idx) { tibble(
      num_zeros = idx,
      num_genes = nrow(x)
    )}) |> 
    list_rbind() |> 
    mutate(num_zeros = factor(num_zeros, levels = zero_levels))

  num_genes_lookup <- paste0(num_genes$num_zeros, paste0(" (n = ", num_genes$num_genes, ")"))
  names(num_genes_lookup) <- num_genes$num_zeros
  
  hist_data <- split(tpms_df, tpms_df$num_zeros) |>
    imap(.f = function(x, q){ calculate_cor_values_vs_all_genes(x, tpms_df, paste0("num-zeros-", q), expt) } ) |> 
    imap(.f = function(x, q){
      tibble(
        range = names(x),
        count = as.integer(x),
        midpoint = seq(-0.95,0.95,0.1),
        num_zeros = q
      )
    }) |> 
    list_rbind() |> 
    # order num_zeros numerically
    mutate(
      num_zeros = factor(num_zeros, levels = zero_levels),
      count = case_when(
        count == 0 ~ 1,
        TRUE ~ count
      )
    )
  
  n <- max(tpms_df$num_zeros)
  if (n > 23) {
    n <- max(tpms_df$num_zeros)
    gap = floor((n - 12)/3)
    zero_subset <- c(
      0:5, 
      (5 + gap + 1 - 3):(5 + gap + 3),
      (5 + 2*gap + 1 - 3):(5 + 2*gap + 3),
      (n - 5):n
    ) |>
      as.character()
  } else {
    zero_subset <- seq(0, n) |>
      as.character()
  }

    
  maxy <- ceiling(max(hist_data$count) |> log10())
  filter_by_zeros_plot <- hist_data |>
    filter(num_zeros %in% zero_subset) |>
    # label top bar
    mutate(bar_cat = case_when(
      range == "(0.9,1]" ~ "high",
      TRUE ~ "other"
    )) |> 
    ggplot(aes(x = midpoint, y = count, fill = bar_cat)) +
    geom_hline(yintercept = high_cor_threshold, colour = "firebrick4") +
    geom_col(show.legend = FALSE) +
    facet_wrap(vars(num_zeros), ncol = 6,
               labeller = labeller(num_zeros = num_genes_lookup)) +
    scale_y_log10(breaks = 10^seq(1,maxy)) +
    scale_fill_manual(values = c("darkblue", "grey50")) +
    lims(x = c(-1,1)) +
    labs(
      title = paste0("Genes with lots of zeros have the majority of spuriously ",
                     "<span style='color:#00008B'><strong>high</strong></span> correlations"),
    ) +
    theme_minimal() +
    theme(
      title = ggtext::element_markdown()
    )
  
  filter_by_zeros_table <- filter(hist_data, midpoint == 0.95) |> 
    inner_join(num_genes, by = join_by(num_zeros == num_zeros)) |> 
    mutate(num_zeros = as.character(num_zeros) |> as.integer(),
           cumul_genes = cumsum(num_genes))
  # find indices where count >= threshold after lowest value
  min_val <- which.min(filter_by_zeros_table$count)
  idx <- which(filter_by_zeros_table$count >= high_cor_threshold & seq_along(filter_by_zeros_table$count) > min_val)
  zeros_above_threshold <- filter_by_zeros_table$num_zeros |>
    magrittr::extract(idx)

  return(
    list(
      plot = filter_by_zeros_plot,
      table = filter_by_zeros_table,
      zeros_above_threshold = zeros_above_threshold
    )
  )
}

filter_by_mean_tpms_plot <- function(tpms_df, num_zeros_cutoff) {
  # run filter_by_mean_tpm for increasing tpm threshold
  quantiles <- seq(0, 1, 0.01)
  tpm_thresholds <- quantile(tpms_df$mean_tpm, quantiles)
  filter_df <- tibble(
    tpm_threshold = tpm_thresholds,
    tpm_quantile = quantiles,
    tn = map_int(tpm_thresholds, \(x) { true_negative(x, tpms_df, num_zeros_cutoff) }),
    fn = map_int(tpm_thresholds, \(x) { false_negative(x, tpms_df, num_zeros_cutoff) }),
    fp = map_int(tpm_thresholds, \(x) { false_positive(x, tpms_df, num_zeros_cutoff) }),
    tp = map_int(tpm_thresholds, \(x) { true_positive(x, tpms_df, num_zeros_cutoff) })
  )
  axis_breaks <- seq(0,1,0.2)
  xaxis_labels <- tpm_thresholds[paste0(as.character(axis_breaks*100), "%")] |> 
    round(digits = 1)
  filter_df |>
    mutate(
      tn_pc = tn / sum(tpms_df$num_zeros > num_zeros_cutoff),
      fn_pc = fn / sum(tpms_df$num_zeros <= num_zeros_cutoff),
      fp_pc = fp / sum(tpms_df$num_zeros > num_zeros_cutoff),
      tp_pc = tp / sum(tpms_df$num_zeros <= num_zeros_cutoff)
    ) |>
    dplyr::select(c(tpm_threshold, tpm_quantile, tn_pc:tp_pc)) |>
    pivot_longer(
      cols = tn_pc:tp_pc,
      names_to = "group"
    ) |>
    ggplot(aes(x = tpm_quantile, y = value, colour = group)) +
    geom_line() +
    scale_y_continuous(breaks = axis_breaks) +
    scale_x_continuous(breaks = axis_breaks, labels = xaxis_labels)
}

true_negative <- function(cutoff, tpms_df, num_zeros_cutoff) {
  sum(tpms_df$mean_tpm < cutoff & tpms_df$num_zeros > num_zeros_cutoff)
}
false_negative <- function(cutoff, tpms_df, num_zeros_cutoff) {
  sum(tpms_df$mean_tpm < cutoff & tpms_df$num_zeros <= num_zeros_cutoff)
}
false_positive <- function(cutoff, tpms_df, num_zeros_cutoff) {
  sum(tpms_df$mean_tpm >= cutoff & tpms_df$num_zeros > num_zeros_cutoff)
}
true_positive <- function(cutoff, tpms_df, num_zeros_cutoff) {
  sum(tpms_df$mean_tpm >= cutoff & tpms_df$num_zeros <= num_zeros_cutoff)
}

cor_matrix_plot_by_num_zeros <- function(tpms_df, expt, samples, num_genes = 400) {
  cor_mat <- tpms_df |>
    dplyr::select(starts_with("zmp")) |>
    t() |> cor(method = "spearman")
  colnames(cor_mat) <- tpms_df$GeneID
  rownames(cor_mat) <- tpms_df$GeneID
  # sample equal number of genes from each low, medium and high num_zeros
  sample_genes <- function(zeros, tpms_df, n) {
    filtered_tpms_df <- tpms_df |>
      dplyr::filter(num_zeros == zeros)
    n <- ifelse(nrow(filtered_tpms_df) < n, nrow(filtered_tpms_df), n)
    filtered_tpms_df |>
      dplyr::slice_sample(n = n)
  }
  # pick subsets to sample
  n_samples <- nrow(samples)
  midpoint <- floor(n_samples/2)
  subsets <- c(0,1,midpoint, midpoint + 1, n_samples - 2, n_samples - 1)
  tpms_sample <- purrr::map(subsets,
                            function(x){ sample_genes(x, tpms_df, num_genes) }) |> 
    list_rbind() |> 
    arrange(num_zeros, mean_tpm)
  cor_mat <- cor_mat[tpms_sample$GeneID, tpms_sample$GeneID]
  # reclaim memory
  gc() |> invisible()
  
  cor_heatmap <- 
    biovisr::matrix_heatmap(
      cor_mat, x_title = "Gene2", y_title = "Gene1", fill_title = "Cor", 
      xaxis_labels = FALSE, yaxis_labels = FALSE, base_size = 30
    ) +
    scale_fill_viridis_b(breaks = c(-1,-0.5,0,0.5,0.9,1), direction = -1) +
    guides(fill = guide_legend(override.aes = list(size = 2)))
  
  num_zeros_y <- tpms_sample |>
    mutate(GeneID = factor(GeneID, levels = rev(GeneID))) |>
    ggplot(aes(x = 1, y = GeneID, fill = num_zeros)) +
    geom_tile() +
    scale_fill_viridis_c(option = "rocket", direction = -1) +
    guides(fill = guide_legend(override.aes = list(size = 2))) +
    theme_void() +
    theme(
      legend.position = "left", 
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 25)
    )
  
  num_zeros_x <- tpms_sample |>
    mutate(GeneID = factor(GeneID, levels = GeneID)) |>
    ggplot(aes(y = 1, x = GeneID, fill = num_zeros)) +
    geom_tile() +
    scale_fill_viridis_c(option = "rocket", direction = -1) +
    theme_void() +
    theme(legend.position = "none")
  
  return(list(
    cor_heatmap = cor_heatmap,
    num_zeros_x = num_zeros_x,
    num_zeros_y = num_zeros_y
  ))
}

effect_of_zeros <- function(expt) {
  tpms <- load_tpms(expt)
  samples <- load_samples(expt)
  tpm_list <- filter_all_zeros_and_transpose_tpms(tpms, expt)
  tpms_no_all_zeros <- add_metadata(tpm_list[[1]], tpm_list[[2]])
  # filter by num_zeros
  num_zeros_cutoff <- ceiling(nrow(samples) * 0.25)
  # no more than 1/4 of samples can be zero
  filtered_by_zeros <- tpms_no_all_zeros |> 
    filter(num_zeros < num_zeros_cutoff)
  cat(glue::glue("Expt {expt}: {nrow(filtered_by_zeros)} genes have zeros in fewer then 1/4 of samples"), "\n")
  
  plots <- num_zero_plots(tpms_no_all_zeros, samples, num_zeros_cutoff)
  pdf(file.path("plots", paste(expt, "num-zero-plots.pdf", sep = "-")))
  walk(plots, print)
  filter_by_tpms <- filter_by_mean_tpms_plot(tpms_no_all_zeros, num_zeros_cutoff)
  print(filter_by_tpms)
  
  # filter by mean_tpm, sample genes and plot hist of cor values
  # sample genes with lots of zeros and plot hist of cor values
  sub_sample <- tpms_no_all_zeros |>
    filter(num_zeros > num_zeros_cutoff) |>
    slice_sample(prop = 0.316)
  cor_hist <- plot_cor_hist(sub_sample) +
    labs(title = "Genes with more zeros are skewed towards high correlation values")
  print(cor_hist)
  # sample genes and plot hist of cor values
  sub_sample <- filtered_by_zeros |> slice_sample(prop = 0.316)
  cor_hist <- plot_cor_hist(sub_sample) +
    labs(title = "Genes filtered by zero don't have an excess of extreme correlation values")
  print(cor_hist)

  cor_by_mean_tpm_list <- cor_by_mean_tpm(tpms_no_all_zeros, expt)
  print(cor_by_mean_tpm_list[[1]])
  
  mean_tpm_dist <- filtered_by_zeros |> 
    ggplot(aes(x = mean_tpm)) +
    geom_histogram() +
    scale_x_log10() +
    labs(title = "Mean TPM distribution after filtering for zeros in fewer than 1/4 of the samples")
  print(mean_tpm_dist)
  
  high_cor_threshold <- 25
  cor_by_num_zeros_list <- cor_by_num_zeros(tpms_no_all_zeros, expt, high_cor_threshold)
  print(cor_by_num_zeros_list$plot)
  
  mean_tpm_filtered_dist <- tpms_no_all_zeros |> 
    filter(!(num_zeros %in% cor_by_num_zeros_list$zeros_above_threshold)) |>
    ggplot(aes(x = mean_tpm)) +
    geom_histogram() +
    scale_x_log10(labels = scales::comma) +
    labs(title = "Mean TPM distribution after filtering for genes above the zeros threshold")
  print(mean_tpm_filtered_dist)
  dev.off()
  
  matrix_heatmap_list <- 
    cor_matrix_plot_by_num_zeros(tpms_no_all_zeros, expt, samples, num_genes = 400)
  pdf(file = glue::glue("plots/cor-vs-num-zeros-heatmap-{expt}.pdf"), 
      width = 30, height = 25)
  print(
    plot_spacer() + 
    matrix_heatmap_list$num_zeros_x + 
    matrix_heatmap_list$num_zeros_y + 
    matrix_heatmap_list$cor_heatmap + 
    patchwork::plot_layout(nrow = 2, widths = c(1,20), heights = c(1,20))
  )
  dev.off()
  
  dplyr::select(cor_by_num_zeros_list$table, -midpoint) |>
    write_tsv(file = file.path(output_dir, glue::glue("cor-by-num-zeros-hist-{expt}.tsv")))
  
  total_genes <- 
    filter(cor_by_num_zeros_list$table,
           num_zeros >= min(as.integer(cor_by_num_zeros_list$zeros_above_threshold))) |> 
    pull(num_genes) |> sum()
  zeros_above_threshold <- 
    glue::glue_collapse(cor_by_num_zeros_list$zeros_above_threshold, ", ", 
                        last = " and ")
  cat(
    glue::glue("Expt {expt}: genes with {zeros_above_threshold} zeros",
               "({total_genes} genes) have more than {high_cor_threshold} spurious",
               "correlations with other genes.", .sep = " "), "\n")
  
  tpms_no_all_zeros |> 
    magrittr::set_colnames(
      sub("(zmp_ph.*)", "\\1 tpm", colnames(tpms_no_all_zeros))
    ) |>
    mutate(
      filter_by_zeros = num_zeros > num_zeros_cutoff,
      filter_by_tpm = mean_tpm < quantile(mean_tpm, 0.25)
    ) |>
    relocate(filter_by_zeros:filter_by_tpm, .after = expr_quantile) |>
    arrange(mean_tpm) |>
    write_tsv(file = file.path(output_dir, paste(expt, "filtering_zeros.tsv", sep = "-")))
}

output_dir <- "spurious-correlations"
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
