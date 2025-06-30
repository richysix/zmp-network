#!/usr/bin/env Rscript

## BASED ON CODE FROM https://github.com/sarbal/CoExpNets

library('optparse')

option_list <- list(
  make_option(
    "--annotation",
    default = "../reference/danio_rerio-e99-annotation.test.txt",
    help="Gene annotation file [default %default]"
  ),
  make_option(
    "--go_annotation",
    default = "../reference/danio_rerio_e109_go.test.txt",
    help="GO annotation file [default %default]"
  ),
  make_option(
    "--expts_file", type = "character",
    default = "expts.txt",
    help = "Name of the overall samples file [default %default]"
  ),
  make_option(
    "--samples_file", type = "character",
    default = "expt-sample-condition-tfap2-plus.tsv",
    help = "Name of the overall samples file [default %default]"
  ),
  make_option(
    "--num_orderings", type = "integer",
    default = 10,
    help = "Num of different orderings of the input files to do [default %default]"
  ),
  make_option(c("-d", "--debug"), action="store_true", default=FALSE,
              help="Print extra output [default %default]")
)

desc <- paste('Script to aggregate coexpression networks', sep = "\n")

# testing
if (interactive()) {
  cmd_args <- c('zmp_ph46-tpm-filtered-orig.test.mat', 
                'zmp_ph250-tpm-filtered-orig.test.mat',
                'zmp_ph204-tpm-filtered-orig.test.mat',
                'zmp_ph192-tpm-filtered-orig.test.mat',
                'zmp_ph238-tpm-filtered-orig.test.mat',
                'zmp_ph213-tpm-filtered-orig.test.mat',
                'zmp_ph71-tpm-filtered-orig.test.mat')
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list,
    description = desc,
    usage = "Usage: %prog [options] input_files..." ),
  args = cmd_args,
  positional_arguments = c(2, Inf)
)

packages <- c('tidyverse', 'vroom', 'pryr')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# function to load matrix
load_matrix <- function(file_name) {
  mat <- vroom::vroom(file_name, delim = "\t", col_names = TRUE,
                      show_col_types = FALSE, progress = FALSE) |>
    as.matrix()
  idx <- mat[, "dummy"]
  mat <- mat[, -1]
  
  tab_file <- sub(".mat", ".tab", file_name)
  tab_info <- read_tsv(tab_file, col_names = c("idx", "GeneID"),
                       show_col_types = FALSE)
  gene_id_lookup <- tab_info$GeneID
  names(gene_id_lookup) <- tab_info$idx
  
  colnames(mat) <- gene_id_lookup[ sub("^X", "", colnames(mat)) ] |> unname()
  rownames(mat) <- colnames(mat)
  return(mat)
}

# function to rank edges and standardise from 0 to 1
rank_and_standardise <- function(mat, absolute = FALSE) {
  # set self-edges to zero
  gene_ids <- colnames(mat)
  mat <- as.matrix(mat)
  diag(mat) <- 0
  if (absolute) {
    mat <- rank(abs(mat), na.last = "keep", ties.method = "min") |> 
      matrix(nrow = nrow(mat), ncol = ncol(mat))
  } else {
    mat <- rank(mat, na.last = "keep", ties.method = "min") |> 
      matrix(nrow = nrow(mat), ncol = ncol(mat))
  }
  rownames(mat) <- gene_ids
  colnames(mat) <- gene_ids
  mat/max(mat, na.rm = TRUE)
}

# # function to run GBA on subset network
# run_gba_on_subset <- function(mat, gene_id2term_id, go_annotation) {
#   annotations <- make_annotations(gene_id2term_id, colnames(mat), 
#                                   unique(go_annotation$TermID))
#   annotations <- filter_network_cols(
#     annotations,
#     min = cmd_line_args$options$min.term.size,
#     max = cmd_line_args$options$max.term.size
#   )
#   
#   Anno_groups_voted <- run_GBA(
#     mat,
#     annotations,
#     min = cmd_line_args$options$min.term.size,
#     max = cmd_line_args$options$max.term.size
#   )
#   return(Anno_groups_voted)
# }

# function to output a network to disk
output_network <- function(network, idx_lookup, filebase) {
  # set self-edges to zero
  diag(network) <- 0
  gene_ids <- colnames(network)
  # network <- 
  network |>
    as_tibble(rownames = "GeneID") |>
    # run GBA script expects node ids. convert with idx_lookup
    magrittr::set_colnames(c("dummy", idx_lookup[gene_ids])) |>
    dplyr::mutate(dummy = idx_lookup[dummy]) |>
    vroom_write(file = paste0(filebase, ".mat.gz"), progress = FALSE)
  # convert to nodes file
  # network |>
  #   dplyr::rename(source = dummy) |>
  #   tidyr::pivot_longer(cols = -source, names_to = "target", values_to = "weight") |>
  #   dplyr::filter(weight > 0) |>
  #   dplyr::mutate(edge_idx = consecutive_id(source, target)) |>
  #   dplyr::relocate(edge_idx) |>
  #   vroom_write(file = paste0(filebase, ".edges.tsv.gz"))
  return(NULL)
}

# function to aggregate 2 networks by addition
# adds network m2 to m1
aggregate_networks <- function(m1, m2, idx_lookup,
                               gene_id2term_id, go_annotation, outfile_base) {
  in_common <- match(colnames(m2), colnames(m1))
  f1 <- !is.na(in_common) # remove NAs
  f2 <- in_common[f1]
  # extract common subset of network
  m1[f2, f2] <- m1[f2, f2] + m2
  
  # re-rank and normalise
  m1 <- rank_and_standardise(m1)
  # output aggregated network
  output_network(m1, idx_lookup, outfile_base)
  
  if (cmd_line_args$options$debug) {
    print(m1[1:5,1:5])
  }
  return(m1)
}

# function to load new network and add to aggregate
add_network_to_aggregate <- function(file_name, agg, idx_lookup,
                                     gene_id2term_id, go_annotation, outfile_base) {
  # load network
  net <- load_matrix(file_name)
  net <- rank_and_standardise(net, absolute = TRUE)
  if (cmd_line_args$options$debug) {
    cat(glue::glue("Mem after loading/ranking input network ({file_name}) = {mem_used()/1024/1024} M"), "\n")
  }
  # add to aggregate network
  agg <- aggregate_networks(agg, net, idx_lookup, gene_id2term_id,
                            go_annotation, outfile_base)
  return(agg)
}

make_aggregate_network <- function(order_idx, ordering, annotation, gene_id2term_id, go_annotation) {
  # make base network
  num_genes <- length(annotation$GeneID)
  network <- matrix(0, nrow = num_genes, ncol = num_genes)
  colnames(network) <- annotation$GeneID
  rownames(network) <- colnames(network)
  
  idx_lookup <- annotation$idx |> 
    magrittr::set_names(annotation$GeneID)

  i <- 1
  for (file in ordering[[order_idx]]) {
    outfile_base <- glue::glue("agg-ordering{order_idx}-networks{i}")
    network <- add_network_to_aggregate(file, network,
                                        idx_lookup, gene_id2term_id, 
                                        go_annotation, outfile_base)
    if (cmd_line_args$options$debug) {
      cat(glue::glue("Mem after {i} files, (Round: {order_idx}) = {mem_used()/1024/1024} M"), "\n")
    }
    i <- i + 1
  }

  # set self-edges to zero
  idx <- seq(1, length(network), ncol(network) + 1)
  network[ idx ] <- 0
  # prune aggregate network and output
  top_decile <- quantile(network, 0.9)
  network[ network < top_decile ] <- 0
  # remove any genes where all values are zero
  to_keep <- colSums(network == 0) != num_genes
  network <- network[ to_keep, to_keep ]
  output_network(network, idx_lookup, 
                 glue::glue("agg-ordering{order_idx}-networks{i}"))
  return(invisible(NULL))
}

# load all gene IDs
annotation <- read_tsv(cmd_line_args$options$annotation,
                       col_names = c("GeneID", "Chr", "Start", "End", "Strand",
                                     "Biotype", "Name", "Description")) |> 
  dplyr::mutate(idx = consecutive_id(GeneID))

tibble(
  node_idx = annotation$idx,
  gene_id = annotation$GeneID
) |> vroom_write(file = paste0("agg-network.nodes.tsv"))

# load GO annotation
go_annotation <- read_tsv(
  cmd_line_args$options$go_annotation,
  col_names = c("GeneID", "TermID", "Domain"),
  show_col_types = TRUE
)
gene_id2term_id <- go_annotation |> 
  dplyr::select(GeneID, TermID) |> 
  as.matrix()

# get numbers of samples from each file
expts <- read_tsv(cmd_line_args$options$expts_file, show_col_types = FALSE, col_names = c("expt"))
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

# sample orderings of files
ordering <- vector("list", length = 10)
# order first by increasing sample size
ordering[[1]] <- purrr::map_chr(levels(sample_counts$expt),
                                \(x) grep(x, cmd_line_args$args, value = TRUE))
# next 8 randomly
ordering[2:9] <- purrr::map(1:8, \(x) sample(cmd_line_args$args))
# last one by decreasing sample size
ordering[[10]] <- rev(ordering[[1]])


cat("Starting..\n")
Sys.time()
mem_used()
purrr::map(seq_along(ordering),
           \(x) make_aggregate_network(x, ordering, annotation,
                                       gene_id2term_id, go_annotation))
Sys.time()
mem_used()

# make_aggregate_network(1, ordering, annotation, gene_id2term_id, go_annotation)
# mem_used()
# Sys.time()

# roc <- purrr::map(seq_along(roc_list),
#                   \(x) {
#                     tibble(
#                       ordering_num = x,
#                       num_networks = factor(c(seq_along(cmd_line_args$args), "agg", "agg-top-decile")),
#                       roc = roc_list[[x]]
#                     ) |> 
#                       mutate(num_networks = fct_inorder(num_networks))
#                   }) |> 
#   list_rbind()
# 
# ## plot step graph of AUC
# ggplot(roc, aes(x = num_networks, y = roc)) + 
#   stat_summary(fun.data = "mean_cl_boot") +
#   stat_summary(fun = mean, geom = "line", aes(group = 1)) +
#   theme_minimal() +
#   labs(x = "Number of aggregated networks", y = "Average AUROC")
# 
# # output auroc data
# write_tsv(roc, file = "AUROC-data-long.tsv")

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
