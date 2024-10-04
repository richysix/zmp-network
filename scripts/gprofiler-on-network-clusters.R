#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option("--min_cluster_size", action = "store", default = 20,
              help = "Minimum size of cluster to test [default %default]"),
  make_option("--output_file", action = "store",
              default = "clusters-go-enrichment.tsv",
              help = "Name of file to output the results to [default %default]"),
  make_option("--rds_file", action = "store", default = NULL,
              help = "Name of file to save the GO enrichment list to [default %default]"),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste(
  "Test members of each cluster in a network for enrichment of GO terms", 
  "Requires an input file with columns 'gene_id' and 'cluster_id'",
  sep = "\n"
)

if (interactive()) {
  cmd_args <- c("~/work/apocrita/data/scratch/bty114/zmp-network/nf/results/expt-zmp_ph192/all-tpm-t44.mci.I40.nodes.tsv") # add test files here
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc,
    usage = "Usage: %prog [options] nodes_file"
  ),
  args = cmd_args,
  positional_arguments = 1
)

packages <- c("tidyverse", "gprofiler2")
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

# run GO enrichment for a cluster
run_gost <- function(cluster_df, background_ids) {
  cluster_num <- cluster_df$cluster_id[1]
  go_results <- tryCatch(
    gost(
      query = cluster_df$gene_id,
      organism = 'drerio',
      significant = TRUE,
      correction_method = 'fdr',
      domain_scope = "custom",
      evcodes = TRUE,
      custom_bg = background_ids,
      sources = c("GO:BP", "GO:CC", "GO:MF"),
      highlight = TRUE
    ),
    message = function(cnd) {
      if (grepl("No results to show", conditionMessage(cnd))) {
        message(glue::glue("cluster-{cluster_num}: No results to show"))
        return(NULL)
      } else {
        message(cnd)
      }
    }
  )
  
  if (!is.null(go_results)) {
    go_results$query <- paste0("cluster-", cluster_num)
  }
  return(go_results)
}

# load nodes file
nodes <- read_tsv(file = cmd_line_args$args[1], show_col_types = FALSE)
# check for genes and cluster columns
if (any(!c("gene_id", "cluster_id") %in% colnames(nodes))) {
  stop(paste("Input file must have columns named",
             paste0(c("'gene_id'", "'cluster_id'"), collapse = " and ")))
}

# filter clusters by size
filtered_nodes <- nodes |>
  group_by(cluster_id) |>
  summarise(size = length(gene_id)) |>
  dplyr::filter(size >= cmd_line_args$options[['min_cluster_size']]) |>
  semi_join(nodes, y = _, by = join_by(cluster_id))

# run GO enrichment
go_enrichment_results <- split(filtered_nodes, filtered_nodes$cluster_id) |>
  map(\(x) run_gost(x, nodes$gene_id)) |>
  rlang::set_names(paste0('cluster-', unique(filtered_nodes$cluster_id))) |>
  keep(\(x) !is.null(x)) # remove ones with no enrichments

# bind results dfs together and write to output file
map(go_enrichment_results, "result") |>
  list_rbind(names_to = "query") |>
  write_tsv(file = cmd_line_args$options[['output_file']])

if (!is.null(cmd_line_args$options[['rds_file']])) {
  write_rds(go_enrichment_results, file = cmd_line_args$options[['rds_file']])
}

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
