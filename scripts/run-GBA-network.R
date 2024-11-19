#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--auc_file", action = "store", default = 'GBA-auc.tsv',
              help = "Name of output file for AUC values of annotations [default %default]"),
  make_option("--scores_file", action = "store", default = 'GBA-gene-scores.tsv',
              help = "Name of output file for gene scores [default %default]"),
  make_option("--plots_file", action = "store", default = 'GBA-plots.pdf',
              help = "Name of plots file [default %default]"),
  make_option("--min.term.size", action = "store", default = 10,
              help = "Minimum number of annotations to each group [default %default]"),
  make_option("--max.term.size", action = "store", default = 1000,
              help = "Maximum number of annotations to each group [default %default]"),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste('Script to run Guilt By Association scoring of annotation terms', 
              'to assess quality of network at clustering similar functions', 
              'Expects a nodes files, an edges file and an annotation file', 
              sep = "\n")
if (interactive()) {
  cmd_args <- c('nodes.csv', 'edges.csv', 'danio_rerio_e85_go.txt')
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
}
cmd_line_args <- parse_args(
  OptionParser(
    option_list = option_list,
    description = desc, prog = 'run-GBA-network.R',
    usage = "Usage: %prog [options] nodes_file edges_file annotation_file" ),
  args = cmd_args,
  positional_arguments = 3
)

packages <- c('EGAD', 'tidyverse')
for (package in packages) {
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

show_cols <- cmd_line_args$options$debug

# load data
# nodes
nodes <- read_tsv(cmd_line_args$args[1], show_col_types = show_cols)
if (nrow(nodes) == 0) {
  message("The graph has no nodes. Skipping GBA analysis")
  file.create(
    cmd_line_args$options$auc_file,
    cmd_line_args$options$scores_file,
    cmd_line_args$options$plots_file
  )
  quit(save = "no", status = 0)
}

# nodes <- filter(nodes, cluster_id >= 80, cluster_id <= 120)
node_idx2gene_id <- dplyr::select(nodes, node_idx, gene_id)
# edges
edges <- read_tsv(cmd_line_args$args[2], show_col_types = show_cols)
# edges <- filter(edges, source %in% nodes$node_idx, target %in% nodes$node_idx)

# convert edges to gene ids
edges_by_gene_id <- left_join(edges, node_idx2gene_id, 
                              by = c("source" = "node_idx"), 
                              relationship = "many-to-one") |> 
  rename(Gene1 = gene_id) |> 
  left_join(node_idx2gene_id, by = c("target" = "node_idx"), 
            relationship = "many-to-one") |> 
  dplyr::select(Gene1, Gene2 = gene_id) |> 
  as.matrix()

# create network as adj matrix
network <- make_gene_network(edges_by_gene_id, nodes$gene_id)

# network topology properties
assort <- assortativity(network)
nd <- node_degree(network)

# annotation
annotation <- read_tsv(
  cmd_line_args$args[3], 
  col_names = c("GeneID", "TermID", "Domain"),
  show_col_types = show_cols
)
gene_id2term_id <- annotation |> 
  dplyr::select(GeneID, TermID) |> 
  as.matrix()
annotations <- make_annotations(gene_id2term_id, nodes$gene_id, annotation$TermID)
annotations <- filter_network_cols(
  annotations,
  min = cmd_line_args$options$min.term.size,
  max = cmd_line_args$options$max.term.size
)

# run GBA
Anno_groups_voted <- run_GBA(
  network,
  annotations,
  min = cmd_line_args$options$min.term.size,
  max = cmd_line_args$options$max.term.size
)

auc_Anno_nv <- Anno_groups_voted[[1]][,1]
Anno_multifunc_assessment <- calculate_multifunc(annotations)
# For genes
optimallist_genes <-  Anno_multifunc_assessment[,4] 
auc_Anno_mf <- auc_multifunc(annotations, optimallist_genes)
# Or for GO groups
Anno_genes_multifunc_assessment <- calculate_multifunc(t(annotations))
optimallist_GO <-  Anno_genes_multifunc_assessment[,4] 
auc_gene_mf <- auc_multifunc(t(annotations), optimallist_GO)

# Output results
as_tibble(Anno_groups_voted[[1]], rownames = "TermID") |>
  arrange(desc(auc)) |>
  write_tsv(file = cmd_line_args$options$auc_file)

as_tibble(Anno_groups_voted[[2]], rownames = "GeneID") |>
  arrange(GeneID) |>
  write_tsv(cmd_line_args$options$scores_file)

pdf(file = cmd_line_args$options$plots_file)
# plot node degree distribution
stack(nd) |> 
  rename(degree = values, Gene = ind) |> 
  ggplot(aes(x = degree)) +
  geom_density()
# plot ROC curve
scores <- neighbor_voting(annotations, network, nFold = 3, output = "AUROC", 
                          FLAG_DRAW = TRUE) 
# plot AUROC distribution
plot_distribution(auc_Anno_nv, xlab = "Neighbor voting AUROC ",
                  ylab = "Number of functional terms",
                  b = 30, col = "gray64",
                  density = FALSE, avg = FALSE, bars = TRUE) |> 
  invisible()
# plot optimal node degree distribution
plot_distribution(auc_Anno_mf, xlab = "Optimal GO Ranking AUROC",
                  ylab = "Number of functional terms",
                  b = 20, col = "gray64", 
                  density = FALSE, avg = FALSE, bars = TRUE) |> 
  invisible()

plot_distribution(auc_gene_mf, xlab = "Optimal Gene Ranking AUROC",
                  ylab = "Number of functional terms",
                  b = 20, xlim = c(0.2,1), ylim = c(0,4400), col = "gray64", 
                  density = FALSE, avg = FALSE, bars = TRUE) |> 
  invisible()

dev.off() |> 
  invisible()

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
