#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE,
              help = "Print extra output [default %default]")
)

desc <- paste("Script to plot the effect of changing the threshold for",
              "filtering spurious correlations", sep = "\n")

if (interactive()) {
  cmd_args <- c() # add test files here
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
  positional_arguments = 0
)

packages <- c("tidyverse", "ggrepel", "ggtext")
for (package in packages) {
  library(package, character.only = TRUE) |>
    suppressWarnings() |>
    suppressPackageStartupMessages()
}

expt_names <- paste0("zmp_ph", c("46", "71", "204", "213", "192", "238", "250"))
cor_data <- purrr::map(expt_names,
                       \(x) read_tsv(paste0("spurious-correlations/cor-by-num-zeros-hist-", x, ".tsv")) |> 
                              filter(count > 25, num_zeros > 10) |> 
                              mutate(expt = x)
                         )

calc_genes_removed <- function(df, threshold){
  tibble(
    expt = df$expt[1],
    threshold = threshold,
    genes_removed = filter(df, count > threshold) |>
      pull(num_genes) |> sum()
  )
}
calc_genes_removed_by_threshold <- function(df, thresholds) {
  purrr::map(thresholds, \(x) calc_genes_removed(df, x)) |> 
    list_rbind()
}

genes_removed_by_threshold <- purrr::map(cor_data,
            \(x) calc_genes_removed_by_threshold(x, seq(20,100,5))) |> 
  list_rbind() |> 
  mutate(expt = factor(expt, levels = expt_names))

label_df <- genes_removed_by_threshold |> 
  filter(threshold == 100) |> 
  mutate(colour = biovisr::cbf_palette(length(expt)),
         label_text = glue::glue("<span style='color:{colour}'>{expt}</span>"))

svglite::svglite(filename = file.path("plots", "filtering-spurious-correlations.svg"),
                 width = 10, height = 6)  
ggplot(genes_removed_by_threshold,
       aes(x = threshold, y = genes_removed, colour = expt)) +
  geom_vline(xintercept = 25, linetype = "dashed") +
  geom_step(linewidth = 1.5,
            direction = "mid",
            show.legend = FALSE) +
  geom_text_repel(data = label_df,
            aes(label = expt),
            hjust = 0, nudge_x = 2,
            show.legend = FALSE,
            min.segment.length = 5) +
  annotate(geom = "text", x = 35, y = 4800, label = "threshold", hjust = 0) +
  annotate(geom = "curve", x = 35, xend = 25, y = 4800, yend = 4800,
           curvature = 0.2,
           arrow = arrow(length = unit(2, "mm"))) +
  scale_colour_manual(values = label_df$colour) +
  scale_x_continuous(limits = c(20,110),
                     breaks = seq(20,100,20)) +
  scale_y_continuous(limits = c(500,4800)) +
  labs(x = "Threshold (number of spurious correlations)",
       y = "Number of genes removed",
       title = "Having a stringent threshold removes less than <b>20%</b> of all genes",
       subtitle = "(<b>Not</b> including removing genes where all values are zero)") +
  theme_minimal() +
  theme(title = element_markdown())
dev.off()

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
