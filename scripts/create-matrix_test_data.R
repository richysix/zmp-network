# create test data
setwd("/data/scratch/bty114/zmp-network")
library(tidyverse)
annotation <- read_tsv("nf/reference/danio_rerio-e99-annotation.txt",
                       col_names = c("GeneID", "Chr", "Start", "End", "Strand",
                                     "Biotype", "Name", "Description"))
go_annotation <- read_tsv(
  "nf/reference/danio_rerio_e109_go.txt", 
  col_names = c("GeneID", "TermID", "Domain"),
  show_col_types = TRUE
)
counts <- table(go_annotation$TermID)
go_term_counts <- tibble(TermID = names(counts), count = unname(counts))
genes_to_use <- go_annotation |>
  dplyr::summarise(count = n(), .by = TermID) |> 
  arrange(count) |>
  dplyr::filter(count > 20) |>
  head(50) |>
  dplyr::select(TermID) |>
  inner_join(go_annotation)

tab1 <- readr::read_tsv("nf/work/39/b9024100ffb71903bbfe639f900cc0/zmp_ph46-tpm-filtered-orig.tab",
                 col_names = c("idx", "GeneID"))
tab2 <- readr::read_tsv("nf/work/b9/4920fcd3277f76d082e3582e9851d2/zmp_ph250-tpm-filtered-orig.tab",
                 col_names = c("idx", "GeneID"))
# see how many are in the selected genes
# 94 keep these
both <- dplyr::inner_join(tab1, genes_to_use) |>
  dplyr::inner_join(tab2, by = join_by("GeneID" == "GeneID")) |>
  dplyr::select(idx.x, GeneID, idx.y) |>
  unique()

# get 25 genes that are in 1 only, 2 only, both and neither
only_tab1 <- anti_join(tab1, tab2, by = join_by("GeneID" == "GeneID")) |>
  anti_join(both, by = join_by("GeneID" == "GeneID")) |>
  head(50)
only_tab2 <- anti_join(tab2, tab1, by = join_by("GeneID" == "GeneID")) |> 
  head(50)
neither <- 
  anti_join(
    annotation,
    full_join(tab1, tab2, by = join_by("GeneID" == "GeneID"))
  ) |> head(39) |>
  dplyr::select(GeneID)
  
genes <- rbind(
  dplyr::select(only_tab1, GeneID),
  dplyr::select(only_tab2, GeneID),
  dplyr::select(both, GeneID),
  dplyr::select(neither, GeneID)
)

mat1 <- vroom::vroom("nf/work/39/b9024100ffb71903bbfe639f900cc0/zmp_ph46-tpm-filtered-orig.mat",
                     delim = "\t", col_names = TRUE) |> 
  as.matrix()

# get only_tab1 and both from matrix and output with correct col nums
tab1_selected <- rbind(
  only_tab1,
  dplyr::select(both, idx = idx.x, GeneID)
)

# filter for the genes from mat1 and write out
idx <- mat1[, 1]
mat1 <- mat1[, -1]
rownames(mat1) <- idx
mat1[ as.character(tab1_selected$idx), as.character(tab1_selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = tab1_selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", tab1_selected$idx)) |> 
  write_tsv("nf/test/zmp_ph46-tpm-filtered-orig.test.mat")

tab1_selected |> write.table(file = "nf/test/zmp_ph46-tpm-filtered-orig.test.tab",
                    quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
rm(mat1)

tab2_selected <- rbind(
  only_tab2,
  dplyr::select(both, idx = idx.y, GeneID)
)

mat2 <- vroom::vroom("nf/work/b9/4920fcd3277f76d082e3582e9851d2/zmp_ph250-tpm-filtered-orig.mat",
                     delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat2[, 1]
mat2 <- mat2[, -1]
rownames(mat2) <- idx
mat2[ as.character(tab2_selected$idx), as.character(tab2_selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = tab2_selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", tab2_selected$idx)) |> 
  write_tsv("nf/test/zmp_ph250-tpm-filtered-orig.test.mat")
tab2_selected |> 
  write.table(file = "nf/test/zmp_ph250-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
rm(mat2)

rbind(
  dplyr::select(only_tab1, GeneID),
  dplyr::select(only_tab2, GeneID),
  dplyr::select(both, GeneID),
  dplyr::select(neither, GeneID)
) |> inner_join(annotation) |> 
  write.table(file = "nf/reference/danio_rerio-e99-annotation.test.txt",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# output go annotation from just those genes
rbind(
  dplyr::select(only_tab1, GeneID),
  dplyr::select(only_tab2, GeneID),
  dplyr::select(both, GeneID),
  dplyr::select(neither, GeneID)
) |> inner_join(go_annotation) |> 
  write.table(file = "nf/reference/danio_rerio_e109_go.test.txt",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# zmp_ph204
tab <- readr::read_tsv("nf/work/52/28255842fe1510aac315f1b048b6b3/zmp_ph204-tpm-filtered.tab",
                        col_names = c("idx", "GeneID"))
selected <- inner_join(tab, genes)
mat <- vroom::vroom("nf/work/52/28255842fe1510aac315f1b048b6b3/zmp_ph204-tpm-filtered-orig.mat",
                     delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat[, 1]
mat <- mat[, -1]
rownames(mat) <- idx
mat[ as.character(selected$idx), as.character(selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", selected$idx)) |> 
  write_tsv("nf/test/zmp_ph204-tpm-filtered-orig.test.mat")
selected |> 
  write.table(file = "nf/test/zmp_ph204-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# zmp_ph192
tab <- readr::read_tsv("nf/work/eb/baf2405e6243d10c45977b1bc1df5c/zmp_ph192-tpm-filtered.tab",
                       col_names = c("idx", "GeneID"))
selected <- inner_join(tab, genes)
mat <- vroom::vroom("nf/work/eb/baf2405e6243d10c45977b1bc1df5c/zmp_ph192-tpm-filtered-orig.mat",
                    delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat[, 1]
mat <- mat[, -1]
rownames(mat) <- idx
mat[ as.character(selected$idx), as.character(selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", selected$idx)) |> 
  write_tsv("nf/test/zmp_ph192-tpm-filtered-orig.test.mat")
selected |> 
  write.table(file = "nf/test/zmp_ph192-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# zmp_ph238
tab <- readr::read_tsv("nf/work/d3/a14340e7425aed21222df3534f04c8/zmp_ph238-tpm-filtered.tab",
                       col_names = c("idx", "GeneID"))
selected <- inner_join(tab, genes)
mat <- vroom::vroom("nf/work/d3/a14340e7425aed21222df3534f04c8/zmp_ph238-tpm-filtered-orig.mat",
                    delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat[, 1]
mat <- mat[, -1]
rownames(mat) <- idx
mat[ as.character(selected$idx), as.character(selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", selected$idx)) |> 
  write_tsv("nf/test/zmp_ph238-tpm-filtered-orig.test.mat")
selected |> 
  write.table(file = "nf/test/zmp_ph238-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# zmp_ph213
tab <- readr::read_tsv("nf/work/4e/5a354574bcacb6e4d4adb711f0c376/zmp_ph213-tpm-filtered.tab",
                       col_names = c("idx", "GeneID"))
selected <- inner_join(tab, genes)
mat <- vroom::vroom("nf/work/4e/5a354574bcacb6e4d4adb711f0c376/zmp_ph213-tpm-filtered-orig.mat",
                    delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat[, 1]
mat <- mat[, -1]
rownames(mat) <- idx
mat[ as.character(selected$idx), as.character(selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", selected$idx)) |> 
  write_tsv("nf/test/zmp_ph213-tpm-filtered-orig.test.mat")
selected |> 
  write.table(file = "nf/test/zmp_ph213-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

# zmp_ph71
tab <- readr::read_tsv("nf/work/b8/8c34182ebc4d6a86f21b8bfc901a3c/zmp_ph71-tpm-filtered.tab",
                       col_names = c("idx", "GeneID"))
selected <- inner_join(tab, genes)
mat <- vroom::vroom("nf/work/b8/8c34182ebc4d6a86f21b8bfc901a3c/zmp_ph71-tpm-filtered-orig.mat",
                    delim = "\t", col_names = TRUE) |> 
  as.matrix()
idx <- mat[, 1]
mat <- mat[, -1]
rownames(mat) <- idx
mat[ as.character(selected$idx), as.character(selected$idx) ] |> 
  as_tibble() |> 
  mutate(idx = selected$idx) |> 
  relocate(idx, ) |> 
  magrittr::set_colnames(c("dummy", selected$idx)) |> 
  write_tsv("nf/test/zmp_ph71-tpm-filtered-orig.test.mat")
selected |> 
  write.table(file = "nf/test/zmp_ph71-tpm-filtered-orig.test.tab",
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
