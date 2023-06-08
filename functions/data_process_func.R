# Functions used to clean and normalize datasets before performing
# co-expression analysis
library(dplyr)
library(Matrix)

clean_ctmat <- function(ctmat, gene_thr = 0.1, sample_thr = 2) {
  # Filter out Genes and Cells from ctmat that do not pass given treshold
  #
  # Args:
  #   ctmat : expression matrix
  #   gene_thr : % of cells in which a gene need to be expressed to be kept
  #   sample_thr : standard deviation of n_genes beyond which cells are filtered out
  #
  # Return:
  #   ctmat : cleaned expression matrix

  # Filter out Genes
  genes <- rownames(ctmat)
  genes <- genes[
    rowSums(ctmat > 0) %>% sapply(function(row_count) row_count > ncol(ctmat) * gene_thr) == T
  ]

  ctmat <- ctmat[rownames(ctmat) %in% genes, ]
  gc()

  # Filter out Cells
  valid_cells <- tibble(cell = colnames(ctmat), n_gene = colSums(ctmat > 0))

  n_genes <- list(m = mean(valid_cells$n_gene), sd = sd(valid_cells$n_gene))

  valid_cells <- valid_cells %>%
    filter(
      n_gene > (n_genes$m - (n_genes$sd * sample_thr)),
      n_gene < (n_genes$m + (n_genes$sd * sample_thr))
    )

  ctmat <- ctmat[, valid_cells$cell]

  return(ctmat)
}

normalize_ctmat <- function(ctmat) {
  # Convert the count expression matrix into a cpm matrix
  #
  # Args:
  #   ctmat : count expression matrix
  #
  # Return:
  #   ctmat : count per million expression matrix
  print("Calc libsize...")
  lib_sizes <- ctmat %>% colSums()

  print("Getting Intervals...")
  # Retrieve columns as intervals in non zero values vector
  intervals <- c(1:(length(ctmat@p) - 1)) %>% lapply(function(i) {
    c(ctmat@p[i] + 1, ctmat@p[i + 1])
  })

  print("Normalize Values...")
  # Normalized the values
  ctmat@x <- c(1:length(lib_sizes)) %>%
    lapply(function(i) {
      ctmat@x[unlist(intervals[i])[1]:unlist(intervals[i])[2]] / lib_sizes[i] * 1e6
    }) %>%
    unlist()

  return(ctmat)
}
