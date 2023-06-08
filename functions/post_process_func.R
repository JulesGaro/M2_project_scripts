library(tibble)
library(parallel)
library(dplyr)
library(Matrix)
library(data.table)

compute_preservation <- function(mat1, mat2) {
  # compute the preservation between 2 binary matrices
  #
  # Args:
  #   mat1/mat2 : A binary matrix, can be sparse
  #
  # Return:
  #   tp_frac : The fraction of common 1 binary between the 2 matrices


  pos1 <- which(mat1[upper.tri(mat1, diag = F)] == 1)
  pos2 <- which(mat2[upper.tri(mat2, diag = F)] == 1)

  return(length(intersect(pos1, pos2)) / length(pos1))
}

get_upper_tri_tbl <- function(mat) {
  # Return a tibble of the pairs of genes name (without doublon due to symmetry)
  #
  # Args:
  #   mat: a genes by genes symmetric matrix
  #
  # Return:
  #   a tibble with all the genes pairs name like 'gene1_gene2'
  #
  up_tri <- upper.tri(mat, diag = FALSE)
  genes_pair <- which(up_tri == T) %>%
    mclapply(function(i) {
      paste(rownames(mat)[i %/% ncol(up_tri) + 1], rownames(mat)[i %% nrow(up_tri)], sep = "_")
    }, mc.cores = 6) %>%
    unlist()
  gc()
  col_row_tibble <- tibble(genes_pair = genes_pair)
  return(col_row_tibble)
}

auroc <- function(ranking, links) {
  # Compute the AUROC for the given ranked genes_pair and list of true positive links
  #
  # Args:
  #   ranking: pair of genes ranked by correlation (descending) for network 1
  #   links: list of true positive pair of genes (top 1% of correlation) for network 2
  #
  # Return:
  #   (vector):
  #     n_true : number of true positives
  #     n_full : number of ranked genes pair
  #     auroc : Fraction of the area under the ROC curve
  #     p_value : p_value of the mann-whitney test
  #
  rank_len <- length(ranking)
  links_len <- length(links)
  rank <- c(1:rank_len)
  in_set <- c(ranking %chin% links)
  setStatus <- tibble(item = ranking, rank = rank, in_set = in_set)
  rm(ranking)
  rm(links)
  gc()
  result <- wilcox.test(rank ~ in_set, data = setStatus)
  ns <- setStatus %>%
    group_by(in_set) %>%
    summarize(n = n()) %>%
    arrange(in_set)
  rm(setStatus)
  gc()

  # auc
  auc <- result$statistic / (as.double(ns$n[1]) * as.double(ns$n[2]))
  return(c(n_true = links_len, n_full = rank_len, auroc = auc, pvalue = result$p.value))
}
