# Function used to compute non-weighted network from coexpression matrices 

compute_perc_mat <- function(coex_mat) {
  # Compute the percentile matrix for the given coexpression matrix
  # 
  # Args:
  #   coex_mat : the coexpression matrix
  # 
  # Return:
  #   perc_mat : the percentile matrix
  
  coex_mat_copy <- coex_mat
  # Set lower triangle to NA so we don't rank pairs twice
  coex_mat_copy[lower.tri(coex_mat_copy,diag = T)] <- NA
  coexVec <- coex_mat_copy %>% as.vector()
  perc_mat <- coexVec %>% percent_rank() %>% matrix(ncol = ncol(coex_mat))
  rownames(perc_mat) <- rownames(coex_mat)
  colnames(perc_mat) <- colnames(coex_mat)
  # Lower triangle is then set to 0, diagonal to 1
  diag(perc_mat) <- 1
  perc_mat[lower.tri(perc_mat,diag = F)] <- 0
  perc_mat <- perc_mat %>% as("dgCMatrix")
  gc()
  
  return(perc_mat)
}

compute_lnks_mat <- function(perc_mat, thresh = 0.99) {
  # Compute the links matrix for the given percentile matrix
  # 
  # Args:
  #   perc_mat : the percentile matrix
  #   thresh : percentile threshold 
  #
  # Return:
  #   lnks_mat : the links matrix
  lnks_mat <- perc_mat
  
  # Value higher than the treshold (top 0.1% best correlation rank) are links
  lnks_mat@x <- replace(lnks_mat@x, lnks_mat@x >= thresh, 1)
  
  lnks_mat@x <- replace(lnks_mat@x, lnks_mat@x < thresh, 0)
  
  # Going back and forth so the 0 values or getting "pushed out" of the values vector
  lnks_mat <- lnks_mat %>% as.matrix() %>% as("dgCMatrix")
  gc()
  return(lnks_mat)
}