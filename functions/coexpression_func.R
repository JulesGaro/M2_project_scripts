# Functions use to compute co-expression matrices at both levels

cell_level_coexpression <- function(expr_mat,meta, min_n_cell,batch_size = 8) {
  # Compute coexpression across cells within each cell pops 
  # 
  # Args:
  #   expr_mat : expression matrix 
  #   meta : tibble with at least cols [cell, cell_type, cell, cell_pop]
  #   min_n_cell : minimum cell pop size
  # 
  # Return:
  #     aggre_coex_mat : the rank aggregated coexpression matrix 
  
  # meta must have columns: cell, cell_type, cell, cell_pop
  pop_dims <- meta %>% dplyr::select(cell_pop, cell_type, cell) %>% unique()
  
  # Filter out cell pops (cell from a subject within a cell type) that have less than
  # 'min_n_cell' cell
  cell_pops <- meta %>% group_by(cell_pop) %>% summarize(n_cel = n())
  cell_pops <- cell_pops %>% left_join(pop_dims, by = "cell_pop")
  cell_pops <- cell_pops %>% mutate(valid = n_cel > min_n_cell)
  valid_cell_pops <- cell_pops %>% filter(valid) %>% as.tibble()
  
  pops <- valid_cell_pops$cell_pop %>% unique() %>% sort()
  # Get expression matrices of each cell pops
  expr_mats <- pops %>% lapply(function(cp) {
    meta_cp <- meta %>% filter(cell_pop == cp)
    expr_mat[,meta_cp$cell]
  })
  names(expr_mats) <- pops

  # Compute coexpression for each cell pop 
  # coexpression of cell pop within a cell types or computed then aggregate (sum)
  # by batch that are computed with paralelization to be more efficient on time
  # and memory
  vec <- c(1:length(pops))
  pop_batchs_num <- split(vec,ceiling(seq_along(vec) / batch_size))
  
  # too avoid any bug if the last batch contains only 1 matrix it is added to the
  # second last batch
  if (length(pop_batchs_num[[length(pop_batchs_num)]]) == 1) {
    s_last <- pop_batchs_num[[(length(pop_batchs_num)-1)]]
    pop_batchs_num[[length(pop_batchs_num)]] <- append(pop_batchs_num[[length(pop_batchs_num)]], s_last[batch_size])
    pop_batchs_num[[(length(pop_batchs_num)-1)]] <- s_last[1:(batch_size-1)]
  }
  
  # Computing intermediate rank aggregated coexpression matrices of all batches
  inter_coex_mats <- pop_batchs_num %>% lapply(function(vec) {
    print("compute batch...")
    compute_inter_coex_mat(vec,pops,expr_mats,batch_size)
  })
  # aggregating the intermediate matrices (mean)
  aggre_coex_mat <- inter_coex_mats %>% aggregate_coex_batchs(n_pop = length(inter_coex_mats))
  
  return(aggre_coex_mat)
}

compute_inter_coex_mat <- function(vec,pops,expr_mats,batch_size) {
  # Compute rank aggregated matrix (sum) for a batch of expression matrices
  # 
  # Args:
  #   vec : vector of position indicating which pops are in the batch
  #   pops : vector of cell pops names
  #   expr_mats : the list of the cell pops expression matrices
  # 
  # Return:
  #   aggregate_mat: rank aggregated matrix of the batch (sum)
  
  # Compute correlation for the batch
  print("Compute coexpression for :")
  print(pops[vec])
  coex_mats <- expr_mats[pops[vec]] %>% mclapply(function(expr_mat) {
    coex_mat <- compute_coex_mat(expr_mat)
  }, mc.cores = batch_size)
  
  
  print("Remove NA and rank gene-gene correlation values...")
  genes <- rownames(coex_mats[[1]])
  # Replace NA by 0, rank values with ties.method = random, return as vectors
  rank_matrix_values <- coex_mats %>% mclapply(function(coex_mat) {
    coex_mat@x[is.na(coex_mat@x)] <- 0 
    set.seed(42)
    return(rank(coex_mat,ties.method = "random"))
  }, mc.cores = batch_size)
  rm(coex_mats)
  gc(verbose = F)
  
  print("Aggregate and reconstruct matrix...")
  # Make the sum of the vector and then rebuild matrices
  aggregate_mat <- Reduce("+", rank_matrix_values) %>% matrix(nrow=length(genes),
                                                              ncol=length(genes)) %>% as("dgCMatrix")
  rownames(aggregate_mat) <- genes
  colnames(aggregate_mat) <- genes
  return(aggregate_mat)
}

aggreg_subj_in_ct <- function(expr_mat, meta) {
  # Aggregate cell expression vectors by subject in an cell type expression matrix
  # 
  # Args:
  #   expr_mat : expression matrix
  #   meta : tibble with at least the cols [cell, subject, cell_type]
  #   
  # Return:
  #   aggre_subj_mat : genes by subjects matrix
  
  expr_mat <- expr_mat[, meta$cell]
  
  # Aggregate cells by subject (get the mean of expression)
  subjectGrps <- meta %>% 
    dplyr::select(cell, subject) %>% 
    mutate(value = 1) %>% 
    spread(subject, value) %>% 
    as.data.frame()
  
  rownames(subjectGrps) <- subjectGrps[, 1]
  subjectGrps <- subjectGrps[, -1]
  subjectGrps <- subjectGrps %>% as.matrix()
  subjectGrps[is.na(subjectGrps)] <- 0
  
  expr_mat <- expr_mat[, rownames(subjectGrps)]
  aggre_subj_mat <- expr_mat %*% subjectGrps
  aggre_subj_mat <- aggre_subj_mat %>% as.matrix()
  
  # compute the average
  nCels <- subjectGrps %>% apply(2, sum)
  aggre_subj_mat <- t( t(aggre_subj_mat) / nCels) %>% as("dgCMatrix")
  
  return(aggre_subj_mat)
}

compute_coex_mat <- function(expr_mat) {
  # Compute the coexpression matrix for the given expression matrix
  # 
  # Args:
  #   expr_mat : the expression matrix
  # 
  # Return:
  #   coex_mat : the coexpression matrix
  
  # Sort the cols and rows
  expr_mat <- expr_mat[sort(rownames(expr_mat)), sort(colnames(expr_mat))]
  
  # Transpose the matrix to get correlation on genes
  expr_mat <- Matrix::t(expr_mat)
  
  # Compute pearson correlation on the matrix
  coex_mat <- WGCNA::cor(expr_mat, method = "pearson") %>% as("dgCMatrix")
  
  return(coex_mat)
}