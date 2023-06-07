#!/usr/bin/env Rscript

source("../functions/data_process_func.R")
source("../functions/coexpression_func.R")
library(dplyr)
library(stringr)
library(tibble)
library(Matrix)
library(parallel)
library(tidyr)
setwd("~")
opt <- list(
    GENE_THR = 0.05,
    SAMPLE_THR = 2,
    PERCENTILE_THR = 0.99,
    MIN_N_CELL = 20
)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1] # path should end by a "/"
dataset <- args[2]
print("computing :")
print(dataset)
print("Load Preprocessed data...")

load(paste(path, dataset, "_preprocess.RData", sep = ""))
data <- get(paste(dataset, "_output", sep = ""))
rm(list = c(paste(dataset, "_output", sep = "")))
gc(verbose = F)

cell_types <- data$meta$cell_type %>%
    unique() %>%
    sort()

# Separate ctmat and meta by cell types
ctmats <- cell_types %>% mclapply(function(ct) {
    filtered_meta <- data$meta %>% filter(cell_type == ct)
    return(data$ctmat[, filtered_meta$cell])
}, mc.cores = 6)
names(ctmats) <- cell_types

metas <- cell_types %>% mclapply(function(ct) {
    return(data$meta %>% filter(cell_type == ct))
}, mc.cores = 6)
names(metas) <- cell_types

# Data cleaning
print("Data cleaning...")
ctmats <- cell_types %>% mclapply(function(ct) {
    return(clean_ctmat(ctmat = ctmats[[ct]], gene_thr = opt$GENE_THR, sample_thr = opt$SAMPLE_THR))
}, mc.cores = 6)
names(ctmats) <- cell_types
# Get rid of the cells in metadata as well
metas <- cell_types %>% mclapply(function(ct) {
    return(metas[[ct]][metas[[ct]]$cell %in% colnames(ctmats[[ct]]), ])
}, mc.cores = 6)
names(metas) <- cell_types

print("Data normalization...")
# Data Normalization
ctmats <- cell_types %>% mclapply(function(ct) {
    return(normalize_ctmat(ctmats[[ct]]))
}, mc.cores = 6)
names(ctmats) <- cell_types

# ==========COMPUTE CO-EXPRESSION FOR 1ST DATASET=================================
print("Compute cell level Co-expression...")

print("    Compute co-expression")
# Divide expr_mat by cell type to speed up process and use les memory
cell_coex_mats <- cell_types %>% lapply(function(ct) {
    cell_coex_mat <- cell_level_coexpression(
        expr_mat = ctmats[[ct]],
        meta = metas[[ct]],
        min_n_cell = opt$MIN_N_CELL,
        batch_size = 10
    )
    gc(verbose = F)
    return(cell_coex_mat)
})
names(cell_coex_mats) <- cell_types

print("Compute subject level Co-expression...")
print("    Aggregate subjects")
aggre_subj_mats <- cell_types %>% mclapply(function(ct) {
    return(aggreg_subj_in_ct(expr_mat = ctmats[[ct]], meta = metas[[ct]]))
}, mc.cores = 6)
names(aggre_subj_mats) <- cell_types

rm(ctmats)
gc(verbose = F)

print("    Compute co-expression")
subj_coex_mats <- cell_types %>% mclapply(function(ct) {
    compute_coex_mat(expr_mat = aggre_subj_mats[[ct]])
}, mc.cores = 6)
names(subj_coex_mats) <- cell_types

rm(aggre_subj_mats)
gc(verbose = F)

print("Saving results...")
output <- list(cell_coex_mats, subj_coex_mats, metas, dataset, cell_types, opt)
names(output) <- c("cell_coex_mats", "subj_coex_mats", "metas", dataset, "cell_types", "opt")
save(output, file = paste(path, dataset, "_coex.RData", sep = ""))
