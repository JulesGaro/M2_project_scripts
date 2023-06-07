#!/usr/bin/env Rscript

source("../functions/coexpression_func.R")
library(dplyr)
library(stringr)
library(tibble)
library(Matrix)
library(parallel)
library(tidyr)
setwd("~/")
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

load(paste(path, dataset, "_coex.RData", sep = ""))
data <- output
rm(output)
gc(verbose = F)

cell_types <- data$cell_types %>% sort()

# ==========COMPUTE PERCENTILE AND LINKS=========================

print("    Compute percentile matrix")
cell_perc_mats <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data$cell_coex_mats[[ct]])
}, mc.cores = 6)
names(cell_perc_mats) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
cell_lnks_mats <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = cell_perc_mats[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(cell_lnks_mats) <- cell_types

rm(cell_perc_mats)
gc(verbose = F)

print("    Compute percentile matrix")
subj_perc_mats <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data$subj_coex_mats[[ct]])
}, mc.cores = 6)
names(subj_perc_mats) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
subj_lnks_mats <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = subj_perc_mats[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(subj_lnks_mats) <- cell_types
rm(subj_perc_mats)
gc(verbose = F)

output <- list(subj_lnks_mats, cell_lnks_mats, data$metas, dataset, cell_types, opt)
names(output) <- c("subj_lnks_mats", "cell_lnks_mats", "metas", "name", "cell_types", "opt")
save(output, file = paste(path, dataset, "_single_process.RData", sep = ""))
