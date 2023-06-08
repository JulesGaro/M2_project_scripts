#!/usr/bin/env Rscript

source("../functions/processing_func.R")
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
dataset_1 <- args[2]
dataset_2 <- args[3]

print("computing :")
print(dataset_1)
print(dataset_2)
print("Load Preprocessed data...")

load(paste(path, dataset_1, "_single_process.RData", sep = ""))
data_1 <- output
rm(output)
load(paste(path, dataset_2, "_single_process.RData", sep = ""))
data_2 <- output
rm(output)
gc(verbose = F)

cell_types <- intersect(data_1$cell_types, data_2$cell_types) %>% sort()
print("Intersect common genes...")

# Intersect common genes
data_1$cell_coex_mats <- cell_types %>% mclapply(function(ct) {
    genes_1 <- rownames(data_1$cell_coex_mats[[ct]])
    genes_2 <- rownames(data_2$cell_coex_mats[[ct]])
    common_genes <- intersect(genes_1, genes_2) %>% sort()
    data_1$cell_coex_mats[[ct]][common_genes, common_genes]
}, mc.cores = 6)
names(data_1$cell_coex_mats) <- cell_types

data_2$cell_coex_mats <- cell_types %>% mclapply(function(ct) {
    genes_1 <- rownames(data_1$cell_coex_mats[[ct]])
    genes_2 <- rownames(data_2$cell_coex_mats[[ct]])
    common_genes <- intersect(genes_1, genes_2) %>% sort()
    data_2$cell_coex_mats[[ct]][common_genes, common_genes]
}, mc.cores = 6)
names(data_2$cell_coex_mats) <- cell_types

data_1$subj_coex_mats <- cell_types %>% mclapply(function(ct) {
    genes_1 <- rownames(data_1$cell_coex_mats[[ct]])
    genes_2 <- rownames(data_2$cell_coex_mats[[ct]])
    common_genes <- intersect(genes_1, genes_2) %>% sort()
    data_1$subj_coex_mats[[ct]][common_genes, common_genes]
}, mc.cores = 6)
names(data_1$subj_coex_mats) <- cell_types

data_2$subj_coex_mats <- cell_types %>% mclapply(function(ct) {
    genes_1 <- rownames(data_1$cell_coex_mats[[ct]])
    genes_2 <- rownames(data_2$cell_coex_mats[[ct]])
    common_genes <- intersect(genes_1, genes_2) %>% sort()
    data_2$subj_coex_mats[[ct]][common_genes, common_genes]
}, mc.cores = 6)
names(data_2$subj_coex_mats) <- cell_types
gc(verbose = F)

# ==========COMPUTE PERCENTILE AND LINKS FOR 1ST DATASET=========================

print("    Compute percentile matrix")
cell_perc_mats_1 <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data_1$cell_coex_mats[[ct]])
}, mc.cores = 6)
names(cell_perc_mats_1) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
cell_lnks_mats_1 <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = cell_perc_mats_1[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(cell_lnks_mats_1) <- cell_types

rm(cell_perc_mats_1)
gc(verbose = F)

print("    Compute percentile matrix")
subj_perc_mats_1 <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data_1$subj_coex_mats[[ct]])
}, mc.cores = 6)
names(subj_perc_mats_1) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
subj_lnks_mats_1 <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = subj_perc_mats_1[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(subj_lnks_mats_1) <- cell_types
rm(subj_perc_mats_1)
gc(verbose = F)

print("Saving 1st dataset...")
output_1 <- list(subj_lnks_mats_1, cell_lnks_mats_1, data_1$meta, dataset_1, cell_types, opt)
names(output_1) <- c("subj_lnks_mats", "cell_lnks_mats", "meta", "name", "cell_types", "opt")
save(output_1, file = paste(paste, dataset_1, "_", dataset_2, ".RData", sep = ""))

# ==========COOMPUTE PERCENTILE AND LINKS FOR 2ND DATASET================================================
print("Compute cell level Co-expression for the 2nd dataset...")

print("    Compute percentile matrix")
cell_perc_mats_2 <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data_2$cell_coex_mats[[ct]])
}, mc.cores = 6)
names(cell_perc_mats_2) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
cell_lnks_mats_2 <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = cell_perc_mats_2[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(cell_lnks_mats_2) <- cell_types

rm(cell_perc_mats_2)
gc(verbose = F)

print("    Compute percentile matrix")
subj_perc_mats_2 <- cell_types %>% mclapply(function(ct) {
    compute_perc_mat(coex_mat = data_2$subj_coex_mats[[ct]])
}, mc.cores = 6)
names(subj_perc_mats_2) <- cell_types
gc(verbose = F)

print("    Compute links matrix")
subj_lnks_mats_2 <- cell_types %>% mclapply(function(ct) {
    compute_lnks_mat(perc_mat = subj_perc_mats_2[[ct]], thresh = opt$PERCENTILE_THR)
}, mc.cores = 6)
names(subj_lnks_mats_2) <- cell_types

rm(subj_perc_mats_2)
gc(verbose = F)

print("Saving 2nd dataset...")

output_2 <- list(subj_lnks_mats_2, cell_lnks_mats_2, data_2$meta, dataset_2, cell_types, opt)
names(output_2) <- c("subj_lnks_mats", "cell_lnks_mats", "meta", "name", "cell_types", "opt")
save(output_2, file = paste(paste, dataset_2, "_", dataset_1, ".RData", sep = ""))
