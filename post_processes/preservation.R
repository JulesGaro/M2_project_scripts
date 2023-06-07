#!/usr/bin/env Rscript

source("../functions/post_process_func.R")
library(tibble)
library(parallel)
library(dplyr)
library(Matrix)
library(data.table)

setwd("~")

args <- commandArgs(trailingOnly = TRUE)
path <- args[1] # path should end by a "/"
dataset <- args[2]
print("computing :")
print(dataset)

load(paste(path, dataset, "_coex.RData", sep = ""))
coex <- list(output$cell_coex_mats, output$subj_coex_mats)
names(coex) <- c("cell", "subj")
rm(output)

load(paste(path, dataset, "_single_process.RData", sep = ""))

cell_types <- output$cell_types

print("compute col_row_tibble")
col_row_tibbles <- cell_types %>% lapply(function(ct) {
    get_upper_tri_tbl(output$cell_lnks_mats[[ct]])
})
names(col_row_tibbles) <- cell_types
gc(verbose = F)

# =============== NETWORKS PRESERVATION =======================
print("Thresholded")
# Thresholded networks comparison
preserv <- cell_types %>%
    mclapply(function(ct) {
        mat1 <- output$cell_lnks_mats[[ct]]
        mat2 <- output$subj_lnks_mats[[ct]]
        compute_preservation(mat1 = mat1, mat2 = mat2)
    }, mc.cores = 6) %>%
    unlist()
preserv <- tibble(cell_type = cell_types, tp_frac = preserv)

print("Not thresholded")
# AUROC networks comparison
# get coex values for cell upper triangle as vector
preserv_vec1 <- cell_types %>% mclapply(function(ct) {
    coex$cell[[ct]][upper.tri(coex$cell[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(preserv_vec1) <- cell_types %>% sapply(function(ct) {
    paste(ct, "cell", sep = "_")
})

# get coex values for subject upper triangle as vector
preserv_vec2 <- cell_types %>% mclapply(function(ct) {
    coex$subj[[ct]][upper.tri(coex$subj[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(preserv_vec2) <- cell_types %>% sapply(function(ct) {
    paste(ct, "subj", sep = "_")
})
rm(coex)
gc(verbose = F)

# build correlation tibble
preserv_corr_tibles <- cell_types %>% lapply(function(ct) {
    preserv_corr_tible <- cbind(col_row_tibbles[[ct]], preserv_vec1[[paste(ct, "cell", sep = "_")]])
    preserv_corr_tible <- cbind(preserv_corr_tible, preserv_vec2[[paste(ct, "subj", sep = "_")]])
    colnames(preserv_corr_tible) <- c("genes_pair", paste(ct, "cell", sep = "_"), paste(ct, "subj", sep = "_"))
    return(preserv_corr_tible)
})
names(preserv_corr_tibles) <- cell_types

rm(col_row_tibbles)
rm(preserv_vec1)
rm(preserv_vec2)
gc(verbose = F)

print("AUROC computation")
# Get cell level network true positive link gene pairs
cell_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = preserv_corr_tibles[[ct]]$genes_pair,
        links = output$cell_lnks_mats[[ct]][upper.tri(output$cell_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(cell_tps) <- cell_types

# Get subject level network true positive link gene pairs
subj_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = preserv_corr_tibles[[ct]]$genes_pair,
        links = output$subj_lnks_mats[[ct]][upper.tri(output$subj_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(subj_tps) <- cell_types
rm(output)
gc()

print("cell ranking, subj tp")
# cell ranking, subj tp
cell_subj_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- preserv_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, "cell", sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, "cell", sep = "_")], decreasing = T), ]$genes_pair

    auroc(
        ranking = ranked_genes_pair,
        links = subj_tps[[ct]]
    )
}, mc.cores = 6)
names(cell_subj_aurocs) <- cell_types
gc()

print("subj ranking, cell tp")
# subj ranking, cell tp
subj_cell_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- preserv_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, "subj", sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, "subj", sep = "_")], decreasing = T), ]$genes_pair

    auroc(
        ranking = ranked_genes_pair,
        links = cell_tps[[ct]]
    )
}, mc.cores = 6)
names(subj_cell_aurocs) <- cell_types
gc()

rm(preserv_corr_tibles)
gc(verbose = F)
print("Saving")
# Save
output <- list(preserv, cell_subj_aurocs, subj_cell_aurocs)
names(output) <- c("thresh", "ROC_1", "ROC_2")
save(output, file = paste(path, dataset, "_preserv.RData", sep = ""))
