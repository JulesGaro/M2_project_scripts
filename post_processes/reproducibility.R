#!/usr/bin/env Rscript

source("../functions/post_process_func.R")
library(parallel)
library(dplyr)
library(tibble)
library(Matrix)
library(data.table)

setwd("~")

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
dataset_1 <- args[2]
dataset_2 <- args[3]

print("computing :")
print(dataset_1)
print(dataset_2)

paired_file1 <- paste(dataset_1, dataset_2, sep = "_")
paired_file2 <- paste(dataset_2, dataset_1, sep = "_")

load(paste(path, paired_file1, ".rvar", sep = ""))
load(paste(path, paired_file2, ".rvar", sep = ""))

cell_types <- intersect(names(output_1$meta), names(output_2$meta)) %>% sort()

gc(verbose = F)

inter_genes <- cell_types %>% lapply(function(ct) {
    rownames(output_1$cell_lnks_mats[[ct]])
})
names(inter_genes) <- cell_types

load(paste("/space/scratch/jgarreau/data/simple_coex_mats/", dataset_1, "_simple_coex.rvar", sep = ""))
coex1 <- list(output$cell_coex_mats, output$subj_coex_mats) %>% lapply(function(level) {
    level_content <- cell_types %>% mclapply(function(ct) {
        mat <- level[[ct]]
        mat[inter_genes[[ct]], inter_genes[[ct]]]
    }, mc.cores = 6)
    names(level_content) <- cell_types
    return(level_content)
})
names(coex1) <- c("cell", "subj")
rm(output)
load(paste("/space/scratch/jgarreau/data/simple_coex_mats/", dataset_2, "_simple_coex.rvar", sep = ""))
coex2 <- list(output$cell_coex_mats, output$subj_coex_mats) %>% lapply(function(level) {
    level_content <- cell_types %>% mclapply(function(ct) {
        mat <- level[[ct]]
        mat[inter_genes[[ct]], inter_genes[[ct]]]
    }, mc.cores = 6)
    names(level_content) <- cell_types
    return(level_content)
})
names(coex2) <- c("cell", "subj")
rm(output)
gc(verbose = F)

print("compute col_row_tibble")
col_row_tibbles <- cell_types %>% lapply(function(ct) {
    get_upper_tri_tbl(output_1$cell_lnks_mats[[ct]])
})
names(col_row_tibbles) <- cell_types
gc(verbose = F)

# ============== WHOLE NETWORK REPRODUCIBILITY AT THE CELL LEVEL ===============
print("WHOLE NETWORK REPRODUCIBILITY AT THE CELL LEVEL")

print("Thresholded")
# Thresholded networks comparison
repro_cell <- cell_types %>%
    mclapply(function(ct) {
        compute_preservation(mat1 = output_1$cell_lnks_mats[[ct]], mat2 = output_2$cell_lnks_mats[[ct]])
    }, mc.cores = 6) %>%
    unlist()
repro_cell <- tibble(cell_type = cell_types, tp_frac = repro_cell)

print("Not thresholded")
# AUROC networks comparison
# get coex values for 1st dataset upper triangle as vector
repro_cell_vec1 <- cell_types %>% mclapply(function(ct) {
    coex1$cell[[ct]][upper.tri(coex1$cell[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(repro_cell_vec1) <- cell_types %>% sapply(function(ct) {
    paste(ct, dataset_1, sep = "_")
})

# get coex values for 2nd dataset upper triangle as vector
repro_cell_vec2 <- cell_types %>% mclapply(function(ct) {
    coex2$cell[[ct]][upper.tri(coex2$cell[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(repro_cell_vec2) <- cell_types %>% sapply(function(ct) {
    paste(ct, dataset_2, sep = "_")
})
# here
# build correlation tibble
repro_cell_corr_tibles <- cell_types %>% lapply(function(ct) {
    repro_cell_corr_tible <- tibble(
        a = col_row_tibbles[[ct]],
        b = repro_cell_vec1[[paste(ct, dataset_1, sep = "_")]],
        c = repro_cell_vec2[[paste(ct, dataset_2, sep = "_")]]
    )
    colnames(repro_cell_corr_tible) <- c("genes_pair", paste(ct, dataset_1, sep = "_"), paste(ct, dataset_2, sep = "_"))
    return(repro_cell_corr_tible)
})
names(repro_cell_corr_tibles) <- cell_types
rm(repro_cell_vec1)
rm(repro_cell_vec2)
gc(verbose = F)

print("AUROC computation")
# Get cell level network true positive link gene pairs for dataset 1
data1_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = col_row_tibbles[[ct]]$genes_pair,
        links = output_1$cell_lnks_mats[[ct]][upper.tri(output_1$cell_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(data1_tps) <- cell_types

# Get cell level network true positive link gene pairs for dataset 2
data2_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = col_row_tibbles[[ct]]$genes_pair,
        links = output_2$cell_lnks_mats[[ct]][upper.tri(output_2$cell_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(data2_tps) <- cell_types
gc(verbose = F)


# data1 ranking, data2 tp
data1_data2_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- repro_cell_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, dataset_1, sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, dataset_1, sep = "_")], decreasing = T), ]$genes_pair
    auroc(
        ranking = ranked_genes_pair$genes_pair,
        links = data2_tps[[ct]]
    )
}, mc.cores = 6)
names(data1_data2_aurocs) <- cell_types


# data2 ranking, data1 tp
data2_data1_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- repro_cell_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, dataset_2, sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, dataset_2, sep = "_")], decreasing = T), ]$genes_pair

    auroc(
        ranking = ranked_genes_pair$genes_pair,
        links = data1_tps[[ct]]
    )
}, mc.cores = 6)
names(data2_data1_aurocs) <- cell_types

rm(repro_cell_corr_tibles)
coex1$cell <- NULL
coex2cell <- NULL
output_1$cell_lnks_mats <- NULL
output_2$cell_lnks_mats <- NULL
gc(verbose = F)
# Save

output <- list(repro_cell, data1_data2_aurocs, data2_data1_aurocs)
names(output) <- c("thresh", "ROC_1_cell", "ROC_2_cell")
save(output, file = paste(path, dataset_1, "_", dataset_2, "_repro_cell.RData", sep = ""))
rm(output)
gc(verbose = F)


# ============ WHOLE NETWORK REPRODUCIBILITY AT THE SUBJECT LEVEL ==============
print("WHOLE NETWORK REPRODUCIBILITY AT THE SUBJECT LEVEL")

print("Thresholded")
# Thresholded
repro_subj <- cell_types %>%
    mclapply(function(ct) {
        compute_preservation(mat1 = output_1$subj_lnks_mats[[ct]], mat2 = output_2$subj_lnks_mats[[ct]])
    }, mc.cores = 6) %>%
    unlist()
repro_subj <- tibble(cell_type = cell_types, tp_frac = repro_subj)

print("Not thresholded")
# AUROC networks comparison
# get coex values for 1st dataset upper triangle as vector
repro_subj_vec1 <- cell_types %>% mclapply(function(ct) {
    coex1$subj[[ct]][upper.tri(coex1$subj[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(repro_subj_vec1) <- cell_types %>% sapply(function(ct) {
    paste(ct, dataset_1, sep = "_")
})

# get coex values for 2nd dataset upper triangle as vector
repro_subj_vec2 <- cell_types %>% mclapply(function(ct) {
    coex2$subj[[ct]][upper.tri(coex2$subj[[ct]], diag = FALSE)]
}, mc.cores = 6)
names(repro_subj_vec2) <- cell_types %>% sapply(function(ct) {
    paste(ct, dataset_2, sep = "_")
})

repro_subj_corr_tibles <- cell_types %>% lapply(function(ct) {
    repro_subj_corr_tible <- cbind(col_row_tibbles[[ct]], repro_subj_vec1[[paste(ct, dataset_1, sep = "_")]])
    repro_subj_corr_tible <- cbind(repro_subj_corr_tible, repro_subj_vec2[[paste(ct, dataset_2, sep = "_")]])
    colnames(repro_subj_corr_tible) <- c("genes_pair", paste(ct, dataset_1, sep = "_"), paste(ct, dataset_2, sep = "_"))
    return(repro_subj_corr_tible)
})
names(repro_subj_corr_tibles) <- cell_types

rm(repro_subj_vec1)
rm(repro_subj_vec2)
gc(verbose = F)


print("AUROC Computation")
# Get subject level network true positive link gene pairs for dataset 1
data1_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = col_row_tibbles[[ct]]$genes_pair,
        links = output_1$subj_lnks_mats[[ct]][upper.tri(output_1$subj_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(data1_tps) <- cell_types

# Get subject level network true positive link gene pairs for dataset 2
data2_tps <- cell_types %>% mclapply(function(ct) {
    tibble(
        genes_pair = col_row_tibbles[[ct]]$genes_pair,
        links = output_2$subj_lnks_mats[[ct]][upper.tri(output_2$subj_lnks_mats[[ct]], diag = FALSE)]
    ) %>%
        filter(links == 1) %>%
        select(genes_pair) %>%
        as.vector() %>%
        unlist() %>%
        unname()
}, mc.cores = 6)
names(data2_tps) <- cell_types
gc(verbose = F)

# data1 ranking, data2 tp
data1_data2_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- repro_subj_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, dataset_1, sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, dataset_1, sep = "_")], decreasing = T), ]$genes_pair

    auroc(
        ranking = ranked_genes_pair,
        links = data2_tps[[ct]]
    )
}, mc.cores = 6)
names(data1_data2_aurocs) <- cell_types

# data2 ranking, data1 tp
data2_data1_aurocs <- cell_types %>% mclapply(function(ct) {
    # ranking genes pair based on their co-expression values
    ct_coex_tible <- repro_subj_corr_tibles[[ct]] %>% select(genes_pair, paste(ct, dataset_2, sep = "_"))
    ranked_genes_pair <- ct_coex_tible[order(ct_coex_tible[, paste(ct, dataset_2, sep = "_")], decreasing = T), ]$genes_pair

    auroc(
        ranking = ranked_genes_pair,
        links = data1_tps[[ct]]
    )
}, mc.cores = 6)
names(data2_data1_aurocs) <- cell_types

rm(repro_subj_corr_tibles)
gc(verbose = F)

output <- list(repro_subj, data1_data2_aurocs, data2_data1_aurocs)
names(output) <- c("thresh", "ROC_1_subj", "ROC_2_subj")
save(output, file = paste(path, dataset_1, "_", dataset_2, "_repro_subj.RData", sep = ""))
