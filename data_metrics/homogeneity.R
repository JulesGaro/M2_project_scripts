source("~/scripts/functions/data_process_func.R")
library(tidyverse)
library(stringr)
library(Matrix)
library(parallel)
library(ggplot2)
opt <- list(
    GENE_THR = 0.05,
    SAMPLE_THR = 2,
    PERCENTILE_THR = 0.99,
    MIN_N_CELL = 20
)

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")

homogeneity_metrics <- matrix(1, nrow = 1, ncol = 6)
colnames(homogeneity_metrics) <- c("dataset", "cell_type", "mean_of_mean", "sd_of_mean", "mean_of_sd", "sd_of_sd")

for (dataset in datasets) {
    load(paste(dataset, "/paired_temp/", dataset, "_simple_preprocess.rvar", sep = ""))
    data <- get(paste(dataset, "_output", sep = ""))
    rm(list = c(paste(dataset, "_output", sep = "")))
    # preprocessing the data
    cell_types <- data$meta$cell_type %>%
        unique() %>%
        sort()

    # split meta by cell types
    data$metas <- cell_types %>% lapply(function(ct) {
        return(data$meta %>% dplyr::filter(cell_type == ct))
    })
    names(data$metas) <- cell_types

    # split ctmat by cell types, clean and normalize
    data$ctmats <- cell_types %>% lapply(function(ct) {
        ctmat <- data$ctmat[, data$metas[[ct]]$cell]
        return(clean_ctmat(ctmat = ctmat, gene_thr = opt$GENE_THR, sample_thr = opt$SAMPLE_THR) %>% normalize_ctmat())
    })
    names(data$ctmats) <- cell_types

    # get rid of filterout cells in metadatas
    data$metas <- cell_types %>% lapply(function(ct) {
        data$metas[[ct]][data$metas[[ct]]$cell %in% colnames(data$ctmats[[ct]]), ]
    })
    names(data$metas) <- cell_types

    data$ctmat <- data$ctmat %>%
        clean_ctmat(gene_thr = opt$GENE_THR, sample_thr = opt$SAMPLE_THR) %>%
        normalize_ctmat()
    data$meta <- data$meta %>% dplyr::filter(cell %in% colnames(data$ctmat))

    max_size <- min(cell_types %>% sapply(function(ct) {
        nrow(data$metas[[ct]])
    }))
    set.seed(42)
    data$meta <- Reduce(rbind, cell_types %>% lapply(function(ct) {
        meta <- data$meta %>% dplyr::filter(cell_type == ct)
        meta <- meta %>% dplyr::filter(cell %in% sample(meta$cell, size = max_size, replace = F))
    }))
    data$ctmat <- data$ctmat[, data$meta$cell]

    cell_types <- append(cell_types, c("all"))

    data$ctmats <- append(data$ctmats, data$ctmat)
    names(data$ctmats) <- cell_types

    data$metas <- c(data$metas, list(data$meta))
    names(data$metas) <- cell_types

    data <- list(data$ctmats, data$metas)
    names(data) <- c("ctmats", "metas")

    names(data$metas)

    set.seed(42)
    cell_samples <- cell_types %>% mclapply(function(ct) {
        print(ct)
        bs <- c(1:100) %>% lapply(function(i) {
            data$metas[[ct]] %>%
                dplyr::select(cell) %>%
                as.vector() %>%
                unlist() %>%
                unname() %>%
                sample(size = 100, replace = F)
        })
        names(bs) <- c(1:100) %>% sapply(as.character)
        return(bs)
    }, mc.cores = length(cell_types))
    names(cell_samples) <- cell_types


    cells_cor_metrics <- cell_types %>% mclapply(function(ct) {
        bs <- c(1:100) %>%
            sapply(as.character) %>%
            lapply(function(i) {
                cor_mat <- data$ctmats[[ct]][, cell_samples[[ct]][[i]]] %>%
                    WGCNA::cor(method = "pearson") %>%
                    as("dgCMatrix")
                cor_vec <- cor_mat %>% as.vector()
                cor_vec[which(is.na(cor_vec))] <- 0
                metrics <- list(sort(cor_vec), mean(cor_vec), sd(cor_vec))
                names(metrics) <- c("rankcor", "mean", "sd")
                rm(cor_mat)
                return(metrics)
            })
        names(bs) <- c(1:100) %>% sapply(as.character)
        return(bs)
    }, mc.cores = length(cell_types))
    names(cells_cor_metrics) <- cell_types
    rm(data)

    for (ct in cell_types) {
        rankcor_mean <- cbind(c(1:100) %>% sapply(as.character) %>% sapply(function(i) {
            cells_cor_metrics[[ct]][[i]]$rankcor
        }) %>% unname()) %>% rowSums()
        rankcor_mean <- rankcor_mean / 100

        mean_of_mean <- mean(c(1:100) %>% sapply(as.character) %>% sapply(function(i) {
            cells_cor_metrics[[ct]][[i]]$mean
        }) %>% unname())
        sd_of_mean <- sd(c(1:100) %>% sapply(as.character) %>% sapply(function(i) {
            cells_cor_metrics[[ct]][[i]]$mean
        }) %>% unname())
        mean_of_sd <- mean(c(1:100) %>% sapply(as.character) %>% sapply(function(i) {
            cells_cor_metrics[[ct]][[i]]$sd
        }) %>% unname())
        sd_of_sd <- sd(c(1:100) %>% sapply(as.character) %>% sapply(function(i) {
            cells_cor_metrics[[ct]][[i]]$sd
        }) %>% unname())

        homogeneity_metrics <- rbind(c(homogeneity_metrics, dataset, ct, mean_of_mean, sd_of_mean, mean_of_sd, sd_of_sd))
    }
    rm(cells_cor_metrics)
    gc()
}
homogeneity_metrics_1 <- homogeneity_metrics

# Format the homogeneity data
homogeneity_metrics <- matrix(homogeneity_metrics, ncol = 6, nrow = 49, byrow = T)
colnames(homogeneity_metrics) <- c("dataset", "cell_type", "mean_of_mean", "sd_of_mean", "mean_of_sd", "sd_of_sd")
homogeneity_metrics <- homogeneity_metrics[2:nrow(homogeneity_metrics), ] %>% as.tibble()

homogeneity_metrics$mean_of_mean <- homogeneity_metrics$mean_of_mean %>% sapply(as.numeric)
homogeneity_metrics$mean_of_sd <- homogeneity_metrics$mean_of_sd %>% sapply(as.numeric)
homogeneity_metrics$sd_of_mean <- homogeneity_metrics$sd_of_mean %>% sapply(as.numeric)

homogeneity_metrics$dataset_cell_type <- c(1:nrow(homogeneity_metrics)) %>% sapply(function(i) {
    paste(homogeneity_metrics$dataset[i], homogeneity_metrics$cell_type[i], sep = "_")
})

homogeneity_metrics <- homogeneity_metrics[order(homogeneity_metrics$cell_type), ]
homogeneity_metrics$order <- c(1:nrow(homogeneity_metrics)) %>% sapply(as.factor)


save(homogeneity_metrics, file = "homogeneity.RData")
