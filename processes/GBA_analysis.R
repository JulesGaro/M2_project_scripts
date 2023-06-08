source("../functions/EGAD_func.R")

library(EGAD)
library(tidyverse)
library(stringr)
library(parallel)
library(Matrix)

print("Building Annotation Set")
### With Custom GO with BP

args <- commandArgs(trailingOnly = TRUE)
go_annotation_file_path <- args[1]
coexpression_mat_path <- args[2] # path should end by a "/"

GO_unique_filtered <- load_go_annotation(file_path = go_annotation_file_path)

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")

for (dataset in datasets) {
    print(paste("Load :", dataset))
    load(paste(coexpression_mat_path, dataset, "_coex.RData", sep = ""))

    cell_types <- output$cell_types %>% sort()

    print("filtering GO terms...")
    # filter GO term related to at least one genes contained in the network
    ct_annotations <- cell_types %>% mclapply(function(ct) {
        GO_filtered <- dplyr::filter(GO, (GO$GO.ID %in% GO_unique_filtered$GO))
        GO_filtered <- dplyr::filter(GO_filtered, DB_Object_Symbol %in% colnames(output$cell_coex_mats[[ct]]))

        annotations <- make_annotations(
            GO_filtered[, c("DB_Object_Symbol", "GO.ID")],
            unique(GO_filtered$DB_Object_Symbol), unique(GO_filtered$GO.ID)
        )
    }, mc.cores = 6)
    names(ct_annotations) <- cell_types

    print(paste("compute cell neighbor voting for", dataset))
    # Compute neighbor voting for cell level networks
    s <- Sys.time()
    cell_aurocs <- cell_types %>% mclapply(function(ct) {
        neighbor_voting(
            genes.labels = ct_annotations[[ct]],
            network = output$cell_coex_mats[[ct]] %>% as.matrix(),
            nFold = 3, output = "AUROC"
        )
    }, mc.cores = 6)
    names(cell_aurocs) <- cell_types
    print(Sys.time() - s)
    output$cell_coex_mats <- NULL

    print(paste("compute subj neighbor voting for", dataset))
    # Compute neighbor voting for cell level networks
    s <- Sys.time()
    subj_aurocs <- cell_types %>% mclapply(function(ct) {
        neighbor_voting(
            genes.labels = ct_annotations[[ct]],
            network = output$subj_coex_mats[[ct]] %>% as.matrix(),
            nFold = 3, output = "AUROC"
        )
    }, mc.cores = 6)
    names(subj_aurocs) <- cell_types
    print(Sys.time() - s)

    rm(output)
    rm(ct_annotations)

    print(paste("Saving", dataset))
    auroc <- list(cell_aurocs, subj_aurocs, cell_types)
    names(auroc) <- c("cell_aurocs", "subj_aurocs", "cell_types")
    save(auroc, file = paste(coexpression_mat_path, "GBA/", dataset, "_GBA_results.RData", sep = ""))
    rm(auroc)
    gc()
}
