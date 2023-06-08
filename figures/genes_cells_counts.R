library(tidyverse)
library(Matrix)
library(tibble)
source("../functions/data_process_func.R")

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")

global_data_metrics <- datasets %>% lapply(function(dataset) {
    load(paste(dataset, "_preprocess.RData", sep = ""))
    data <- get(paste(dataset, "_output", sep = ""))
    rm(list = c(paste(dataset, "_output", sep = "")))

    data$meta$subject <- data$meta$sample
    data$meta <- data$meta %>% dplyr::select(cell, subject, cell_type)

    cell_types <- data$meta$cell_type %>%
        unique() %>%
        sort()
    metas <- cell_types %>% lapply(function(ct) {
        return(data$meta %>% dplyr::filter(cell_type == ct))
    })
    names(metas) <- cell_types
    ctmats <- cell_types %>% lapply(function(ct) {
        return(data$ctmat[, metas[[ct]]$cell])
    })
    names(ctmats) <- cell_types

    data <- list(ctmats, metas)
    names(data) <- c("ctmats", "metas")

    data$ctmats <- cell_types %>% lapply(function(ct) {
        return(normalize_ctmat(ctmat = clean_ctmat(ctmat = data$ctmats[[ct]])))
    })
    names(data$ctmats) <- cell_types

    data$metas <- cell_types %>% lapply(function(ct) {
        return(data$metas[[ct]] %>% dplyr::filter(cell %in% colnames(data$ctmats[[ct]])))
    })
    names(data$metas) <- cell_types

    n_subject <- data$metas$astrocyte$subject %>%
        unique() %>%
        length()
    n_genes <- cell_types %>% sapply(function(ct) {
        return(nrow(data$ctmats[[ct]]))
    })
    names(n_genes) <- cell_types

    genes <- cell_types %>% sapply(function(ct) {
        return(rownames(data$ctmats[[ct]]))
    })
    names(genes) <- cell_types

    n_cells <- cell_types %>% sapply(function(ct) {
        return(nrow(data$metas[[ct]]))
    })
    names(n_cells) <- cell_types

    sparsity <- cell_types %>% sapply(function(ct) {
        return(length(data$ctmats[[ct]]@x) / length(data$ctmats) * 100)
    })
    names(sparsity) <- cell_types

    metrics <- list(n_genes, n_cells, n_subject, sparsity, genes)
    names(metrics) <- c("n_genes", "n_cells", "n_subject", "sparsity", "genes")
    return(metrics)
})
names(global_data_metrics) <- datasets


all_genes <- Reduce(c, datasets %>% lapply(function(ds) {
    n_genes <- global_data_metrics[[ds]]$n_genes
    names(n_genes) <- paste(rep(ds, length(n_genes)), names(n_genes), sep = "_")
    return(n_genes)
}))

all_ct <- Reduce(c, datasets %>% lapply(function(ds) {
    return(names(global_data_metrics[[ds]]$n_genes))
}))
all_ct

all_genes <- tibble(id = names(all_genes), cell_type = all_ct, n_genes = all_genes)
all_genes <- arrange(all_genes, cell_type)

ggplot(data = all_genes, aes(x = id, y = n_genes, fill = cell_type)) +
    geom_bar(stat = "identity") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "none"
    ) +
    labs(fill = "Cell types") +
    xlab("Dataset/Cell types") +
    ylab("Number of genes") +
    scale_x_discrete(limits = reorder(all_genes$id, all_genes$cell_type)) +
    coord_flip()

all_cells <- Reduce(c, datasets %>% lapply(function(ds) {
    n_cells <- global_data_metrics[[ds]]$n_cells
    names(n_cells) <- paste(rep(ds, length(n_cells)), names(n_cells), sep = "_")
    return(n_cells)
}))

all_cells <- tibble(id = names(all_cells), cell_type = all_ct, n_cells = all_cells)
all_cells <- arrange(all_cells, cell_type)

ggplot(data = all_cells, aes(x = id, y = n_cells, fill = cell_type)) +
    geom_bar(stat = "identity") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "none"
    ) +
    labs(fill = "Cell types") +
    xlab("Dataset/Cell types") +
    ylab("Number of cells") +
    scale_x_discrete(limits = reorder(all_cells$id, all_cells$cell_type)) +
    coord_flip()
