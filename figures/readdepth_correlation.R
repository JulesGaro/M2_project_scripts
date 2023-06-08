source("../functions/data_process_func.R")
library(dplyr)
library(tibble)
library(Matrix)
library(stringr)
library(ggplot2)

load("global_results_cell.RData")
dataset_tibbles_cell <- dataset_tibbles
load("global_results_subj.RData")
dataset_tibbles_subj <- dataset_tibbles
rm(dataset_tibbles)

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")
cell_types <- c("astrocyte", "excitatory", "inhibitory", "microglia", "oligodendrocyte", "opc")
full_cell_repro <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                dataset_tibbles_cell[[dataset]]$tp_frac %>%
                    dplyr::filter(cell_type == ct) %>%
                    dplyr::select(!c("cell_type", "preserv")) %>%
                    as.vector() %>%
                    unname() %>%
                    sapply(as.numeric)
            }) %>%
            as.vector()
    }) %>%
    as.vector()


full_subj_repro <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                dataset_tibbles_subj[[dataset]]$tp_frac %>%
                    dplyr::filter(cell_type == ct) %>%
                    dplyr::select(!c("cell_type", "preserv")) %>%
                    as.vector() %>%
                    unname() %>%
                    sapply(as.numeric)
            }) %>%
            as.vector()
    }) %>%
    as.vector()

full_cell_types_vec <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                rep(ct, ncol(dataset_tibbles_cell[[dataset]]$tp_frac) - 2)
            }) %>%
            as.vector()
    }) %>%
    as.vector()

full_datasets_vec1 <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(ds) {
                rep(ds, length(datasets) - 1)
            }) %>%
            as.vector()
    }) %>%
    as.vector()

full_datasets_vec2 <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                dataset_tibbles_cell[[dataset]]$tp_frac %>%
                    dplyr::filter(cell_type == ct) %>%
                    select(!c("cell_type", "preserv")) %>%
                    colnames() %>%
                    as.vector() %>%
                    unname() %>%
                    sapply(function(name) {
                        unname(unlist(strsplit(name, "_")))[2]
                    })
            }) %>%
            as.vector()
    }) %>%
    as.vector()

full_data <- tibble(
    subject_repro = full_subj_repro,
    cell_repro = full_cell_repro,
    cell_type = full_cell_types_vec,
    dataset1 = full_datasets_vec1,
    dataset2 = full_datasets_vec2
)

full_dataset_vec_unique <- c(1:nrow(full_data)) %>%
    sapply(function(i) {
        vec <- c(full_data$cell_type[i], full_data$dataset1[i], full_data$dataset2[i]) %>% sort()
        return(paste(vec[1], vec[2], vec[3]))
    }) %>%
    as.vector()

full_data <- full_data[!duplicated(data.frame(full_dataset_vec_unique)), ]

# Read depths metrics
mean_read_depths <- datasets %>% lapply(function(dataset) {
    load(paste(dataset, "/paired_temp/", dataset, "_simple_preprocess.rvar", sep = ""))
    data <- get(paste(dataset, "_output", sep = ""))
    rm(list = c(paste(dataset, "_output", sep = "")))

    data$meta$subject <- data$meta$sample
    data$meta <- data$meta %>% dplyr::select(cell, subject, cell_type)

    # separate ctmat by cell types
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

    # clean ctmat by cell types
    ctmats <- cell_types %>% lapply(function(ct) {
        return(clean_ctmat(ctmats[[ct]]))
    })
    names(ctmats) <- cell_types
    gc()
    # Get rid of cells in metadata as well :
    metas <- cell_types %>% lapply(function(ct) {
        return(metas[[ct]] %>% dplyr::filter(cell %in% colnames(ctmats[[ct]])))
    })
    names(metas) <- cell_types

    data <- list(ctmats, metas)
    names(data) <- c("ctmats", "metas")

    # keep only cell for subjects that pass the minimum cell threshold used
    filtered_ctmats <- cell_types %>% lapply(function(ct) {
        meta <- data$metas[[ct]]
        dims <- meta %>%
            dplyr::select(subject, cell_type, cell) %>%
            unique()
        subjects <- meta %>%
            group_by(subject) %>%
            summarize(n_cel = n())
        subjects <- subjects %>% left_join(dims, by = "subject")
        subjects <- subjects %>% mutate(valid = n_cel > 20)
        valid_subjects <- subjects %>%
            filter(valid) %>%
            as.tibble()
        subjects <- valid_subjects$subject %>%
            unique() %>%
            sort()
        meta <- meta %>% filter(subject %in% subjects)
        return(data$ctmats[[ct]][, meta$cell])
    })
    names(filtered_ctmats) <- cell_types

    mean_read_depth <- cell_types %>% sapply(function(ct) {
        mean(colSums(filtered_ctmats[[ct]]))
    })
    return(tibble(mean_read_depth = mean_read_depth, cell_types = cell_types, dataset = rep(dataset, length(cell_types))))
})
names(mean_read_depths) <- datasets

mean_read_depths <- Reduce(rbind, mean_read_depths)

min_read_depth_mean <- c(1:nrow(full_data)) %>% sapply(function(i) {
    rd_mean1 <- filter(mean_read_depths, cell_types == full_data$cell_type[i]) %>% filter(dataset == full_data$dataset1[i])
    rd_mean1 <- rd_mean1$mean_read_depth

    rd_mean2 <- filter(mean_read_depths, cell_types == full_data$cell_type[i]) %>% filter(dataset == full_data$dataset2[i])
    rd_mean2 <- rd_mean2$mean_read_depth
    return(min(rd_mean1, rd_mean2))
})
min_read_depth_mean

full_data$read_depth <- min_read_depth_mean %>% log()

# get rid of Lau fake OPC (0 reproducibility values)
full_data$filtering <- paste(full_data$dataset1, full_data$cell_type)
full_data <- full_data %>% filter(filtering != "lau opc")

ggplot(data = full_data, aes(x = read_depth, y = cell_repro, color = cell_type)) +
    geom_point() +
    ggtitle(label = paste("Pearson correlation :", round(cor(full_data$cell_repro, full_data$read_depth, method = "pearson"), 2))) +
    geom_smooth(method = "lm", color = "black") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    ) +
    ylab("Cell level reproducibility") +
    xlab("Logarithm of the mean of read depth")
