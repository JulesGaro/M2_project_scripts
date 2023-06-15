#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(Matrix)
library(stringr)

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")
path_repro <- ""
path_preserv <- ""

dataset_tibbles <- datasets %>% lapply(function(ds) {
    # Load preserv data
    setwd(path_preserv)
    load(paste(ds, "preserv.rvar", sep = "_"))
    preserv_data <- output
    rm(output)

    # start building tibbles from preservation data
    cell_types <- preserv_data$thresh$cell_type
    tp_frac <- tibble(cell_type = cell_types, preserv = preserv_data$thresh$tp_frac)
    auroc1 <- tibble(
        cell_type = cell_types,
        preserv_roc1 = cell_types %>% sapply(function(ct) {
            preserv_data$ROC_1[[ct]][["auroc.W"]]
        }) %>% unname()
    )
    auroc2 <- tibble(
        cell_type = cell_types,
        preserv_roc2 = cell_types %>% sapply(function(ct) {
            preserv_data$ROC_2[[ct]][["auroc.W"]]
        }) %>% unname()
    )
    rm(preserv_data)

    # handle Lau lack of OPC
    if (ds == "lau") {
        tp_frac <- rbind(tp_frac, c("opc", 0))
        auroc1 <- rbind(auroc1, c("opc", 0))
        auroc2 <- rbind(auroc2, c("opc", 0))
    }

    # add repro data to the created tibbles
    setwd(path_repro)
    files <- list.files()
    files <- files[grepl(ds, files)]
    files <- files[grepl("repro_subj", files)]

    # retrieve other datasets names
    other_datasets <- files %>% sapply(function(file) {
        split_file <- str_split(file, "_") %>% unlist()
        ds1 <- split_file[1]
        ds2 <- split_file[2]
        if (ds == ds1) {
            return(ds2)
        } else {
            return(ds1)
        }
    })

    # get repro data for each other datasets
    for (file in files) {
        load(file)
        repro_data <- output
        rm(output)

        # Get values
        other_tp_frac <- repro_data$thresh$tp_frac
        other_auroc1 <- cell_types %>%
            sapply(function(ct) {
                repro_data$ROC_1_subj[[ct]][["auroc.W"]]
            }) %>%
            unname()
        other_auroc2 <- cell_types %>%
            sapply(function(ct) {
                repro_data$ROC_2_subj[[ct]][["auroc.W"]]
            }) %>%
            unname()

        # Handle Lau lack of OPC
        if (length(repro_data$thresh$cell_type) == 5) {
            other_tp_frac <- append(other_tp_frac, 0)
            other_auroc1 <- append(other_auroc1, 0) %>% unlist()
            other_auroc2 <- append(other_auroc2, 0) %>% unlist()
        }
        # Add values to the tibbles
        tp_frac <- cbind(tp_frac, other_tp_frac)
        auroc1 <- cbind(auroc1, other_auroc1)
        auroc2 <- cbind(auroc2, other_auroc2)
    }

    # names the tibbles columns
    colnames(tp_frac)[3:8] <- paste(rep(c("repro"), 6), other_datasets, sep = "_")
    colnames(auroc1)[3:8] <- paste(rep(c("repro_roc1"), 6), other_datasets, sep = "_")
    colnames(auroc2)[3:8] <- paste(rep(c("repro_roc2"), 6), other_datasets, sep = "_")
    dataset_list <- list(tp_frac, auroc1, auroc2)
    names(dataset_list) <- c("tp_frac", "auroc1", "auroc2")

    return(dataset_list)
})
names(dataset_tibbles) <- datasets

save(dataset_tibbles, file = "global_results_subj.RData")
