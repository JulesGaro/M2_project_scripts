library(dplyr)
library(tibble)
library(Matrix)
library(stringr)
library(ggplot2)

datasets <- c("lau", "lim", "nagy", "pineda", "ramos", "rosmap", "velmeshev")
cell_types <- c("astrocyte", "excitatory", "inhibitory", "microglia", "oligodendrocyte", "opc")

cell_types <- dataset_tibbles_cell$velmeshev$auroc1$cell_type

load("homogeneity.RData")

# homogeneity plot
ggplot(homogeneity_metrics, aes(x = order, y = mean_of_mean, fill = cell_type)) +
    geom_bar(stat = "identity") +
    xlab("Datasets") +
    ylab("Homogeneity") +
    labs(fill = "Cell types") +
    scale_x_discrete(labels = homogeneity_metrics$dataset) +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)
    ) +
    geom_errorbar(aes(x = order, ymin = mean_of_mean - sd_of_mean, ymax = mean_of_mean + sd_of_mean), width = 0.5, colour = "black", alpha = 0.9, size = 0.75) +
    geom_errorbar(aes(x = order, ymin = mean_of_mean - mean_of_sd, ymax = mean_of_mean + mean_of_sd), width = 0.5, colour = "#414141", alpha = 0.9, size = 0.25)


# Correlation with Reproducibility
load("global_results_cell.RData")

full_cell_repro <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                dataset_tibbles[[dataset]]$tp_frac %>%
                    filter(cell_type == ct) %>%
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
                rep(ct, ncol(dataset_tibbles[[dataset]]$tp_frac) - 2)
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
                dataset_tibbles[[dataset]]$tp_frac %>%
                    filter(cell_type == ct) %>%
                    dplyr::select(!c("cell_type", "preserv")) %>%
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

# adding a row for lau missing opc too avoid bug ...
homogeneity_metrics <- rbind(homogeneity_metrics, c("lau", "opc", 0, 0, 0, 0, "lau_opc", "49"))

full_data$homogeneity <- c(1:nrow(full_data)) %>% sapply(function(i) {
    query <- homogeneity_metrics %>% dplyr::filter((dataset %in% c(full_data$dataset1[i], full_data$dataset2[i])))
    query <- query %>% dplyr::filter((cell_type == full_data$cell_type[i]))
    return(mean(as.numeric(query$mean_of_mean[1]), as.numeric(query$mean_of_mean[2])))
})

# filterout Lau OPC row :
full_data$filtering <- paste(full_data$cell_type, full_data$dataset1)
full_data <- full_data %>% dplyr::filter(filtering != "opc lau")

ggplot(full_data, aes(x = homogeneity, y = cell_repro)) +
    geom_point(aes(color = cell_type)) +
    ggtitle(paste("Pearson correlation :", round(cor(full_data$cell_repro, full_data$homogeneity), 2))) +
    geom_smooth(method = "lm", color = "black") +
    ylab("Cell level reproducibility") +
    xlab("Homogeneity") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    ) +
    labs(color = "Cell types")
