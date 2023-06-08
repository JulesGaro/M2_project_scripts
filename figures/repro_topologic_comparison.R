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
                    filter(cell_type == ct) %>%
                    select(!c("cell_type", "preserv")) %>%
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
                    filter(cell_type == ct) %>%
                    select(!c("cell_type", "preserv")) %>%
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

full_datasets_vec <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                rep(dataset, ncol(dataset_tibbles_cell[[dataset]]$tp_frac) - 2)
            }) %>%
            as.vector()
    }) %>%
    as.vector()


full_data <- tibble(
    subj = full_subj_repro,
    cell = full_cell_repro,
    cell_type = full_cell_types_vec,
    dataset = full_datasets_vec
)

# filterout Lau OPC row :
full_data <- full_data %>% dplyr::filter(cell > 0)

# remove doublon and set ds1 and ds2
full_data <- full_data[order(full_data$cell, decreasing = T), ]
ds2 <- full_data[seq(2, nrow(full_data), 2), ]$dataset
full_data <- full_data[seq(1, nrow(full_data), 2), ]
full_data$dataset <- paste(full_data$dataset, "/", ds2, sep = "")


full_data$mean <- c(1:nrow(full_data)) %>% sapply(function(i) {
    return(mean(filter(full_data, cell_type == full_data$cell_type[i])$cell))
})

full_data$sd <- c(1:nrow(full_data)) %>% sapply(function(i) {
    return(sd(filter(full_data, cell_type == full_data$cell_type[i])$cell))
})

mean_data <- full_data %>%
    select(mean, sd, cell_type) %>%
    distinct()

full_data$mean <- c(1:nrow(full_data)) %>% sapply(function(i) {
    return(mean(filter(full_data, cell_type == full_data$cell_type[i])$subj))
})

full_data$sd <- c(1:nrow(full_data)) %>% sapply(function(i) {
    return(sd(filter(full_data, cell_type == full_data$cell_type[i])$subj))
})

subj_mean_data <- full_data %>%
    select(mean, sd, cell_type) %>%
    distinct()

mean_data$level <- rep(c("cell"), 6)
subj_mean_data$level <- rep(c("subject"), 6)
mean_data <- rbind(mean_data, subj_mean_data)
mean_data$factor <- paste(mean_data$cell_type, mean_data$level)

max_repro <- max(union(full_subj_repro, full_cell_repro))
plotlim <- round(max_repro + (max_repro / 10), digits = 1)

wilcox.test(full_data$cell, full_data$subj, alternative = "greater", paired = T)$p.value

full_data$dataset <- factor(full_data$dataset)

ggplot(full_data, aes(x = subj, y = cell, color = cell_type)) +
    geom_point(size = 3) +
    xlim(c(0, plotlim)) +
    ylim(c(0, plotlim)) +
    geom_hline(yintercept = 0.01, color = "red", linetype = "dotted") +
    geom_vline(xintercept = 0.01, color = "red", linetype = "dotted") +
    xlab("Subject level reproducibility") +
    ylab("Cell level reproducibility") +
    labs(color = "Cell-types") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
    ) +
    geom_abline()


ggplot(mean_data, aes(x = cell_type, y = mean, fill = level)) +
    xlab("Cell types") +
    ylab("Mean of reproducibility") +
    labs(fill = "Level") +
    geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
    geom_errorbar(position = position_dodge(width = 0.9), aes(x = cell_type, ymin = mean - sd, ymax = mean + sd), width = 0.4, colour = "black", size = 0.5) +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)
    )
