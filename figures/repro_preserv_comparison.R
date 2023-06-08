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
                print(paste(dataset, ct))
                mean(dataset_tibbles_cell[[dataset]]$tp_frac %>% filter(cell_type == ct) %>%
                    select(!c("cell_type", "preserv")) %>% as.vector() %>% unname() %>% sapply(as.numeric))
            }) %>%
            as.vector()
    }) %>%
    as.vector()


full_preserv <- cell_types %>%
    sapply(function(ct) {
        datasets %>%
            sapply(function(dataset) {
                dataset_tibbles_subj[[dataset]]$tp_frac %>%
                    filter(cell_type == ct) %>%
                    select(c("preserv")) %>%
                    as.vector() %>%
                    unname() %>%
                    sapply(as.numeric)
            }) %>%
            as.vector()
    }) %>%
    as.vector()

cell_types_vec <- cell_types %>%
    sapply(function(ct) {
        return(rep(c(ct), length(datasets)))
    }) %>%
    as.vector()

datasets_vec <- rep(datasets, length(cell_types))

full_data <- tibble(
    repro = full_cell_repro,
    preserv = full_preserv,
    cell_type = cell_types_vec,
    dataset = datasets_vec
)

# filterout Lau OPC row :
full_data <- full_data %>% dplyr::filter(preserv > 0)

full_data$dataset <- factor(full_data$dataset)

wilcox.test(full_data$repro, full_data$preserv, alternative = "greater", paired = T)$p.value

ggplot(full_data, aes(x = preserv, y = repro, color = cell_type, shape = dataset)) +
    geom_point(size = 3, stroke = 1) +
    xlim(c(0, 0.45)) +
    ylim(c(0, 0.45)) +
    scale_shape_manual(values = 1:nlevels(full_data$dataset)) +
    geom_hline(yintercept = 0.01, color = "red", linetype = "dotted") +
    geom_vline(xintercept = 0.01, color = "red", linetype = "dotted") +
    xlab("Within-dataset preservation across levels") +
    ylab("Average Cell level reproducibility across other datasets") +
    labs(color = "Cell types", shape = "Datasets") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
    ) +
    geom_abline()
