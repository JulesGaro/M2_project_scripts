library(tidyverse)
library(stringr)
library(parallel)
library(Matrix)
library(ggplot2)
library(GO.db)

# Load EGAD results
setwd("GBA/")
files <- list.files()
# files <- c("aggregate_EGAD_aurocs.rvar")

aurocs <- files %>% lapply(function(file) {
    load(file)
    names(auroc$cell_aurocs) <- auroc$cell_type
    names(auroc$subj_aurocs) <- auroc$cell_type
    add <- auroc
    rm(auroc)
    return(add)
})
datasets <- files %>%
    sapply(function(file) {
        unlist(str_split(file, pattern = "_"))[1]
    }) %>%
    unname()
names(aurocs) <- datasets

cell_types <- aurocs$velmeshev$cell_types %>% sort()

# Build general information table
average_auc_cell <- datasets %>%
    sapply(function(dataset) {
        cell_types %>% sapply(function(ct) {
            values <- aurocs[[dataset]]$cell_aurocs[[ct]][, "auc"] %>% unname()
            values[is.na(values)] <- 0
            # values <- values[values > 0.6]
            return(mean(values))
        })
    }) %>%
    unname() %>%
    as.vector()

percent_above_abline_cell <- datasets %>%
    sapply(function(dataset) {
        cell_types %>% sapply(function(ct) {
            values <- aurocs[[dataset]]$cell_aurocs[[ct]][, "auc"] %>% unname()
            values[is.na(values)] <- 0
            return(length(which(values > 0.6)) / length(values) * 100)
        })
    }) %>%
    unname() %>%
    as.vector()

average_auc_subj <- datasets %>%
    sapply(function(dataset) {
        cell_types %>% sapply(function(ct) {
            values <- aurocs[[dataset]]$subj_aurocs[[ct]][, "auc"] %>% unname()
            values[is.na(values)] <- 0
            # values <- values[values > 0.6]
            return(mean(values))
        })
    }) %>%
    unname() %>%
    as.vector()

EGAD_metrics <- tibble(
    cell_type = rep(cell_types, length(datasets)),
    dataset = datasets %>% sapply(function(dataset) {
        rep(dataset, length(cell_types))
    }) %>% as.vector(),
    average_auc_cell = average_auc_cell,
    average_auc_subj = average_auc_subj
)

EGAD_metrics$dataset <- factor(EGAD_metrics$dataset)

# -----------Plot AUROC for all the data------------------------
wilcox.test(EGAD_metrics$average_auc_cell, EGAD_metrics$average_auc_subj, alternative = "greater", paired = T)$p.value

ggplot(data = EGAD_metrics, aes(x = average_auc_subj, y = average_auc_cell, color = cell_type, shape = dataset)) +
    scale_shape_manual(values = 1:nlevels(EGAD_metrics$dataset)) +
    geom_point(size = 2, stroke = 1) +
    geom_abline() +
    geom_hline(yintercept = 0.6, color = "purple", linetype = "dotted") +
    geom_vline(xintercept = 0.6, color = "purple", linetype = "dotted") +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    ylim(0.45, 0.6) +
    xlim(0.45, 0.6) +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#eeeeee"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
    ) +
    ylab("Average AUROC at the cell level") +
    xlab("Average AUROC at the subject level") +
    labs(color = "Cell types", shape = "Dataset")

# -----------Functions type prediction differences btw cell and subj level-----------
go <- keys(GO.db, keytype = "GOID")
GO <- select(GO.db, columns = c("GOID", "TERM"), keys = go, keytype = "GOID")
GOID <- GO$GOID
GO <- GO$TERM
names(GO) <- GOID

cell_types_datasets <- cell_types %>% lapply(function(ct) {
    ct_datasets <- c()
    for (dataset in datasets) {
        if (ct %in% aurocs[[dataset]]$cell_types) {
            ct_datasets <- append(ct_datasets, dataset)
        }
    }
    return(ct_datasets)
})
names(cell_types_datasets) <- cell_types

average_auc <- c("cell_aurocs", "subj_aurocs") %>% lapply(function(level) {
    ct_average_auc <- cell_types %>% lapply(function(ct) {
        datasets <- cell_types_datasets[[ct]]
        ct_aucs <- datasets %>% lapply(function(ds) {
            unified_go <- Reduce(union, datasets %>% lapply(function(ds) {
                rownames(aurocs[[ds]][[level]][[ct]])
            }))
            aucs <- aurocs[[ds]][[level]][[ct]]
            aucs <- cbind(aucs, rownames(aucs))
            colnames(aucs) <- c("auc", "avg_node_degree", "degree_nul_auc", "go")
            aucs <- aucs %>% as.tibble()
            aucs <- aucs %>% dplyr::select(auc, go)
            inter_go <- unified_go[!(unified_go %in% aucs$go)]
            added_go <- tibble(go = inter_go, auc = rep(c(0), length(inter_go)))
            aucs <- rbind(aucs, added_go)
            aucs <- aucs[order(aucs$go), ]
            return(aucs)
        })
        names(ct_aucs) <- datasets
        ct_go <- ct_aucs$lim$go
        ct_aucs <- datasets %>% lapply(function(ds) {
            return(ct_aucs[[ds]]$auc)
        })
        ct_aucs <- ct_aucs %>% as.data.frame()

        ct_average_auc <- c(1:nrow(ct_aucs)) %>% sapply(function(i) {
            ct_aucs[i, ] %>%
                unname() %>%
                as.vector() %>%
                unlist() %>%
                as.numeric() %>%
                mean()
        })
        go_terms <- ct_go %>% sapply(function(id) {
            GO[[id]]
        })
        ct_average_auc <- tibble(term = go_terms, go = ct_go, auc = ct_average_auc)
        return(ct_average_auc)
    })
    names(ct_average_auc) <- cell_types
    return(ct_average_auc)
})
names(average_auc) <- c("cell_aurocs", "subj_aurocs")

# Excitatory
top_cell <- average_auc$cell_aurocs$excitatory %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$excitatory %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$excitatory %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$excitatory %>% dplyr::filter(auc >= 0.6)

ex_top_left <- average_auc$cell_aurocs$excitatory %>% filter(go %in% intersect(top_cell$go, top_subj$go))
ex_bot_right <- average_auc$subj_aurocs$excitatory %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
ex_top_right <- average_auc$cell_aurocs$excitatory %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))

# Inhibitory
top_cell <- average_auc$cell_aurocs$inhibitory %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$inhibitory %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$inhibitory %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$inhibitory %>% dplyr::filter(auc >= 0.6)

in_top_left <- average_auc$cell_aurocs$inhibitory %>% filter(go %in% intersect(top_cell$go, top_subj$go))
in_bot_right <- average_auc$subj_aurocs$inhibitory %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
in_top_right <- average_auc$cell_aurocs$inhibitory %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))


# Microglia
top_cell <- average_auc$cell_aurocs$microglia %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$microglia %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$microglia %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$microglia %>% dplyr::filter(auc >= 0.6)

mi_top_left <- average_auc$cell_aurocs$microglia %>% filter(go %in% intersect(top_cell$go, top_subj$go))
mi_bot_right <- average_auc$subj_aurocs$microglia %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
mi_top_right <- average_auc$cell_aurocs$microglia %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))


# Astrocytes
top_cell <- average_auc$cell_aurocs$astrocyte %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$astrocyte %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$astrocyte %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$astrocyte %>% dplyr::filter(auc >= 0.6)

as_top_left <- average_auc$cell_aurocs$astrocyte %>% filter(go %in% intersect(top_cell$go, top_subj$go))
as_bot_right <- average_auc$subj_aurocs$astrocyte %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
as_top_right <- average_auc$cell_aurocs$astrocyte %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))


# Oligodendrocytes
top_cell <- average_auc$cell_aurocs$oligodendrocyte %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$oligodendrocyte %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$oligodendrocyte %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$oligodendrocyte %>% dplyr::filter(auc >= 0.6)

ol_top_left <- average_auc$cell_aurocs$oligodendrocyte %>% filter(go %in% intersect(top_cell$go, top_subj$go))
ol_bot_right <- average_auc$subj_aurocs$oligodendrocyte %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
ol_top_right <- average_auc$cell_aurocs$oligodendrocyte %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))


# OPC
top_cell <- average_auc$cell_aurocs$opc %>% dplyr::filter(auc >= 0.6)
top_subj <- average_auc$subj_aurocs$opc %>% dplyr::filter(auc < 0.5)
bottom_cell <- average_auc$cell_aurocs$opc %>% dplyr::filter(auc < 0.5)
bottom_subj <- average_auc$subj_aurocs$opc %>% dplyr::filter(auc >= 0.6)

opc_top_left <- average_auc$cell_aurocs$opc %>% filter(go %in% intersect(top_cell$go, top_subj$go))
opc_bot_right <- average_auc$subj_aurocs$opc %>% filter(go %in% intersect(bottom_cell$go, bottom_subj$go))
opc_top_right <- average_auc$cell_aurocs$opc %>% filter(go %in% intersect(top_cell$go, bottom_subj$go))


# Focus on top right term (well retrieve in both levels)
common_top_right <- c(
    ex_top_right$go,
    in_top_right$go,
    mi_top_right$go,
    as_top_right$go,
    ol_top_right$go,
    opc_top_right$go
)

common_top_right <- table(common_top_right) %>%
    as.data.frame() %>%
    as.tibble()
common_top_right <- common_top_right %>% dplyr::filter(Freq == 6)
common_top_right$go <- common_top_right$common_top_right %>% as.vector()
common_top_right <- common_top_right %>% dplyr::select(go, Freq)
common_top_right$term <- common_top_right$go %>% sapply(function(id) {
    GO[[id]]
})

bot_right <- list(as_bot_right, ex_bot_right, in_bot_right, mi_bot_right, ol_bot_right, opc_bot_right)
names(bot_right) <- cell_types

top_left <- list(as_top_left, ex_top_left, in_top_left, mi_top_left, ol_top_left, opc_top_left)
names(top_left) <- cell_types

# Label position of GO terms
for (ct in cell_types) {
    average_auc$cell_aurocs[[ct]]$position <- average_auc$cell_aurocs[[ct]]$go %>%
        sapply(function(id) {
            if (id %in% bot_right[[ct]]$go) {
                return("subject")
            } else {
                return("other")
            }
        })
    average_auc$cell_aurocs[[ct]]$position <- c(1:length(average_auc$cell_aurocs[[ct]]$go)) %>%
        sapply(function(i) {
            if (average_auc$cell_aurocs[[ct]]$go[i] %in% top_left[[ct]]$go) {
                return("cell")
            } else {
                return(average_auc$cell_aurocs[[ct]]$position[i])
            }
        })
    average_auc$cell_aurocs[[ct]]$position <- c(1:length(average_auc$cell_aurocs[[ct]]$go)) %>%
        sapply(function(i) {
            if (average_auc$cell_aurocs[[ct]]$go[i] %in% common_top_right$go) {
                return("common")
            } else {
                return(average_auc$cell_aurocs[[ct]]$position[i])
            }
        })
}

# excitatory  inhibitory  astrocyte microglia oligodendrocyte opc
cell_type <- "astrocyte"

ggplot(mapping = aes(
    x = average_auc$subj_aurocs[[cell_type]]$auc,
    y = average_auc$cell_aurocs[[cell_type]]$auc,
    color = average_auc$cell_aurocs[[cell_type]]$position,
)) +
    geom_point(size = 5) +
    ggtitle(cell_type) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0.6, color = "purple", linetype = "dotted") +
    geom_vline(xintercept = 0.6, color = "purple", linetype = "dotted") +
    xlim(0.2, 1) +
    ylim(0.2, 1) +
    labs(color = "Specificity") +
    ylab("Average AUROC at the cell level") +
    xlab("Average AUROC at the subject level") +
    theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
    ) +
    scale_color_manual(values = c("blue", "purple", "black", "red")) +
    geom_abline()
