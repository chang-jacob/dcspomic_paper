## ---------------------------
##
## Script name: crc_dcspomic.R
##
## Author: Jake Chang
##
## Date Modified: 2025-08-19
##
## ---------------------------
##
#' Description:
#' This script performs differential colocalization analysis on a colorectal
#' cancer TMA previously published in Cell, 2020 by Schurch et al:
#' (https://shorturl.at/jMCXx).
#'
#' Data processing of the TMA follows the notebooks published here:
#' (https://github.com/nickcee/spicyRPaper/blob/main/CODEX_ColonCancer_Schurch.Rmd)
##
## ---------------------------

## load libraries
library(spomic)
library(dcspomic)

library(metafor)
library(tidyverse)
library(tidyplots)
library(pbapply)
library(beepr)
library(janitor)

scratch <- config::get("scratch")
oak <- config::get("oak")
date <- format(Sys.Date(), "%Y%m%d")

source("analysis/plotting_utils.R")

# Group 1: Crohn's like reaction (CLR)
# Group 2: diffuse inflammatory infiltration (DII)

## ---------------------------

## run DC-SPOMIC

group1_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group1_spomics_Lcross_global_envelope.rds"))
group2_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group2_spomics_Lcross_global_envelope.rds"))

dcspomic <- create_dcspomic(
  group1_name = "CLR",
  group1_spomics = group1_spomics,
  group2_name = "DII",
  group2_spomics = group2_spomics)
dcspomic <- set_dcspomic_hyperparameters(
  dcspomic,
  r = 75, # 75um
  colocalization_type = "Lcross",
  tau_estimator = "SJ")

dcspomic <- run_dcspomic(dcspomic)

# Remove the cell pairs containing undefined cells and recompute adjusted p-vals
results_without_unknowns <- dcspomic@results$differential_testing[
  !grepl("\\(undefined\\)", dcspomic@results$differential_testing$i_j),
]
results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
results_without_unknowns$holm <- p.adjust(results_without_unknowns$pval, method="holm")
dcspomic@results$differential_testing <- results_without_unknowns

alpha <- 0.05
plot_volcano(dcspomic, alpha = alpha)

tidyplots::save_plot(filename = "output/crc_tma/lcross_75um/volcano.png", units = "in")
tidyplots::save_plot(filename = "output/crc_tma/lcross_75um/volcano.pdf", units = "in")
tidyplots::save_plot(filename = "output/crc_tma/lcross_75um/volcano.svg", units = "in")

## ---------------------------

generate_forest_plots <- function(dcspomic, i, j) {
  get_forest_data <- function(model, group_name) {
    yi <- as.vector(model$yi)
    vi <- model$vi
    slabs <- model$slab
    weights_vec <- weights(model)

    data.frame(
      yi = yi,
      vi = vi,
      ci.lb = yi - 1.96 * sqrt(vi),
      ci.ub = yi + 1.96 * sqrt(vi),
      study = slabs,
      group = group_name,
      weight = weights_vec,
      is_summary = FALSE,
      stringsAsFactors = FALSE
    )
  }

  get_summary_row <- function(model, group_name) {
    est <- as.numeric(model$b)
    se <- model$se
    data.frame(
      yi = est,
      vi = NA,
      ci.lb = est - 1.96 * se,
      ci.ub = est + 1.96 * se,
      study = "Pooled Estimate",
      group = group_name,
      weight = NA,
      is_summary = TRUE,
      stringsAsFactors = FALSE
    )
  }
  pair <- paste0(i, "_",  j)
  mod1 <- dcspomic@results$group1_models[[pair]]
  mod2 <- dcspomic@results$group2_models[[pair]]

  df1 <- get_forest_data(mod1, "Group 1")
  df2 <- get_forest_data(mod2, "Group 2")
  sum1 <- get_summary_row(mod1, "Group 1")
  sum2 <- get_summary_row(mod2, "Group 2")

  df <- bind_rows(df1, df2, sum1, sum2)
  # Plot position
  df <- df %>%
    mutate(study_id = paste(study, group),
           y_pos = rev(seq_len(n()))
           )
  df <- df %>%
    mutate(weight = ifelse(is.na(weight), 0, weight))  # or assign fixed like 1 or 5
  g1_pool <- df |> filter(study == "Pooled Estimate", group == "Group 1") |> pull(yi)
  g2_pool <- df |> filter(study == "Pooled Estimate", group == "Group 2") |> pull(yi)

  ggplot(df, aes(x = yi, y = y_pos)) +

    # Error bars for individual studies (gray)
    geom_errorbarh(
      data = df %>% filter(!is_summary),
      aes(xmin = ci.lb, xmax = ci.ub),
      height = 0.2,
      color = "gray50"
    ) +

    # Error bars for summary points (colored)
    geom_errorbarh(
      data = df %>% filter(is_summary),
      aes(xmin = ci.lb, xmax = ci.ub, color = group),
      height = 0.2,
      linewidth = 0.9
    ) +

    # Open dots for studies
    geom_point(
      data = df %>% filter(!is_summary),
      aes(size = weight, fill = group),
      shape = 21, stroke = 0.6, color = "gray40"
    ) +

    # Diamonds for summary points
    geom_point(
      data = df %>% filter(is_summary),
      aes(color = group),
      shape = 18, size = 6
    ) +
    geom_vline(xintercept=g1_pool,color = "blue", linetype = "dotted") +
    geom_vline(xintercept=g2_pool,color = "red", linetype = "dotted") +
    geom_vline(xintercept=dcspomic@details$hyperparameters$r,color = "gray", linetype = "dotted") +
    scale_fill_manual(values = c("Group 1" = "lightblue", "Group 2" = "lightcoral"),   guide = "none"
) +
    scale_color_manual(values = c("Group 1" = "blue", "Group 2" = "red"),
                       labels = c("CLR", "DII")) +

    scale_size_continuous(range = c(0.5, 5)) +

    scale_y_continuous(
      breaks = df$y_pos,
      labels = df$study,
      expand = expansion(mult = c(0.01, 0.01))
    ) +

    labs(
      x = "Colocalization (Lcross: r=75 um)",
      y = NULL,
      color = "Group",
      fill = "Group",
      size = "Weight",
      title = paste(pair, "Forest Plot")
    ) +

    ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=5),
        legend.title = element_text(size = 5),
        title = element_text(size = 7),
        legend.position = "right")
}

i = "(b_cells)"
j = "(b_cells)"
generate_forest_plots(dcspomic, i, j)
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.pdf"), height = 5, width = 4, units = "in")
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.png"), height = 5, width = 4, units = "in")
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.svg"), height = 5, width = 4, units = "in")

p1 <- plot_cell_pair(group1_spomics[[6]], i, j)
p2 <- plot_cell_pair(group2_spomics[[10]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")

tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_bcells_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")

###########

i = "(cd4_t_cells_cd45ro)"
j = "(cd4_t_cells_cd45ro)"
generate_forest_plots(dcspomic, i, j)
ggsave("output/crc_tma/lcross_75um/cd4_cd4_forestplot.pdf", height = 5, width = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/cd4_cd4_forestplot.png", height = 5, width = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/cd4_cd4_forestplot.svg", height = 5, width = 4, units = "in")

p1 <- plot_cell_pair(group1_spomics[[3]], i, j)
p2 <- plot_cell_pair(group2_spomics[[11]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd4tcells_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")
#####
i = "(cd4_t_cells_cd45ro)"
j = "(tregs)"
generate_forest_plots(dcspomic, i, j)
ggsave("output/crc_tma/lcross_75um/cd4_tregs_forestplot.pdf", height = 5, width = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/cd4_tregs_forestplot.png", height = 5, width = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/cd4_tregs_forestplot.svg", height = 5, width = 4, units = "in")


p1 <- plot_cell_pair(group1_spomics[[2]], i, j)
p2 <- plot_cell_pair(group2_spomics[[2]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_tregs_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")
#####
i = "(cd4_t_cells_cd45ro)"
j = "(cd8_t_cells)"
generate_forest_plots(dcspomic, i, j)
p1 <- plot_cell_pair(group1_spomics[[12]], i, j)
p2 <- plot_cell_pair(group2_spomics[[11]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")

tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")
#####
i = "(b_cells)"
j = "(cd8_t_cells)"
generate_forest_plots(dcspomic, i, j)
p1 <- plot_cell_pair(group1_spomics[[6]], i, j)
p2 <- plot_cell_pair(group2_spomics[[18]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")

tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/bcells_cd8tcells_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")

#####
i = "(b_cells)"
j = "(cd4_t_cells_cd45ro)"
generate_forest_plots(dcspomic, i, j)
p1 <- plot_cell_pair(group1_spomics[[8]], i, j)
p2 <- plot_cell_pair(group2_spomics[[1]], i, j)

tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p1,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group1.svg",
                     width = 4,
                     height = 4,
                     units = "in")

tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.png",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.pdf",
                     width = 4,
                     height = 4,
                     units = "in")
tidyplots::save_plot(plot = p2,
                     filename = "output/crc_tma/lcross_75um/cd4tcells_cd8tcells_group2.svg",
                     width = 4,
                     height = 4,
                     units = "in")

## ---------------------------