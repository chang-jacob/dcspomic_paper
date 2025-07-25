## ---------------------------
##
## Script name: crc_dcspomic.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-10
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
# devtools::load_all("/Users/jacobchang/Lab/spomic")
# devtools::load_all("/Users/jacobchang/Lab/dcspomic")

devtools::load_all("/home/groups/plevriti/jachang4/spomic")
devtools::load_all("/home/groups/plevriti/jachang4/dcspomic")


library(metafor)
library(tidyverse)
library(tidyplots)
library(pbapply)
library(beepr)
library(janitor)

source("dcspomic_paper/analysis/plotting_utils.R")

is_sherlock <- grepl("sh", Sys.info()[["nodename"]]) || Sys.getenv("SHERLOCK") != ""


output_path <- ifelse(is_sherlock,
                      "/oak/stanford/groups/plevriti/jachang4/dcspomic_output",
                      "output")
plot_path <- ifelse(is_sherlock,
                    "/home/groups/plevriti/jachang4/dcspomic_paper/output",
                    "output")

# Schurch et al. has a curated TMA with two groups which we will perform
# differential colocalization analysis on

# Group 1: Crohn's like reaction (CLR)
# Group 2: diffuse inflammatory infiltration (DII)

## ---------------------------

## Create spomic objects for each TMA core from Schurch et al.

# The data is processed following the spicyR manuscript for fair comparison
# https://github.com/nickcee/spicyRPaper/blob/main/CODEX_ColonCancer_Schurch.Rmd

df <- read.csv("data/crc_tma/processed_crc_tma.csv") |>
  dplyr::select(-X) |>
  janitor::clean_names()
head(df)

codex_df <- df |>
  dplyr::select(cell_id, cell_type, group, histology, patient, x, y) |>
  dplyr::mutate(cell_type = paste0("(", cell_type, ")")) |>
  # convert pixels to um
  dplyr::mutate(x = x*0.37742,
                y = y*0.37742)

# generate sample overview plot
codex_df |>
  dplyr::mutate(patient = as.factor(patient)) |>
  dplyr::group_by(patient, cell_type, group) |>
  dplyr::summarise(count = n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(group_label = ifelse(group == 1, "CLR", "DII")) |>
  tidyplots::tidyplot(x = patient, y = count, color = cell_type) |>
  tidyplots::add_barstack_relative() |>
  tidyplots::adjust_x_axis_title("Patients") |>
  tidyplots::adjust_y_axis_title("Cell type proportion") |>
  tidyplots::adjust_font(fontsize=5) |>
  tidyplots::remove_legend_title() |>
  tidyplots::adjust_colors(new_colors = crc_tma_palette) |>
  tidyplots::remove_x_axis_labels() |>
  tidyplots::remove_x_axis_ticks() |>
  tidyplots::split_plot(by="group_label") |>
  tidyplots::save_plot("output/crc_tma/cell_proportions.png", units = "in", height = 1.5, width = 6) |>
  tidyplots::save_plot("output/crc_tma/cell_proportions.pdf", units = "in",  height = 9, width = 6) |>
  tidyplots::save_plot("output/crc_tma/cell_proportions.svg", units = "in", height = 9, width = 6.5)

codex_df |>
  dplyr::mutate(patient = as.factor(patient)) |>
  dplyr::group_by(patient, cell_type, group) |>
  dplyr::summarise(count = n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(group_label = ifelse(group == 1, "CLR", "DII")) |>
  tidyplots::tidyplot(x = patient, y = count, color = cell_type) |>
  tidyplots::add_barstack_relative() |>
  tidyplots::adjust_x_axis_title("Patients") |>
  tidyplots::adjust_y_axis_title("Cell type proportion") |>
  tidyplots::adjust_font(fontsize=5) |>
  tidyplots::remove_legend_title() |>
  tidyplots::adjust_colors(new_colors = crc_tma_palette) |>
  tidyplots::remove_x_axis_labels() |>
  tidyplots::remove_x_axis_ticks() |>
  tidyplots::remove_legend() |>
  # tidyplots::adjust_size(width = 1, unit = "in") |>
  tidyplots::split_plot(by="group_label") |>
  tidyplots::save_plot("output/crc_tma/cell_proportions_wo_legend.png", units = "in", height = 1.5, width = 6) |>
  tidyplots::save_plot("output/crc_tma/cell_proportions_wo_legend.pdf", units = "in",  height = 1.5, width = 6) |>
  tidyplots::save_plot("output/crc_tma/cell_proportions_wo_legend.svg", units = "in", height = 9, width = 6.5)

group1_codex_df <- codex_df |>
  filter(group == 1)
group1_patient_ids <- unique(group1_codex_df$patient)
group1_spomics <- list()
for(pat in group1_patient_ids) {
  group1_spomics[[paste0("patient", pat)]] <- group1_codex_df |>
    filter(patient == pat) |>
    mutate(sample = paste0("patient", pat)) |>
    create_spomic()
}

group2_codex_df <- codex_df |>
  filter(group == 2)
group2_patient_ids <- unique(group2_codex_df$patient)
group2_spomics <- list()
for(pat in group2_patient_ids) {
  group2_spomics[[paste0("patient", pat)]] <- group2_codex_df |>
    filter(patient == pat) |>
    mutate(sample = paste0("patient", pat)) |>
    create_spomic()
}

## ---------------------------

## process the spomics and compute spatial statistics

set.seed(123)
alert <- TRUE
# Lcross (global envelope)
# # Outputs are saved, no need to run again unless you want a different colocalization statistic
group1_spomics <- pblapply(seq_along(group1_spomics), function(i) {
  spomic <- group1_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
})
saveRDS(group1_spomics, "output/crc_tma/group1_spomics_Lcross_global_envelope.rds")
if(alert) beep()

group2_spomics <- pblapply(seq_along(group2_spomics), function(i) {
  spomic <- group2_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
})
saveRDS(group2_spomics, "output/crc_tma/group2_spomics_Lcross_global_envelope.rds")
if(alert) beep()

## Lcross.inhom
set.seed(123)
alert <- TRUE
# Lcross (global envelope)
# # Outputs are saved, no need to run again unless you want a different colocalization statistic
group1_spomics <- pblapply(seq_along(group1_spomics), function(i) {
  spomic <- group1_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross.inhom", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
}, cl=parallel::detectCores()-1)
saveRDS(group1_spomics, "output/crc_tma/group1_spomics_Lcross_inhom_global_envelope.rds")
if(alert) beep()

group2_spomics <- pblapply(seq_along(group2_spomics), function(i) {
  spomic <- group2_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross.inhom", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
}, cl=parallel::detectCores()-1)
saveRDS(group2_spomics, "output/crc_tma/group2_spomics_Lcross_inhom_global_envelope.rds")

## CLQ
# # Outputs are saved, no need to run again unless you want a different colocalization statistic
group1_spomics <- pblapply(seq_along(group1_spomics), function(i) {
  spomic <- group1_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "CLQ", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
}, cl=parallel::detectCores()-1)
saveRDS(group1_spomics, "output/crc_tma/group1_spomics_CLQ.rds")

group2_spomics <- pblapply(seq_along(group2_spomics), function(i) {
  spomic <- group2_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "CLQ", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
}, cl=parallel::detectCores()-1)
saveRDS(group2_spomics, "output/crc_tma/group2_spomics_CLQ.rds")

## ---------------------------

## run DC-SPOMIC

group1_spomics <- readRDS(file.path(output_path, "crc_tma/group1_spomics_Lcross_global_envelope.rds"))
group2_spomics <- readRDS(file.path(output_path, "crc_tma/group2_spomics_Lcross_global_envelope.rds"))

dcspomic <- create_dcspomic(
  group1_name = "CLR",
  group1_spomics = group1_spomics,
  group2_name = "DII",
  group2_spomics = group2_spomics)
dcspomic <- set_dcspomic_hyperparameters(
  dcspomic,
  # r = 198.939, # 75um
  r = 75,
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


###########################
generate_forest_plots <- function(dcspomic, i, j) {
  get_forest_data <- function(model, group_name) {
    yi <- as.vector(model$yi)
    vi <- model$vi
    slabs <- model$slab
    # slabs <- model$slab[-10]
    # slabs <- model$slab[-c(15,16)]
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

    # scale_size_continuous(range = c(0.5, 3)) +
    scale_size_continuous(range = c(0.5, 5)) +

    scale_y_continuous(
      breaks = df$y_pos,
      labels = df$study,
      expand = expansion(mult = c(0.01, 0.01))
    ) +

    labs(
      x = "Colocalization (Lcross: r=75 um)",
      # x = paste0("Colocalization (", checkmate_obj@details$hyperparameters$colocalization_type,", r=", checkmate_obj@details$hyperparameters$r, ")"),
      y = NULL,
      color = "Group",
      fill = "Group",
      size = "Weight",
      title = paste(pair, "Forest Plot")
    ) +

    ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        # axis.text = element_text(size = 10),
        legend.text = element_text(size=5),

        legend.title = element_text(size = 5),
        title = element_text(size = 7),
        legend.position = "right")


}


#############
i = "(b_cells)"
j = "(b_cells)"
generate_forest_plots(dcspomic, i, j)
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.pdf"), height = 5, width = 4, units = "in")
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.png"), height = 5, width = 4, units = "in")
ggsave(file.path(plot_path, "crc_tma/lcross_75um/bcells_bcells_forestplot.svg"), height = 5, width = 4, units = "in")

# ggsave("output/crc_tma/lcross_75um/bcells_bcells_forestplot.pdf", height = 5, width = 4, units = "in")
# ggsave("output/crc_tma/lcross_75um/bcells_bcells_forestplot.png", height = 5, width = 4, units = "in")
# ggsave("output/crc_tma/lcross_75um/bcells_bcells_forestplot.svg", height = 5, width = 4, units = "in")
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

## Colocalization network analysis
library(igraph)
library(tidyr)
library(dplyr)

# Assuming your dataframe is named results_without_unknowns

# 1. Split i_j into 'row' and 'col'
results_sep <- results_without_unknowns %>%
  separate(i_j, into = c("row", "col"), sep = "\\)_\\(", remove = FALSE) %>%
  mutate(
    row = gsub("\\(", "", row),
    col = gsub("\\)", "", col)
  )

# 2. Collect all unique labels
all_labels <- sort(unique(c(results_sep$row, results_sep$col)))

# 3. Create empty square matrix
adj_matrix <- matrix(NA, nrow = length(all_labels), ncol = length(all_labels),
                     dimnames = list(all_labels, all_labels))

# 4. Fill with z_score values
for(i in seq_len(nrow(results_sep))) {
  r <- results_sep$row[i]
  c <- results_sep$col[i]
  adj_matrix[r, c] <- results_sep$z_score[i]
}

# 5. Convert to data.frame if you wish
adj_matrix_df <- as.data.frame(adj_matrix)
adj_matrix_neg <- adj_matrix
adj_matrix_neg[adj_matrix_neg > 0 | is.na(adj_matrix_neg)] <- 0
adj_matrix_neg <- abs(adj_matrix_neg)

adj_matrix_df <- as.data.frame(adj_matrix_neg)
sym_adj_matrix <- (adj_matrix_neg + t(adj_matrix_neg)) / 2



colocalization_matrix <- results_without_unknowns

# Convert to graph from adjacency matrix
g <- graph_from_adjacency_matrix(
  sym_adj_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Optional: customize edge width by weight
E(g)$width <- E(g)$weight * 2  # scale for visibility

# Optional: set node size by degree or other metric
V(g)$size <- degree(g, mode = "in") + 5

# Optional: add color
E(g)$color <- ifelse(E(g)$weight > 1.96, "red", "gray")

# Plot the graph
plot(g, edge.arrow.size = 0.5, layout = layout_with_fr)

# Run Louvain clustering
cl <- cluster_louvain(g, resolution = 1.25)

# Assign cluster to each node
V(g)$group <- membership(cl)

# Plot with clusters
plot(g, vertex.color = V(g)$group, vertex.label = V(g)$name, layout = layout_with_fr)

library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)

# Convert igraph to tidygraph
tg <- as_tbl_graph(g)

# Create ggraph plot with label repulsion
set.seed(123)
ggraph(tg, layout = "fr") +
  # Lighter, semi-transparent edges weighted by value
  geom_edge_link(aes(width = weight), color = "#88888833") +  # RGBA hex for light gray, alpha=33
  scale_edge_width(range = c(0.2, 2.5)) +  # Adjust min/max thickness as desired

  # Nodes and labels
  geom_node_point(aes(color = as.factor(group)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +

  # A clean, minimal theme
  theme_void() +
  guides(edge_width = "none")  # Remove edge width legend

ggsave("output/crc_tma/lcross_75um/clr_network.png", plot = last_plot())
ggsave("output/crc_tma/lcross_75um/clr_network.pdf", plot = last_plot())
ggsave("output/crc_tma/lcross_75um/clr_network.svg", plot = last_plot())


three_colors <- tidyplots::new_color_scheme(x=c("#d7d2cb", "#0371b2", "#d36027"))
spomic <- dcspomic@group1_spomics[[1]]
spomic@df |>
  dplyr::mutate(cell_type2 = ifelse(cell_type == i, i, ifelse(cell_type == j, j, "(other)"))) |>
  dplyr::mutate(cell_type2 = factor(cell_type2, levels = c("(other)", i, j))) |>
  dplyr::arrange(cell_type2) |>
  tidyplots::tidyplot(x = x, y = y, color = cell_type2) |>
  tidyplots::add_data_points(size = 0.1) |>
  tidyplots::adjust_colors(new_colors = three_colors) |>
  tidyplots::remove_padding() |>
  tidyplots::adjust_title(spomic@details$sample) |>
  tidyplots::adjust_legend_title("Cell type") |>
  # tidyplots::add_caption(paste(count_i, i, "\n", count_j, j)) |>
  spomic_style()

# Define your 4 main cell types
cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")

# 5-color scheme: replace/add your preferred hex codes
five_colors <- tidyplots::new_color_scheme(
  x = c("#d7d2cb", "#0371b2", "#d36027", "#17b872", "#bb27d3")
)

# Map cell types in the data frame
spomic <- dcspomic@group1_spomics[[1]]
plot_quad <- function(spomic){
spomic@df |>
  dplyr::mutate(
    cell_type2 = dplyr::case_when(
      cell_type == "(b_cells)" ~ "(b_cells)",
      cell_type == "(tregs)" ~ "(tregs)",
      cell_type == "(cd8_t_cells)" ~ "(cd8_t_cells)",
      cell_type == "(cd4_t_cells_cd45ro)" ~ "(cd4_t_cells_cd45ro)",
      TRUE ~ "(other)"
    )
  ) |>
  dplyr::mutate(
    cell_type2 = factor(cell_type2, levels = c("(other)", "(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)"))
  ) |>
  dplyr::arrange(cell_type2) |>
  tidyplots::tidyplot(x = x, y = y, color = cell_type2) |>
  tidyplots::add_data_points(size = 0.1) |>
  tidyplots::adjust_colors(new_colors = five_colors) |>
  tidyplots::remove_padding() |>
  tidyplots::adjust_title(spomic@details$sample) |>
  tidyplots::adjust_legend_title("Cell type") |>
  # tidyplots::add_caption(
  #   paste(count_i, i, "\n", count_j, j, "\n", count_k, k, "\n", count_l, l)
  # ) |>
  spomic_style()
  }

plot_quad(dcspomic@group1_spomics[[3]]) |>
  tidyplots::save_plot(filename = "output/crc_tma/tls_cluster_group1.pdf")

plot_quad(dcspomic@group2_spomics[[3]]) |>
  tidyplots::save_plot(filename = "output/crc_tma/tls_cluster_group2.pdf")

dcspomic@group1_spomics[[1]]@df %>% filter(cell_type %in% c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)"))


library(dplyr)
library(RANN)

cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")
df_sub <- dcspomic@group1_spomics[[1]]@df %>% filter(cell_type %in% cell_types)

# Initialize result dataframe
result <- data.frame(
  cell_id = df_sub$cell_id,
  cell_type = df_sub$cell_type,
  nn_b_cells = NA,
  d_b_cells = NA,
  nn_tregs = NA,
  d_tregs = NA,
  nn_cd8 = NA,
  d_cd8 = NA,
  nn_cd4 = NA,
  d_cd4 = NA,
  dist_sum = NA
)

for (i in 1:nrow(df_sub)) {
  this_cell <- df_sub[i, ]
  this_xy <- as.numeric(this_cell[, c("x", "y")])
  dists <- c()
  ids <- c()

  for (target_type in cell_types) {
    if (this_cell$cell_type == target_type) {
      # Find nearest neighbor of *same* type, excluding self
      others <- df_sub %>% filter(cell_type == target_type, cell_id != this_cell$cell_id)
    } else {
      others <- df_sub %>% filter(cell_type == target_type)
    }
    # KD-tree is overkill for 1 query, use base distance
    if (nrow(others) == 0) {
      ids <- c(ids, NA)
      dists <- c(dists, NA)
    } else {
      xy_mat <- as.matrix(others[, c("x", "y")])
      diff <- sweep(xy_mat, 2, this_xy, "-")
      euclid <- sqrt(rowSums(diff^2))
      min_idx <- which.min(euclid)
      ids <- c(ids, others$cell_id[min_idx])
      dists <- c(dists, euclid[min_idx])
    }
  }
  # Assign results
  result$nn_b_cells[i] <- ids[1]
  result$d_b_cells[i] <- dists[1]
  result$nn_tregs[i] <- ids[2]
  result$d_tregs[i] <- dists[2]
  result$nn_cd8[i] <- ids[3]
  result$d_cd8[i] <- dists[3]
  result$nn_cd4[i] <- ids[4]
  result$d_cd4[i] <- dists[4]
  result$dist_sum[i] <- sum(dists, na.rm = TRUE)
}

result %>% group_by(cell_type) %>% summarise(mean(dist_sum, na.rm = TRUE))
mean(result$dist_sum, na.rm=TRUE)
hist(result$dist_sum, na.rm=TRUE)

grab_immune_cluster_dist <- function(spomic) {
  cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")
  df_sub <- spomic@df %>% filter(cell_type %in% cell_types)

  # Initialize result dataframe
  result <- data.frame(
    cell_id = df_sub$cell_id,
    cell_type = df_sub$cell_type,
    nn_b_cells = NA,
    d_b_cells = NA,
    nn_tregs = NA,
    d_tregs = NA,
    nn_cd8 = NA,
    d_cd8 = NA,
    nn_cd4 = NA,
    d_cd4 = NA,
    dist_sum = NA
  )

  for (i in 1:nrow(df_sub)) {
    this_cell <- df_sub[i, ]
    this_xy <- as.numeric(this_cell[, c("x", "y")])
    dists <- c()
    ids <- c()

    for (target_type in cell_types) {
      if (this_cell$cell_type == target_type) {
        # Find nearest neighbor of *same* type, excluding self
        others <- df_sub %>% filter(cell_type == target_type, cell_id != this_cell$cell_id)
      } else {
        others <- df_sub %>% filter(cell_type == target_type)
      }
      # KD-tree is overkill for 1 query, use base distance
      if (nrow(others) == 0) {
        ids <- c(ids, NA)
        dists <- c(dists, NA)
      } else {
        xy_mat <- as.matrix(others[, c("x", "y")])
        diff <- sweep(xy_mat, 2, this_xy, "-")
        euclid <- sqrt(rowSums(diff^2))
        min_idx <- which.min(euclid)
        ids <- c(ids, others$cell_id[min_idx])
        dists <- c(dists, euclid[min_idx])
      }
    }
    # Assign results
    result$nn_b_cells[i] <- ids[1]
    result$d_b_cells[i] <- dists[1]
    result$nn_tregs[i] <- ids[2]
    result$d_tregs[i] <- dists[2]
    result$nn_cd8[i] <- ids[3]
    result$d_cd8[i] <- dists[3]
    result$nn_cd4[i] <- ids[4]
    result$d_cd4[i] <- dists[4]
    result$dist_sum[i] <- sum(dists, na.rm = TRUE)
    result$dist_avg[i] <- mean(dists, na.rm = TRUE)

  }

  # mean(result$dist_sum, na.rm=TRUE)
  mean(result$dist_avg, na.rm=TRUE)

}

g1_immune_dist <- c()
for(i in 1:length(dcspomic@group1_spomics)) {
  print(i)
  g1_immune_dist[i] <- grab_immune_cluster_dist(dcspomic@group1_spomics[[i]])
}

g2_immune_dist <- c()
for(i in 1:length(dcspomic@group2_spomics)) {
  print(i)
  g2_immune_dist[i] <- grab_immune_cluster_dist(dcspomic@group2_spomics[[i]])
}

wilcox.test(g1_immune_dist, g2_immune_dist)
### ===========================

## Compare to t-test
retrieve_colocalization_stats <- function(dcspomic) {
  list1 <- list()
  for(i in 1:length(dcspomic@group1_spomics)) {
    list1[[i]] <- dcspomic@group1_spomics[[i]] |> get_spatial_summary() |> mutate(group = 1)
  }

  list2 <- list()
  for(i in 1:length(dcspomic@group2_spomics)) {
    list2[[i]] <- dcspomic@group2_spomics[[i]] |> get_spatial_summary() |> mutate(group = 2)
  }

  df <- rbind(bind_rows(list1), bind_rows(list2))
  return(df)
}

fix_dcspomic_names <- function(dcspomic) {
  for(i in 1:length(dcspomic@group1_spomics)) {
    names(dcspomic@group1_spomics)[i] <- dcspomic@group1_spomics[[i]]@details$sample
  }
  for(i in 1:length(dcspomic@group2_spomics)) {
    names(dcspomic@group2_spomics)[i] <- dcspomic@group2_spomics[[i]]@details$sample
  }
}

colocalization_stats <- retrieve_colocalization_stats(dcspomic)
head(colocalization_stats)

g1_df <- colocalization_stats |> filter(group == 1)
g2_df <- colocalization_stats |> filter(group == 2)

g1_samples <- unique(g1_df$sample)
g2_samples <- unique(g2_df$sample)

set.seed(567)
g1_samp_order <- sample(g1_samples, replace = FALSE)
g2_samp_order <- sample(g2_samples, replace = FALSE)

dcspomic_num_significant <- c()
num_significant <- c()
sig_list <- list()
dcspomic_sig_list <- list()
for(i in 2:17) {
  print(i)
  g1_df_filtered <- g1_df |> filter(sample %in% g1_samp_order[1:i])
  g2_df_filtered <- g2_df |> filter(sample %in% g2_samp_order[1:i])

  cell_pairs <- unique(c(g1_df_filtered$i_j, g2_df_filtered$i_j))
  cell_pairs <- cell_pairs[!grepl("\\(undefined\\)", cell_pairs)] # remove the pairs containing undefined cells

  ttest_res_list <- list()

  temp_dcspomic <- create_dcspomic(group1_name = "CLR", group2_name = "DII",
                                   group1_spomics = Filter(Negate(is.null), dcspomic@group1_spomics[g1_samp_order[1:i]]),
                                   group2_spomics = Filter(Negate(is.null), dcspomic@group2_spomics[g2_samp_order[1:i]]))
  temp_dcspomic <- set_dcspomic_hyperparameters(temp_dcspomic,
                                           r = 198.939, # 75um
                                           colocalization_type = "Lcross",
                                           tau_estimator = "SJ")
  temp_dcspomic <- run_dcspomic(temp_dcspomic)
  results_without_unknowns <- temp_dcspomic@results$differential_testing[
    !grepl("\\(undefined\\)", temp_dcspomic@results$differential_testing$i_j),
  ]
  results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
  dcspomic_sigpairs <- results_without_unknowns |> filter(FDR < 0.05) |> pull(i_j)
  dcspomic_num_significant <- c(dcspomic_num_significant, length(dcspomic_sigpairs))

  for(p in cell_pairs) {
    ttest_res <- t.test(x = g1_df_filtered |> filter(i_j == p) |> pull(colocalization_stat),
                        y = g2_df_filtered |> filter(i_j == p) |> pull(colocalization_stat))
    ttest_res_list[[p]] <- data.frame(i_j = p,
                                      t_stat = ttest_res$statistic,
                                      pval = ttest_res$p.value)
  }
  ttest_res_df <- bind_rows(ttest_res_list) |>
    arrange(pval) |>
    mutate(FDR = p.adjust(pval, method = "fdr"))
  stat_sig_pairs <- ttest_res_df |> filter(FDR < 0.05) |> pull(i_j)
  print(stat_sig_pairs)
  num_significant <- c(num_significant, length(stat_sig_pairs))
  sig_list[[i]] <- stat_sig_pairs
  dcspomic_sig_list[[i]] <- dcspomic_sigpairs


}

data.frame(num_samples_per_group = 2:17,
           num_samples_total = 2*2:17,
           num_significant = num_significant) |>
  tidyplot(x = num_samples_total, y = num_significant) |>
  add_data_points() |>
  add_line()

data.frame(num_samples_per_group = 2:17,
           num_samples_total = 2*2:17,
           num_significant = num_significant,
           dcspomic_num_significant = dcspomic_num_significant) |>
  tidyplot(x = num_samples_total, y = dcspomic_num_significant) |>
  add_data_points() |>
  add_line()

## survival

crc_metadata <- read.csv("data/crc_tma/tma_metadata.csv") |>
  janitor::clean_names()
head(crc_metadata)

fix_dcspomic_names <- function(dcspomic) {
  for(i in 1:length(dcspomic@group1_spomics)) {
    names(dcspomic@group1_spomics)[i] <- dcspomic@group1_spomics[[i]]@details$sample
  }
  for(i in 1:length(dcspomic@group2_spomics)) {
    names(dcspomic@group2_spomics)[i] <- dcspomic@group2_spomics[[i]]@details$sample
  }
}

retrieve_colocalization_stats <- function(dcspomic) {
  list1 <- list()
  for(i in 1:length(dcspomic@group1_spomics)) {
    list1[[i]] <- dcspomic@group1_spomics[[i]] |> get_spatial_summary()
  }

  list2 <- list()
  for(i in 1:length(dcspomic@group2_spomics)) {
    list2[[i]] <- dcspomic@group2_spomics[[i]] |> get_spatial_summary()
  }

  df <- rbind(bind_rows(list1), bind_rows(list2))
  return(df)
}

df <- retrieve_colocalization_stats(dcspomic)
df_wide <- df %>%
  select(sample, i_j, colocalization_stat) %>%
  pivot_wider(
    names_from = i_j,
    values_from = colocalization_stat
  )
df_wide$patient <- as.numeric(sub("patient", "", df_wide$sample))
# normalize <- function(x) (x - min(x)) / (max(x) - min(x))
normalize <- function(x, ...) {
  # normalization code, e.g.:
  (x - min(x)) / (max(x) - min(x))
}

df_normalized <- df_wide
df_normalized[ , -1] <- lapply(df_wide[ , -1], normalize)

# View result
print(df_wide)

library(survival)

survival_df <- inner_join(df_wide, crc_metadata, by = c("patient"="patient"))
survival_df <- inner_join(df_normalized, crc_metadata, by = c("patient"="patient"))

# Fit Cox proportional hazards model
cox_model <- coxph(Surv(dfs, dfs_censor) ~ age +
                     # sex +
                     # cp_tnm_simple +
                     group,
                   data = survival_df)
summary(cox_model)


spatial_cox_model <- coxph(Surv(dfs, dfs_censor) ~ age +
                     # sex +
                     # group +
                       # cp_tnm_simple +
                       # p_tnm +
                     # p_t +
                     # p_n +
                     `(b_cells)_(b_cells)` +
                     `(cd4_t_cells_cd45ro)_(cd4_t_cells_cd45ro)` +
                     `(cd4_t_cells_cd45ro)_(tregs)`,
                   data = survival_df)

spatial_cox_model <- coxph(Surv(dfs, dfs_censor) ~ age +
                             # sex +
                             # group +
                             # cp_tnm_simple +
                             # p_tnm +
                             # p_t +
                             # p_n +
                             `(b_cells)_(b_cells)`
                           # +
                           #   `(cd4_t_cells_cd45ro)_(cd4_t_cells_cd45ro)` +
                           #   `(cd4_t_cells_cd45ro)_(tregs)`
                           ,
                           data = survival_df)

# Show the results
summary(spatial_cox_model)





