## ---------------------------
##
## Script name: crc_data_process.R
##
## Author: Jake Chang
##
## Date Modified: 2025-08-20
##
## ---------------------------
##
#' Description:
#' This script takes the dataset from Schurch et al., with data-processing following
#' the spicyR paper. 
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
# Lcross (global envelope)
group1_spomics <- pblapply(seq_along(group1_spomics), function(i) {
  spomic <- group1_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
})
saveRDS(group1_spomics, file.path(oak, "dcspomic_output/crc_tma/group1_spomics_Lcross_global_envelope.rds"))

group2_spomics <- pblapply(seq_along(group2_spomics), function(i) {
  spomic <- group2_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  return(spomic)
})
saveRDS(group2_spomics, file.path(oak, "dcspomic_output/crc_tma/group2_spomics_Lcross_global_envelope.rds"))

## Lcross.inhom
set.seed(123)
# Lcross (global envelope)
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
# group1_spomics <- pblapply(seq_along(group1_spomics), function(i) {
#   spomic <- group1_spomics[[i]]
#   spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "CLQ", fixed_distance = FALSE)
#   spomic <- get_spatial_stats(spomic)
#   return(spomic)
# }, cl=parallel::detectCores()-1)
# saveRDS(group1_spomics, "output/crc_tma/group1_spomics_CLQ.rds")
# 
# group2_spomics <- pblapply(seq_along(group2_spomics), function(i) {
#   spomic <- group2_spomics[[i]]
#   spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "CLQ", fixed_distance = FALSE)
#   spomic <- get_spatial_stats(spomic)
#   return(spomic)
# }, cl=parallel::detectCores()-1)
# saveRDS(group2_spomics, "output/crc_tma/group2_spomics_CLQ.rds")
