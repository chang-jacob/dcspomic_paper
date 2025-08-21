## ---------------------------
##
## Script name: plotting_utils.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-28
##
## ---------------------------
##
#' Description:
#' This script contains useful functions and color schemes for creating plots
#' for the DC-SPOMIC manuscript
##
## ---------------------------

## load libraries

library(tidyplots)
library(ggplot2)

## ---------------------------

## color scheme for simulations 
simulation_colors <- tidyplots::new_color_scheme(
  x = c("(A)" = "#56B4E9",
        "(B)" = "#E69F00",
        "(C)" = "#D3D3D3",
        "(Rare 1)" = "#009E73",
        "(Rare 2)" = "#D55E00"),
  name = "simulation_color_scheme")

## color scheme for the CRC TMA

crc_tma_palette <- tidyplots::new_color_scheme(c(
  # Immune
  "(b_cells)"="#CF6B97",
  "(cd4_t_cells_cd45ro)"="#56BAE9",
  "(cd68_cd163_macrophages)"="#1171B8",
  "(cd68_macrophages)"="#2596D3",
  "(cd8_t_cells)"="#12B283",
  "(granulocytes)"="#38A1B6",
  "(plasma_cells)"="#6174A4",

  "(tregs)"="#89B64A",

  # Stroma / Endo
  "(smooth_muscle)"="#F6C414",
  "(stroma)"="#F3AE16",
  "(vasculature)"="#F09C1C",

  # Malignant
  "(tumor_cells)"="#D44F1C",

  # Other
  "(undefined)"="#D3D3D3"
),
name = "crc_tma_palette"
)


