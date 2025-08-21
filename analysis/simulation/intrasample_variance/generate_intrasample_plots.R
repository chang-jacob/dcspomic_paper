## ---------------------------
##
## Script name: genarate_intrasample_plots.R
##
## Author: Jake Chang
##
## Date Modified: 2025-08-06
##
## ---------------------------
##
#' Description:
#' This script takes the outputs that are saved in intrasample_variance_simulation.R
#' and generates plots that are used in the figures of the manuscript. 
#' 
#' Note:
#' The files from intrasample_variance_simulation.R are saved to SCRATCH and/or OAK 
#' with the date that the script was run. You will need to adjust the input path 
#' to adjust for this.
##
## ---------------------------

library(spomic)
library(tidyplots)
library(dplyr)
library(tidyr)

## setup 
scratch <- config::get("scratch")
oak <- config::get("oak")
date <- format(Sys.Date(), "%Y%m%d")
# date <- "20250728"

source("analysis/plotting_utils.R")

homogeneous_spomics <- readRDS(file.path(scratch, 
                                         "dcspomic_output", 
                                         date,
                                         "intrasample_simulation", 
                                         "homogeneous_pattern_spomics.rds"))

plot_spomic(homogeneous_spomics[[1]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.pdf") 
# |>
#   save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.svg")

plot_spomic(homogeneous_spomics[[2]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example2.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example2.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example2.svg")

plot_spomic(homogeneous_spomics[[3]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example3.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example3.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example3.svg")

plot_spomic(homogeneous_spomics[[4]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example4.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example4.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example4.svg")

# ---------------

divergent_spomics <- readRDS(file.path(oak, 
                                         "dcspomic_output", 
                                         date,
                                         "intrasample_simulation", 
                                         "divergent_pattern_spomics.rds"))

plot_spomic(divergent_spomics[[1]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example1.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example1.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example1.svg")

plot_spomic(divergent_spomics[[2]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example2.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example2.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example2.svg")

plot_spomic(divergent_spomics[[3]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example3.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example3.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example3.svg")

plot_spomic(divergent_spomics[[4]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example4.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example4.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/example4.svg")

# ---------------

cluster_spomics <- readRDS(file.path(oak, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "cluster_pattern_spomics.rds"))

plot_spomic(cluster_spomics[[1]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example1.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example1.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example1.svg")

plot_spomic(cluster_spomics[[2]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example2.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example2.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example2.svg")

plot_spomic(cluster_spomics[[3]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example3.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example3.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example3.svg")

plot_spomic(cluster_spomics[[4]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks() |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example4.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example4.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/example4.svg")

# ---------------
homogeneous_stats <- readRDS(file.path(scratch, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "homogeneous_pattern_stats.rds"))
divergent_stats <- readRDS(file.path(scratch, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "divergent_pattern_stats.rds"))
cluster_stats <- readRDS(file.path(scratch, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "cluster_pattern_stats.rds"))

homogeneous_corr <- homogeneous_stats %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>% 
  mutate(pattern = "Homogeneous")
divergent_corr <- divergent_stats %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>% 
  mutate(pattern = "Divergent")
cluster_corr <- cluster_stats %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>% 
  mutate(pattern = "Cluster")

corr_df <- rbind(homogeneous_corr, divergent_corr, cluster_corr)

corr_df %>% 
  mutate(pattern = factor(pattern, levels = c("Cluster", "Divergent", "Homogeneous"))) %>% 
  tidyplot(y = pattern, x = correlation2, color = pattern) %>% 
  add_violin(draw_quantiles = 0.5) %>% 
  add_median_value(accuracy = 0.01) %>% 
  add_reference_lines(x = 0) %>% 
  adjust_legend_position(position = "bottom") %>% 
  remove_legend_title() %>% 
  adjust_size(width=3, height = 2, unit = "in") %>% 
  # remove_x_axis_title() %>% 
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations.png") |>
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations.pdf") |>
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations.svg")

## repeat analysis without rare cell types

homogeneous_corr_wo_rare <- homogeneous_stats %>% 
  filter(!grepl("Rare", i_j)) %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>% 
  mutate(pattern = "Homogeneous")
divergent_corr_wo_rare <- divergent_stats %>% 
  filter(!grepl("Rare", i_j)) %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>%   mutate(pattern = "Divergent")
cluster_corr_wo_rare <- cluster_stats %>% 
  filter(!grepl("Rare", i_j)) %>% 
  group_by(round, sim) %>% 
  summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ), 
            correlation2 = cor(theta_k_hat, colocalization_stat, use = "complete.obs" )) %>%   mutate(pattern = "Cluster")

corr_df_wo_rare <- rbind(homogeneous_corr_wo_rare , divergent_corr_wo_rare , cluster_corr_wo_rare )

corr_df_wo_rare %>% 
  mutate(pattern = factor(pattern, levels = c("Cluster", "Divergent", "Homogeneous"))) %>% 
  tidyplot(y = pattern, x = correlation, color = pattern) %>% 
  add_violin(draw_quantiles = 0.5) %>% 
  add_median_value(accuracy = 0.01) %>% 
  add_reference_lines(x = 0) %>% 
  adjust_legend_position(position = "bottom") %>% 
  remove_legend_title() %>% 
  adjust_size(width=3, height = 2, unit = "in") %>% 
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations_wo_rare.png") |>
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations_wo_rare.pdf") |>
  save_plot("output/simulation/intrasample_variance/intrasample_variance_correlations_wo_rare.svg")

plot_simulation_results <- function(df, title) {
  # Reduce to first slice per (i_j, round)
  foo <- df |>
    group_by(i_j, round) |>
    slice(1) |>
    ungroup()

  corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

  foo_processed <- foo %>%
    filter(is.finite(colocalization_var), is.finite(sigma2_k)) %>%
    separate(i_j, into = c("i", "j"), sep = "\\)_\\(", remove = FALSE) %>%
    mutate(
      i = paste0(i, ")"),
      j = paste0("(", j),
      prevalence = case_when(
        grepl("Rare", i, ignore.case = TRUE) & grepl("Rare", j, ignore.case = TRUE) ~ "both rare",
        grepl("Rare", i, ignore.case = TRUE) | grepl("Rare", j, ignore.case = TRUE) ~ "common and rare",
        TRUE ~ "both common"
      )
    )

  color_map <- c(
    "both rare" = "#d62728",
    "common and rare" = "#1f77b4",
    "both common" = "#2ca02c"
  )

  foo_processed <- foo_processed %>%
    mutate(color_code = color_map[prevalence])

  p <- foo_processed |>
    tidyplot(x = colocalization_var, y = sigma2_k) |>
    add_data_points(color = "gray80", size = 0.1) |>
    add_curve_fit(method = "lm", se = FALSE, color = "black", alpha = 0.6)

  p |>
    add_data_points(color = foo_processed$color_code, size = 1) |>
    add_caption(paste("r =", round(corr, 2))) |>
    adjust_x_axis_title("$hat(sigma)[k]^2$") |>
    adjust_y_axis_title("$sigma[k]^2$") |>
    adjust_font(fontsize = 7) |>
    adjust_title(title, fontsize = 7)
}

homogeneous_stats <- readRDS(file.path(scratch, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "homogeneous_pattern_stats.rds"))

homogeneous_stats %>% group_by(sim) %>% summarise(correlation = cor(sigma))

foo <- homogeneous_stats %>% filter(sim == 1) 
foo <- homogeneous_stats |>
  group_by(i_j, sim) |>
  slice(1) |>
  ungroup()

foo <- homogeneous_stats %>% filter(round == 1)
foo2 <- foo %>% group_by(sim) %>% summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ))
corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")


foo <- homogeneous_stats %>% group_by(round, sim) %>% summarise(correlation = cor(sigma2_k, colocalization_var, use = "complete.obs" ))
foo %>% 
  mutate(type = "homogeneous") %>% 
  tidyplot(x = type, y = correlation) %>% 
  add_violin(draw_quantiles = 0.5)
# %>% 
#   add_data_points_beeswarm()
  
  ggplot(aes(y=correlation)) + 
  geom_boxplot() + 
  theme_minimal()


plot_simulation_results(homogeneous_stats, "") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.svg")

divergent_stats <- readRDS(file.path(oak, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "divergent_pattern_stats.rds"))

plot_simulation_results(divergent_stats, "") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.svg")

cluster_stats <- readRDS(file.path(oak, 
                                       "dcspomic_output", 
                                       date,
                                       "intrasample_simulation", 
                                       "cluster_pattern_stats.rds"))

plot_simulation_results(cluster_stats, "") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison.svg")

# ------------------------------------------------------------------------------
# generate version of the variance plots with legend for Illustrator file
df <- cluster_stats
foo <- df |>
  group_by(i_j, round) |>
  slice(1) |>
  ungroup()

corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

foo_processed <- foo %>%
  filter(is.finite(colocalization_var), is.finite(sigma2_k)) %>%
  separate(i_j, into = c("i", "j"), sep = "\\)_\\(", remove = FALSE) %>%
  mutate(
    i = paste0(i, ")"),
    j = paste0("(", j),
    prevalence = case_when(
      grepl("Rare", i, ignore.case = TRUE) & grepl("Rare", j, ignore.case = TRUE) ~ "both rare",
      grepl("Rare", i, ignore.case = TRUE) | grepl("Rare", j, ignore.case = TRUE) ~ "common and rare",
      TRUE ~ "both common"
    )
  )

color_map <- c(
  "both rare" = "#d62728",
  "common and rare" = "#1f77b4",
  "both common" = "#2ca02c"
)
prevalence_colors <- new_color_scheme(x = c(  "both rare" = "#d62728",
                                              "common and rare" = "#1f77b4",
                                              "both common" = "#2ca02c"),
                                      name = "prevalence_color_scheme")


foo_processed <- foo_processed %>%
  mutate(color_code = color_map[prevalence])

p <- foo_processed |>
  tidyplot(x = colocalization_var, y = sigma2_k, color = prevalence) |>
  add_data_points() |>
  adjust_colors(new_colors = prevalence_colors)

p |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison_w_legend.png") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison_w_legend.pdf") |>
  save_plot("output/simulation/intrasample_variance/cluster_process/intrasample_variance_comparison_w_legend.svg")

