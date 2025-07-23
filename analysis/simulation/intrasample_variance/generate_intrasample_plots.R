# devtools::load_all("/home/groups/plevriti/jachang4/spomic")
devtools::load_all("/Users/jacobchang/Lab/spomic")

library(tidyplots)
library(dplyr)
library(tidyr)

oak_path <- "/oak/stanford/groups/plevriti"
simulation_colors <- new_color_scheme(x = c("(A)" = "#56B4E9",
                                            "(B)" = "#E69F00",
                                            "(C)" = "#D3D3D3",
                                            "(Rare 1)" = "#009E73",
                                            "(Rare 2)" = "#D55E00"),
                                      name = "simulation_color_scheme")

homogeneous_spomics <- readRDS("output/simulation/intrasample_variance/homogeneous_pattern_spomics.rds")

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
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.svg")

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

divergent_spomics <- readRDS("output/simulation/intrasample_variance/divergent_pattern_spomics.rds")

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

cluster_spomics <- readRDS("output/simulation/intrasample_variance/cluster_pattern_spomics.rds")

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

homogeneous_stats <- readRDS("output/simulation/intrasample_variance/homogeneous_pattern_stats.rds")

plot_simulation_results(homogeneous_stats, "") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.png") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.pdf") |>
  save_plot("output/simulation/intrasample_variance/homogeneous_process/intrasample_variance_comparison.svg")

divergent_stats <- readRDS("output/simulation/intrasample_variance/divergent_pattern_stats.rds")

plot_simulation_results(divergent_stats, "") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.png") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.pdf") |>
  save_plot("output/simulation/intrasample_variance/divergent_process/intrasample_variance_comparison.svg")

cluster_stats <- readRDS("output/simulation/intrasample_variance/cluster_pattern_stats.rds")

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

