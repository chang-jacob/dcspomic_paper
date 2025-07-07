## ---------------------------
##
## Script name: intrasample_variance_simulation.R
##
## Author: Jake Chang
##
## Date Modified: 2025-06-20
##
## ---------------------------
##
#' Description:
#' This script simulates three different spatial patterns from spatial point
#' processes. Multiple realizations of each pattern allow for empirical calculation
#' of intra-sample variance, $\sigma_k^2$. This can then be compared to intra-sample
#' variance calculated through spatial bootstrapping when only a single point
#' process realization is available.
##
## ---------------------------

## Load libraries
devtools::load_all("/Users/jacobchang/Lab/spomic")
library(spatstat)
library(dplyr)
library(progress)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)

## Set global parameters
square_size <- 2000
window <- owin(square(square_size))
cell_types <- c("(A)", "(B)", "(C)", "(Rare 1)", "(Rare 2)")

# LAMBDA_A and LAMBDA_B are different for the three different patterns, but
# LAMBDA_C, LAMBDA_RARE1, and LAMBDA_RARE2 are constant throughout the script.
LAMBDA_C <- 0.1/1000
LAMBDA_RARE1 <- 0.01/1000
LAMBDA_RARE2 <- 0.005/1000

## ---------------------------

## Define function to generate consistent plots and save them
generate_variance_plot <- function(spomic_colocalization_df, spomic1, title, output_name) {
  # spomic_colocalization_df <- bind_rows(cluster_spomic_colocalization_results)
  ground_truth_estimates <- spomic_colocalization_df |>
    group_by(i_j) |>
    summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
              sigma2_k = var(colocalization_stat, na.rm=TRUE))

  df <- ground_truth_estimates |> inner_join(spomic1 |> get_spatial_summary())
  corr <- cor(df$colocalization_var, df$sigma2_k, use = "complete.obs")

  df |>
    tidyplot(x = colocalization_var, y = sigma2_k) |>
    add_data_points() |>
    add_data_labels_repel(label = i_j, max.overlaps = 10, fontsize = 5) |>
    add_curve_fit(method = "lm", alpha = 0.1) |>
    add_caption(paste("r =", round(corr, 2))) |>
    adjust_x_axis_title("$hat(sigma)[k]^2$") |>
    adjust_y_axis_title("$sigma[k]^2$") |>
    adjust_font(fontsize=5) |>
    adjust_title(title = title, fontsize = 7)

  save_plot(plot = ggplot2::last_plot(),
            filename = output_name)
}

simulation_colors <-
  new_color_scheme(c("(A)" = "#56B4E9",
                     "(B)" = "#E69F00",
                     "(C)" = "#D3D3D3",
                     "(Rare 1)" = "#009E73",
                     "(Rare 2)" = "#D55E00"),
                   name = "simulation_color_scheme")

## ---------------------------

## Homogeneous spatial pattern (aka complete spatial randomness (CSR))
LAMBDA_A <- 0.2/1000
LAMBDA_B <- 0.15/1000

set.seed(123)
homogeneous_spomic <- list()
homogeneous_intensity_plots <- list()
homogeneous_spomic_colocalization_results <- list()
for(i in 1:100) {
  print(i)
  homogeneous_intensity <- c(LAMBDA_A,
                             LAMBDA_B,
                             LAMBDA_C,
                             LAMBDA_RARE1,
                             LAMBDA_RARE2)
  pp <- rmpoispp(lambda = homogeneous_intensity, win = window, types = cell_types)
  df <- pp |> as.data.frame() |> rename(cell_type = marks) |> mutate(sample = "")
  spomic <- create_spomic(df)
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, colocalization_type="Lcross")
  spomic <- get_spatial_stats(spomic)
  homogeneous_spomic_colocalization_results[[i]] <- get_spatial_summary(spomic)

  if(i %in% 1:10) homogeneous_spomic[[i]] <- spomic

  homogeneous_intensity_plots[[i]] <- plot_spomic(spomic)
}
saveRDS(homogeneous_spomic_colocalization_results, "output/simulation/intrasample_variance/homogeneous_simulation.rds")
saveRDS(homogeneous_spomic, "output/simulation/intrasample_variance/homogeneous_spomics.rds")

homogeneous_spomic_colocalization <- bind_rows(homogeneous_spomic_colocalization_results)
generate_variance_plot(spomic_colocalization_df = homogeneous_spomic_colocalization,
                       spomic1 = homogeneous_spomic[[1]],
                       title = "Homogeneous Pattern",
                       output_name = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.png")
generate_variance_plot(spomic_colocalization_df = homogeneous_spomic_colocalization,
                       spomic1 = homogeneous_spomic[[1]],
                       title = "Homogeneous Pattern",
                       output_name = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.pdf")

homogeneous_intensity_plots[[1]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous1.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous1.pdf")

homogeneous_intensity_plots[[2]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous2.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous2.pdf")

homogeneous_intensity_plots[[3]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous3.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous3.pdf")

homogeneous_intensity_plots[[4]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous4.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/homogeneous4.pdf")

## ---------------------------

## Diverging spatial pattern (aka complete spatial randomness (CSR))

LAMBDA_A <- function(x, y, rate) { dexp(x, rate = rate) }
LAMBDA_B <- function(x, y, rate) {
  max_value <- 2000
  dexp(max_value - x, rate = rate)
}


set.seed(123)
divergent_spomic <- list()
divergent_intensity_plots <- list()
divergent_spomic_colocalization_results <- list()
rate <- runif(1, min=1e-10, max = 0.005)
for(i in 1:100) {
  print(i)
  # Generate each component point process
  pp1 <- rpoispp(lambda = function(x,y) LAMBDA_A(x, y, rate), win = window)
  pp1$marks <- factor(rep("(A)", pp1$n), levels = cell_types)

  pp2 <- rpoispp(lambda = function(x,y) LAMBDA_B(x, y, rate), win = window)
  pp2$marks <- factor(rep("(B)", pp2$n), levels = cell_types)

  pp3 <- rpoispp(lambda = LAMBDA_C, win = window)
  pp3$marks <- factor(rep("(C)", pp3$n), levels = cell_types)

  pp4 <- rpoispp(lambda = LAMBDA_RARE1, win = window)
  pp4$marks <- factor(rep("(Rare 1)", pp4$n), levels = cell_types)

  pp5 <- rpoispp(lambda = LAMBDA_RARE2, win = window)
  pp5$marks <- factor(rep("(Rare 2)", pp5$n), levels = cell_types)

  # Combine all types
  pp_all <- superimpose(pp1, pp2, pp3, pp4, pp5, W = window)
  pp_all$marks <- factor(pp_all$marks, levels = cell_types)

  df <- as.data.frame(pp_all) %>%
    rename(cell_type = marks) %>%
    mutate(sample = paste0("sample_", i))

  spomic <- create_spomic(df)
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = TRUE, colocalization_type="Lcross")
  spomic <- get_spatial_stats(spomic)


  divergent_spomic_colocalization_results[[i]] <- get_spatial_summary(spomic)

  if(i %in% 1:10) divergent_spomic[[i]] <- spomic

  divergent_intensity_plots[[i]] <- plot_spomic(spomic)
}

saveRDS(divergent_spomic_colocalization_results, "output/simulation/intrasample_variance/divergent_simulation.rds")
saveRDS(divergent_spomic, "output/simulation/intrasample_variance/divergent_spomics.rds")

divergent_spomic_colocalization <- bind_rows(divergent_spomic_colocalization_results)
generate_variance_plot(spomic_colocalization_df = divergent_spomic_colocalization,
                       spomic1 = divergent_spomic[[1]],
                       title = "Diverging Pattern",
                       output_name = "output/simulation/intrasample_variance/divergent_process_sigma2k.png")
generate_variance_plot(spomic_colocalization_df = divergent_spomic_colocalization,
                       spomic1 = divergent_spomic[[1]],
                       title = "Diverging Pattern",
                       output_name = "output/simulation/intrasample_variance/divergent_process_sigma2k.pdf")

divergent_intensity_plots[[1]] |>
  adjust_colors(new_colors = simulation_colors) |>
  remove_title()
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent1.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent1.pdf")

divergent_intensity_plots[[2]] |>
  adjust_colors(new_colors = simulation_colors) |>
  remove_title()
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent2.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent2.pdf")

divergent_intensity_plots[[3]] |>
  adjust_colors(new_colors = simulation_colors) |>
  remove_title()
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent3.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent3.pdf")

divergent_intensity_plots[[4]] |>
  adjust_colors(new_colors = simulation_colors) |>
  remove_title()
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent4.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/divergent4.pdf")

## ---------------------------

## Cluster spatial pattern

set.seed(123)
cluster_spomic <- list()
cluster_intensity_plots <- list()
cluster_spomic_colocalization_results <- list()

parent_lambda0 = 0.00001 + runif(1, 0, 0.0001)
offspring_n0 = 3 + rnbinom(1, size=10, mu=2)
offspring_dispersion0 = 175

for(i in 1:100) {
print(i)
    parent_lambda <- parent_lambda0
    offspring_n <- offspring_n0
    offspring_dispersion <- offspring_dispersion0

    parents <- rpoispp(lambda = parent_lambda, win = square(window_size))
    offspring_x <- numeric()
    offspring_y <- numeric()

    if (parents$n > 0) {
      num_offspring <- rpois(parents$n, offspring_n) + 1
      new_x <- rnorm(sum(num_offspring),
                     mean = rep(parents$x, num_offspring),
                     sd = offspring_dispersion/2)
      new_y <- rnorm(sum(num_offspring),
                     mean = rep(parents$y, num_offspring),
                     sd = offspring_dispersion/2)

      valid_idx <- (new_x >= 0) & (new_x <= window_size) & (new_y >= 0) & (new_y <= window_size)
      offspring_x <- new_x[valid_idx]
      offspring_y <- new_y[valid_idx]
    }

    parent_ppp <- ppp(parents$x, parents$y, window = square(window_size), marks = rep("(A)", parents$n))
    offspring_ppp <- ppp(offspring_x, offspring_y, window = square(window_size), marks = rep("(B)", length(offspring_x)))

    pp3 <- rpoispp(lambda = LAMBDA_C, win = window)
    pp3$marks <- factor(rep("(C)", pp3$n), levels = cell_types)

    pp4 <- rpoispp(lambda = LAMBDA_RARE1, win = window)
    pp4$marks <- factor(rep("(Rare 1)", pp4$n), levels = cell_types)

    pp5 <- rpoispp(lambda = LAMBDA_RARE2, win = window)
    pp5$marks <- factor(rep("(Rare 2)", pp5$n), levels = cell_types)

    marked_process <- superimpose(parent_ppp, offspring_ppp, pp3, pp4, pp5)
    marks(marked_process) <- as.factor(marks(marked_process))
    spomic <- create_spomic(ppp_to_df(marked_process, sample_name = ""))
    spomic <- set_spomic_hyperparameters(spomic = spomic,
                                         colocalization_type = "Lcross",
                                         fixed_distance = FALSE,
                                         r = R)
    spomic <- get_spatial_stats(spomic)

   cluster_spomic_colocalization_results[[i]] <- get_spatial_summary(spomic)
   if(i %in% 1:10) cluster_spomic[[i]] <- spomic
   cluster_intensity_plots[[i]] <- plot_spomic(spomic)
}

saveRDS(cluster_spomic_colocalization_results, "output/simulation/intrasample_variance/cluster_simulation.rds")
saveRDS(cluster_spomic, "output/simulation/intrasample_variance/cluster_spomics.rds")

cluster_spomic_colocalization <- bind_rows(cluster_spomic_colocalization_results)
generate_variance_plot(spomic_colocalization_df = cluster_spomic_colocalization,
                       spomic1 = cluster_spomic[[1]],
                       title = "Cluster Pattern",
                       output_name = "output/simulation/intrasample_variance/cluster_process_sigma2k.png")
generate_variance_plot(spomic_colocalization_df = cluster_spomic_colocalization,
                       spomic1 = cluster_spomic[[1]],
                       title = "Cluster Pattern",
                       output_name = "output/simulation/intrasample_variance/cluster_process_sigma2k.pdf")

cluster_intensity_plots[[1]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster1.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster1.pdf")

cluster_intensity_plots[[2]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster2.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster2.pdf")

cluster_intensity_plots[[3]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster3.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster3.pdf")

cluster_intensity_plots[[4]] |>
  adjust_colors(new_colors = simulation_colors)
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster4.png")
save_plot(plot = ggplot2::last_plot(),
          filename = "output/simulation/intrasample_variance/cluster4.pdf")
