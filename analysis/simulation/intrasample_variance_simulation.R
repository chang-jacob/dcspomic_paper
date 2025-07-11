## ---------------------------
##
## Script name: intrasample_variance_simulation.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-09
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
# devtools::load_all("/Users/jacobchang/Lab/spomic")
devtools::load_all("/home/groups/plevriti/jachang4/spomic")
library(spatstat)
library(dplyr)
library(progress)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)
library(tidyplots)

library(pbapply)
library(purrr)


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

run_sim <- function(i, j, lambda, win, types) {
  homogeneous_intensity <- lambda
  pp <- rmpoispp(lambda = homogeneous_intensity, win = win, types = types)
  df <- pp |> as.data.frame() |> rename(cell_type = marks) |> mutate(sample = "")
  spomic <- create_spomic(df)
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = FALSE, colocalization_type="Lcross")
  spomic <- get_spatial_stats(spomic)
  result <- get_spatial_summary(spomic) |> mutate(round = j, sim = i)
  plot_obj <- plot_spomic(spomic)
  # Only return what you want to keep. For now, return summary & optionally plots.
  list(summary = result, plot = if(i <= 10) plot_obj else NULL)
}

set.seed(123)
n_rounds <- 10
n_sims <- 100
lambda <- c(LAMBDA_A, LAMBDA_B, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2)

homogeneous_list <- vector("list", n_rounds)

for (j in 1:n_rounds) {
  cat("Running round", j, "\n")
  
  # Use pblapply for parallel simulation within each round
  sims <- pblapply(1:n_sims, function(i) {
    run_sim(i, j, lambda, window, cell_types)
  }, cl = parallel::detectCores() - 1)  # Use all but one core
  
  # Extract summaries and (if needed) plots
  homogeneous_spomic_colocalization_results <- map(sims, "summary")
  homogeneous_spomic_plots <- map(sims[1:10], "plot")
  
  homogeneous_list[[j]] <- bind_rows(homogeneous_spomic_colocalization_results)
  # Optionally: store plots as well if you need them:
  # homogeneous_plots_list[[j]] <- homogeneous_spomic_plots
}

homogeneous_rounds <- bind_rows(homogeneous_list)

ground_truth_estimates <- homogeneous_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(homogeneous_rounds)
# corr <- cor(df$colocalization_var, df$sigma2_k, use = "complete.obs")
foo <- df |>
  group_by(i_j, round) |>
  slice(1)
corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

foo %>% 
  # mutate(round = factor(round)) %>% 
  tidyplot(x = colocalization_var, y = sigma2_k) |>
  add_data_points() |>
  add_data_labels_repel(label = i_j, max.overlaps = 10, fontsize = 5) |>
  add_curve_fit(method = "lm", alpha = 0) |>
  add_caption(paste("r =", round(corr, 2))) |>
  adjust_x_axis_title("$hat(sigma)[k]^2$") |>
  adjust_y_axis_title("$sigma[k]^2$") |>
  adjust_font(fontsize=5) |>
  adjust_title(title = "Homogeneous pattern", fontsize = 7) %>% 
  save_plot(filename = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.png") %>% 
  save_plot(filename = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.pdf") %>% 
  save_plot(filename = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.svg") 




# homogeneous_spomic <- list()
# homogeneous_intensity_plots <- list()
# homogeneous_spomic_colocalization_results <- list()
homogeneous_list <- list()
for(j in 1:10){
  print(j)
  homogeneous_spomic <- list()
  homogeneous_intensity_plots <- list()
  homogeneous_spomic_colocalization_results <- list()
  for(i in 1:100) {
    # print(i)
    homogeneous_intensity <- c(LAMBDA_A,
                               LAMBDA_B,
                               LAMBDA_C,
                               LAMBDA_RARE1,
                               LAMBDA_RARE2)
    pp <- rmpoispp(lambda = homogeneous_intensity, win = window, types = cell_types)
    df <- pp |> as.data.frame() |> rename(cell_type = marks) |> mutate(sample = "")
    spomic <- create_spomic(df)
    spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = TRUE, colocalization_type="Lcross")
    spomic <- get_spatial_stats(spomic)
    homogeneous_spomic_colocalization_results[[i]] <- get_spatial_summary(spomic) |>
      mutate(round = j)

    if(i %in% 1:10) homogeneous_spomic[[i]] <- spomic

    homogeneous_intensity_plots[[i]] <- plot_spomic(spomic)
  }
  homogeneous_list[[j]] <- bind_rows(homogeneous_spomic_colocalization_results)
}
homogeneous_rounds <- bind_rows(homogeneous_list)
saveRDS(homogeneous_rounds, "output/simulation/intrasample_variance/homogeneous_simulation_10rounds.rds")

generate_variance_plot2 <- function(spomic_colocalization_df, spomic1, title, output_name) {
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
# 
#   save_plot(plot = ggplot2::last_plot(),
#             filename = output_name)
}

head(homogeneous_rounds)

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
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = TRUE, colocalization_type="Lcross")
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


n_rounds <- 10
n_sims <- 100



set.seed(123)
rate <- runif(1, min=1e-10, max = 0.005)

LAMBDA_A <- function(x, y) { dexp(x, rate =rate) }
LAMBDA_B <- function(x, y) {
  max_value <- 2000
  dexp(max_value - x, rate = rate)
}

lambda <- c(LAMBDA_A, LAMBDA_B, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2)

divergent_list <- vector("list", n_rounds)

for (j in 1:n_rounds) {
  cat("Running round", j, "\n")
  
  # Use pblapply for parallel simulation within each round
  sims <- pblapply(1:n_sims, function(i) {
    run_sim(i, j, lambda, window, cell_types)
  }, cl = parallel::detectCores() - 1)  # Use all but one core
  
  # Extract summaries and (if needed) plots
  divergent_spomic_colocalization_results <- map(sims, "summary")
  divergent_spomic_plots <- map(sims[1:10], "plot")
  
  divergent_list[[j]] <- bind_rows(divergent_spomic_colocalization_results)
  # Optionally: store plots as well if you need them:
  # homogeneous_plots_list[[j]] <- homogeneous_spomic_plots
}

divergent_rounds <- bind_rows(divergent_list)


ground_truth_estimates <- divergent_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(divergent_rounds)
# corr <- cor(df$colocalization_var, df$sigma2_k, use = "complete.obs")

foo <- df |>
  group_by(i_j, round) |>
  slice(1) %>% 
  # summarise(colocalization_var = mean(colocalization_var, na.rm=TRUE), 
  #           sigma2_k = mean(sigma2_k, na.rm=TRUE)) |> 
  ungroup() 

corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

foo %>% 
  # mutate(round = factor(round)) %>% 
  tidyplot(x = colocalization_var, y = sigma2_k) |>
  add_data_points() |>
  add_data_labels_repel(label = i_j, max.overlaps = 10, fontsize = 5) |>
  add_curve_fit(method = "lm", alpha = 0) |>
  add_caption(paste("r =", round(corr, 2))) |>
  adjust_x_axis_title("$hat(sigma)[k]^2$") |>
  adjust_y_axis_title("$sigma[k]^2$") |>
  adjust_font(fontsize=5) |>
  adjust_title(title = "Divergent pattern", fontsize = 7) %>%
  save_plot(filename = "output/simulation/intrasample_variance/divergent_process_sigma2k.png") %>%
  save_plot(filename = "output/simulation/intrasample_variance/divergent_process_sigma2k.pdf") %>%
  save_plot(filename = "output/simulation/intrasample_variance/divergent_process_sigma2k.svg")
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









#########


run_cluster_sim <- function(i, parent_lambda, offspring_n, offspring_dispersion, window_size, square_size, cell_types, square_size, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2, R) {
  library(spatstat)
  # Generate parents
  parents <- spatstat.random::rpoispp(lambda = parent_lambda, win = window_size)
  offspring_x <- numeric()
  offspring_y <- numeric()
  
  if (parents$n > 0) {
    num_offspring <- rpois(parents$n, offspring_n) + 1
    new_x <- rnorm(sum(num_offspring), mean = rep(parents$x, num_offspring), sd = offspring_dispersion/2)
    new_y <- rnorm(sum(num_offspring), mean = rep(parents$y, num_offspring), sd = offspring_dispersion/2)
    valid_idx <- (new_x >= 0) & (new_x <= square_size) & (new_y >= 0) & (new_y <= square_size)
    offspring_x <- new_x[valid_idx]
    offspring_y <- new_y[valid_idx]
  }
  
  parent_ppp <- spatstat.geom::ppp(parents$x, parents$y, window = window_size, marks = rep("(A)", parents$n))
  offspring_ppp <- spatstat.geom::ppp(offspring_x, offspring_y, window = window_size, marks = rep("(B)", length(offspring_x)))
  
  # Add other types
  pp3 <- spatstat.random::rpoispp(lambda = LAMBDA_C, win = window)
  pp3$marks <- factor(rep("(C)", pp3$n), levels = cell_types)
  pp4 <- spatstat.random::rpoispp(lambda = LAMBDA_RARE1, win = window)
  pp4$marks <- factor(rep("(Rare 1)", pp4$n), levels = cell_types)
  pp5 <- spatstat.random::rpoispp(lambda = LAMBDA_RARE2, win = window)
  pp5$marks <- factor(rep("(Rare 2)", pp5$n), levels = cell_types)
  
  marked_process <- spatstat.geom::superimpose(parent_ppp, offspring_ppp, pp3, pp4, pp5)
  spatstat.geom::marks(marked_process) <- as.factor(spatstat.geom::marks(marked_process))
  
  spomic <- create_spomic(ppp_to_df(marked_process, sample_name = ""))
  spomic <- set_spomic_hyperparameters(spomic = spomic, colocalization_type = "Lcross", fixed_distance = FALSE, r = 200)
  spomic <- get_spatial_stats(spomic)
  
  list(
    summary = get_spatial_summary(spomic),
    plot = plot_spomic(spomic)
  )
}

library(pbapply)
library(parallel)

set.seed(123)
parent_lambda0 <- 0.00001 + runif(1, 0, 0.0001)
offspring_n0 <- 3 + rnbinom(1, size=10, mu=2)
offspring_dispersion0 <- 175

# Adjust as needed: if you want these random per simulation, move into run_cluster_sim
n_sims <- 100

sims <- pblapply(
  1:n_sims,
  function(i) run_cluster_sim(
    i,
    parent_lambda = parent_lambda0,
    offspring_n = offspring_n0,
    offspring_dispersion = offspring_dispersion0,
    window_size = window_size,
    cell_types = cell_types,
    window = window,
    square_size = square_size,
    LAMBDA_C = LAMBDA_C,
    LAMBDA_RARE1 = LAMBDA_RARE1,
    LAMBDA_RARE2 = LAMBDA_RARE2,
    R = R
  ),
  cl = detectCores() - 1
)

# Results
cluster_spomic_colocalization_results <- lapply(sims, `[[`, "summary")
cluster_spomic_plots <- lapply(sims[1:10], `[[`, "plot")


##########
run_cluster_sim <- function(
    i, round, parent_lambda, offspring_n, offspring_dispersion, window_size, 
    square_size, cell_types, window, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2, R
) {
  suppressPackageStartupMessages({
    library(spatstat.geom)
    library(spatstat.random)
  })
  tryCatch({
    print(sprintf("SIM %d: starting...", i))
    # --- Parent process ---
    parents <- rpoispp(lambda = parent_lambda, win = window)
    offspring_x <- numeric()
    offspring_y <- numeric()
    if (parents$n > 0) {
      num_offspring <- rpois(parents$n, offspring_n) + 1
      new_x <- rnorm(sum(num_offspring), mean = rep(parents$x, num_offspring), sd = offspring_dispersion/2)
      new_y <- rnorm(sum(num_offspring), mean = rep(parents$y, num_offspring), sd = offspring_dispersion/2)
      valid_idx <- (new_x >= 0) & (new_x <= square_size) & (new_y >= 0) & (new_y <= square_size)
      offspring_x <- new_x[valid_idx]
      offspring_y <- new_y[valid_idx]
    }
    parent_ppp <- ppp(parents$x, parents$y, window = window, marks = rep("(A)", parents$n))
    offspring_ppp <- ppp(offspring_x, offspring_y, window = window, marks = rep("(B)", length(offspring_x)))
    # --- Add other types ---
    pp3 <- rpoispp(lambda = LAMBDA_C, win = window)
    pp3$marks <- factor(rep("(C)", pp3$n), levels = cell_types)
    pp4 <- rpoispp(lambda = LAMBDA_RARE1, win = window)
    pp4$marks <- factor(rep("(Rare 1)", pp4$n), levels = cell_types)
    pp5 <- rpoispp(lambda = LAMBDA_RARE2, win = window)
    pp5$marks <- factor(rep("(Rare 2)", pp5$n), levels = cell_types)
    marked_process <- superimpose(parent_ppp, offspring_ppp, pp3, pp4, pp5)
    marks(marked_process) <- as.factor(marks(marked_process))
    spomic <- create_spomic(ppp_to_df(marked_process, sample_name = ""))
    spomic <- set_spomic_hyperparameters(spomic = spomic, colocalization_type = "Lcross", fixed_distance = FALSE, r = 200)
    spomic <- get_spatial_stats(spomic)
    temp <- get_spatial_summary(spomic)
    temp$round = round
    list(
      summary = temp,
      plot = plot_spomic(spomic)
    )
  }, error = function(e) {
    print(sprintf("SIM %d ERROR: %s", i, conditionMessage(e)))
    NULL
  })
}


# Main simulation wrapper
library(pbapply)
library(parallel)

set.seed(123)
parent_lambda0 <- 0.00001 + runif(1, 0, 0.0001)
offspring_n0 <- 3 + rnbinom(1, size=10, mu=2)
offspring_dispersion0 <- 175
n_sims <- 100

# --- Ensure all required objects are defined ---
# window_size, square_size, cell_types, window, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2, R
# ... define those before running

# --- Run the parallel simulation safely ---
sims <- pbapply::pblapply(
  1:n_sims,
  function(i) run_cluster_sim(
    i,
    parent_lambda = parent_lambda0,
    offspring_n = offspring_n0,
    offspring_dispersion = offspring_dispersion0,
    window_size = window_size,
    square_size = square_size,
    cell_types = cell_types,
    window = window,
    LAMBDA_C = LAMBDA_C,
    LAMBDA_RARE1 = LAMBDA_RARE1,
    LAMBDA_RARE2 = LAMBDA_RARE2,
    R = R
  ),
  cl = parallel::detectCores() - 1
)

cluster_rounds <- bind_rows(lapply(sims, `[[`, "summary"))
ground_truth_estimates <- cluster_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(cluster_rounds)
# corr <- cor(df$colocalization_var, df$sigma2_k, use = "complete.obs")

foo <- df |>
  group_by(i_j, round) |>
  slice(1) %>% 
  # summarise(colocalization_var = mean(colocalization_var, na.rm=TRUE), 
  #           sigma2_k = mean(sigma2_k, na.rm=TRUE)) |> 
  ungroup() 

corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

foo %>% 
  # mutate(round = factor(round)) %>% 
  tidyplot(x = colocalization_var, y = sigma2_k) |>
  add_data_points() |>
  add_data_labels_repel(label = i_j, max.overlaps = 10, fontsize = 5) |>
  add_curve_fit(method = "lm", alpha = 0) |>
  add_caption(paste("r =", round(corr, 2))) |>
  adjust_x_axis_title("$hat(sigma)[k]^2$") |>
  adjust_y_axis_title("$sigma[k]^2$") |>
  adjust_font(fontsize=5) |>
  adjust_title(title = "Divergent pattern", fontsize = 7)


n_rounds <- 10
n_sims <- 100

all_sims <- vector("list", n_rounds)

for (j in 1:n_rounds) {
  cat("Running round", j, "\n")
  
  # ... set or randomize your parameters for each round if needed ...
  sims <- pbapply::pblapply(
    1:n_sims,
    function(i) run_cluster_sim(
      i = i,
      round = j,   # <--- Pass the round number here!
      parent_lambda = parent_lambda0,
      offspring_n = offspring_n0,
      offspring_dispersion = offspring_dispersion0,
      window_size = window_size,
      square_size = square_size,
      cell_types = cell_types,
      window = window,
      LAMBDA_C = LAMBDA_C,
      LAMBDA_RARE1 = LAMBDA_RARE1,
      LAMBDA_RARE2 = LAMBDA_RARE2,
      R = R
    ),
    cl = parallel::detectCores() - 1
  )
  all_sims[[j]] <- lapply(sims, function(x) x$summary)
}

# Combine all results into a single dataframe:
library(dplyr)
summaries_list <- unlist(all_sims, recursive = FALSE)
summaries_list <- summaries_list[!sapply(summaries_list, is.null)]
cluster_rounds <- bind_rows(summaries_list)


ground_truth_estimates <- cluster_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(cluster_rounds)
# corr <- cor(df$colocalization_var, df$sigma2_k, use = "complete.obs")

foo <- df |>
  group_by(i_j, round) |>
  slice(1) %>% 
  # summarise(colocalization_var = mean(colocalization_var, na.rm=TRUE), 
  #           sigma2_k = mean(sigma2_k, na.rm=TRUE)) |> 
  ungroup() 

corr <- cor(foo$colocalization_var, foo$sigma2_k, use = "complete.obs")

foo %>% 
  # mutate(round = factor(round)) %>% 
  tidyplot(x = colocalization_var, y = sigma2_k) |>
  add_data_points() |>
  add_data_labels_repel(label = i_j, max.overlaps = 10, fontsize = 5) |>
  add_curve_fit(method = "lm", alpha = 0) |>
  add_caption(paste("r =", round(corr, 2))) |>
  adjust_x_axis_title("$hat(sigma)[k]^2$") |>
  adjust_y_axis_title("$sigma[k]^2$") |>
  adjust_font(fontsize=5) |>
  adjust_title(title = "Cluster pattern", fontsize = 7) %>% 
  save_plot(filename = "output/simulation/intrasample_variance/cluster_process_sigma2k.png") %>%
  save_plot(filename = "output/simulation/intrasample_variance/cluster_process_sigma2k.pdf") %>%
  save_plot(filename = "output/simulation/intrasample_variance/cluster_process_sigma2k.svg")




