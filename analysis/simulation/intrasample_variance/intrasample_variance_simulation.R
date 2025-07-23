## ---------------------------
##
## Script name: intrasample_variance_simulation.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-22
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
library(tidyr)
library(progress)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)
library(tidyplots)
library(pbapply)
library(purrr)
library(pbapply)
library(parallel)

## Set global parameters
square_size <- 2000
window <- owin(square(square_size))
cell_types <- c("(A)", "(B)", "(C)", "(Rare 1)", "(Rare 2)")

n_rounds <- 10
n_sims <- 100

oak_path <- "/oak/stanford/groups/plevriti"

# LAMBDA_A and LAMBDA_B are different for the three different patterns, but
# LAMBDA_C, LAMBDA_RARE1, and LAMBDA_RARE2 are constant throughout the script.
LAMBDA_C <- 0.1/1000
LAMBDA_RARE1 <- 0.01/1000
LAMBDA_RARE2 <- 0.005/1000

simulation_colors <- new_color_scheme(x = c("(A)" = "#56B4E9",
                                            "(B)" = "#E69F00",
                                            "(C)" = "#D3D3D3",
                                            "(Rare 1)" = "#009E73",
                                            "(Rare 2)" = "#D55E00"),
                                      name = "simulation_color_scheme")

## ---------------------------

safe_saveRDS <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = file)
}

## Define function to generate a point pattern and compute spatial statistics 

run_sim <- function(i, j, lambda, win, types = c("(A)", "(B)", "(C)", "(Rare 1)", "(Rare 2)")) {
  intensity <- lambda
  pp <- rmpoispp(lambda = intensity, win = win, types = types)
  df <- pp |> as.data.frame() |> rename(cell_type = marks) |> mutate(sample = "")
  spomic <- create_spomic(df)
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = FALSE, colocalization_type="Lcross")
  spomic <- get_spatial_stats(spomic)
  result <- get_spatial_summary(spomic) |> mutate(round = j, sim = i)
  
  list(summary = result, spomic = if(i <= 10) spomic else NULL)
}

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

## ---------------------------

## Homogeneous spatial pattern (aka complete spatial randomness (CSR))

LAMBDA_A <- 0.2/1000
LAMBDA_B <- 0.15/1000

set.seed(123)

lambda <- c(LAMBDA_A, LAMBDA_B, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2)

homogeneous_list <- vector("list", n_rounds)

for (j in 1:n_rounds) {
  cat("Running round", j, "\n")

  sims <- pblapply(1:n_sims, function(i) {
    run_sim(i, j, lambda, window, cell_types)
  }, cl = parallel::detectCores() - 1)  # Use all but one core
  
  homogeneous_spomic_colocalization_results <- map(sims, "summary")
  homogeneous_spomics <- map(sims[1:10], "spomic")
  
  homogeneous_list[[j]] <- bind_rows(homogeneous_spomic_colocalization_results)
}
homogeneous_rounds <- bind_rows(homogeneous_list)

ground_truth_estimates <- homogeneous_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(homogeneous_rounds)

safe_saveRDS(object = df, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "homogeneous_pattern_stats.rds"))

safe_saveRDS(object = homogeneous_spomics, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "homogeneous_pattern_spomics.rds"))

# plot_simulation_results(df, title = "Homogeneous pattern") |>
#   save_plot(filename = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.png") %>% 
#     save_plot(filename = "output/simulation/intrasample_variance/homogeneous_process_sigma2k.pdf")

## ---------------------------

## Diverging spatial pattern (aka complete spatial randomness (CSR))

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
  
  sims <- pblapply(1:n_sims, function(i) {
    run_sim(i, j, lambda, window, cell_types)
  }, cl = parallel::detectCores() - 1)  # Use all but one core
  
  divergent_spomic_colocalization_results <- map(sims, "summary")
  divergent_spomics <- map(sims[1:10], "spomic")
  
  divergent_list[[j]] <- bind_rows(divergent_spomic_colocalization_results)
}

divergent_rounds <- bind_rows(divergent_list)

ground_truth_estimates <- divergent_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(divergent_rounds)

safe_saveRDS(object = df, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "divergent_pattern_stats.rds"))

safe_saveRDS(object = divergent_spomics, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "divergent_pattern_spomics.rds"))

# 
# plot_simulation_results(df, title = "Divergent pattern") %>%  
#   save_plot(filename = "output/simulation/intrasample_variance/divergent_process_sigma2k.png") %>%
#   save_plot(filename = "output/simulation/intrasample_variance/divergent_process_sigma2k.pdf")


## ---------------------------

## Cluster spatial pattern

run_cluster_sim <- function(
    i, round, parent_lambda, offspring_n, offspring_dispersion, window, 
    square_size, cell_types, LAMBDA_C, LAMBDA_RARE1, LAMBDA_RARE2
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
      spomic = if(i <= 10) spomic else NULL
    )
  }, error = function(e) {
    print(sprintf("SIM %d ERROR: %s", i, conditionMessage(e)))
    NULL
  })
}

set.seed(123)
parent_lambda0 <- 0.00001 + runif(1, 0, 0.0001)
offspring_n0 <- 3 + rnbinom(1, size=10, mu=2)
offspring_dispersion0 <- 175


all_sims <- vector("list", n_rounds)
all_spomics <- vector("list", n_rounds)

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
      window = window,
      square_size = square_size,
      cell_types = cell_types,
      LAMBDA_C = LAMBDA_C,
      LAMBDA_RARE1 = LAMBDA_RARE1,
      LAMBDA_RARE2 = LAMBDA_RARE2
    ),
    cl = parallel::detectCores() - 1
  )
  all_sims[[j]] <- lapply(sims, function(x) x$summary)
  all_spomics[[j]] <- lapply(sims, function(x) x$spomic)
}

summaries_list <- unlist(all_sims, recursive = FALSE)
summaries_list <- summaries_list[!sapply(summaries_list, is.null)]
cluster_rounds <- bind_rows(summaries_list)

spomics_list <- unlist(all_spomics, recursive = FALSE)
cluster_spomics <- spomics_list[!sapply(spomics_list, is.null)]

ground_truth_estimates <- cluster_rounds |>
  group_by(round, i_j) |>
  summarise(theta_k_hat = mean(colocalization_stat, na.rm=TRUE),
            sigma2_k = var(colocalization_stat, na.rm=TRUE))

df <- ground_truth_estimates |> inner_join(cluster_rounds)

safe_saveRDS(object = df, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "cluster_pattern_stats.rds"))

safe_saveRDS(object = cluster_spomics, 
             file = file.path(oak_path, 
                              "jachang4", 
                              "dcspomic_output", 
                              "intrasample_simulation", 
                              "cluster_pattern_spomics.rds"))


# 
# 
# plot_simulation_results(df, title = "Cluster pattern") %>%  
#   save_plot(filename = "output/simulation/intrasample_variance/cluster_process_sigma2k.png") %>%
#   save_plot(filename = "output/simulation/intrasample_variance/clyuster_process_sigma2k.pdf")
