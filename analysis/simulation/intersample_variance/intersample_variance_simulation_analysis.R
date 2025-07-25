## ---------------------------
##
## Script name: intersample_variance_simulation_analysis.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-07
##
## ---------------------------
##
#' Description:
#' This script reads in the doubly stochastic point process simulation data
#' and calculates the RMSE between DC-SPOMIC estimates and ground truth estimates
#' for both mu and tau2.
#'
#' Note:
#' The input files to this script are generated through the script
#' cluster_sbatch.R. This script was run on the Sherlock HPC at Stanford.
##
## ---------------------------

## Load libraries

devtools::load_all("/Users/jacobchang/Lab/spomic")
library(spatstat)
library(dplyr)
library(tidyr)
library(purrr)
library(pbmcapply)
library(metafor)
library(pbapply)
library(concaveman)

## Set hyperparameters
N_SIM <- 100
WINDOW_SIZE <- 2000
R <- 200
BASE_DIR <- "/Users/jacobchang/Lab/dcspomic_paper"
data_path <- "output/simulation/intersample_variance"

## ---------------------------

## Read simulation data files and join into a single data frame
## The simulations were split into four files for faster processing

df1 <- readRDS(file.path(data_path, "cluster_process_df1.rds"))
df2 <- readRDS(file.path(data_path, "cluster_process_df2.rds"))
df3 <- readRDS(file.path(data_path, "cluster_process_df3.rds"))
df4 <- readRDS(file.path(data_path, "cluster_process_df4.rds"))
results <- rbind(df1, df2, df3, df4)

# results <- readRDS(file.path(data_path, "cluster_process_test.rds"))
results$i_j <- gsub("([A-Z])_([A-Z])", "(\\1)_(\\2)", results$i_j)


## Aggregate simulation results at the "patient" level

intermediate_level <- results |>
  dplyr::group_by(patient, i_j) |>
  dplyr::summarise(patient_mean = mean(colocalization_stat, na.rm = TRUE),
                   patient_sigma2 = var(colocalization_stat, na.rm = TRUE))

## Aggregate to get population estimates of mu ad tau2

population_colocalization <- intermediate_level |>
  dplyr::group_by(i_j) |>
  dplyr::summarise(mu = mean(patient_mean, na.rm = TRUE),
                   tau2 = var(patient_mean, na.rm = TRUE))

## ---------------------------

## Function for RMSE

get_rmse <- function(truth, estimates) {
  idx <- !is.na(truth) & !is.na(estimates)
  sqrt(mean((estimates[idx] - truth[idx])^2))
}
# get_rmse <- function(truth, estimates) {
#   return(sqrt(sum((estimates-truth)^2, na.rm=TRUE) / length(estimates)))
# }

get_rmse_plots <- function(results, pair, seed = 123) {
  n_inter <- 100
  n_intra <- 100
  results_ab <- results |> dplyr::filter(i_j == pair)

  vec <- 1:n_inter
  half <- length(vec) / 2
  idxs <- as.vector(rbind(vec[1:half], vec[length(vec):(half + 1)]))

  # This sampling order is used downstream so that for every number of
  # samples that we give DC-SPOMIC access to, they are as different as possible
  # This essentially makes it as difficult as possible to assess tau2.
  patient_sample_order <- results_ab |>
    dplyr::group_by(patient) |>
    dplyr::summarise(v = var(colocalization_stat, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(v)) |>
    dplyr::mutate(id = dplyr::row_number()) |>
    dplyr::select(patient, id)

  results_ab <- dplyr::inner_join(results_ab, patient_sample_order)

  population_colocalization_ab <- population_colocalization |> dplyr::filter(i_j == pair)

  set.seed(seed)
  random_intra_sample <- sample(1:100, replace = FALSE)

  mu_hats <- c()
  tau2_hats <- c()
  naive_mu_hats <- c()
  naive_tau2_hats <- c()

  df_list <- list()

  counter <- 1
  for(k in random_intra_sample){
    print(counter)
    counter <- counter + 1
    for(i in 2:n_inter) {
      reduced_result <- results_ab |>
        dplyr::filter(id %in% idxs[1:i]) |>
        dplyr::filter(intra_sample == k)

      model <- metafor::rma.uni(yi = colocalization_stat,
                                vi = colocalization_var,
                                method = "SJ",
                                data = reduced_result)

      mu_hats[i] <- coef(summary(model))$estimate
      tau2_hats[i] <- summary(model)$tau2
      naive_mu_hats[i] <- mean(reduced_result$colocalization_stat, na.rm = TRUE)
      naive_tau2_hats[i] <- var(reduced_result$colocalization_stat, na.rm = TRUE)

    }

    foo <- data.frame(n = 1:n_inter,
                      mu_hat = mu_hats,
                      naive_mu_hat = naive_mu_hats,
                      true_mu = population_colocalization_ab$mu,
                      tau2_hat = tau2_hats,
                      true_tau2 = population_colocalization_ab$tau2,
                      naive_tau2_hat = naive_tau2_hats,
                      intra_n = k)
    df_list[[k]] <- foo
  }

  dcspomic_df <- dplyr::bind_rows(df_list)
  dcspomic_df <- dcspomic_df |>
    dplyr::mutate(n = as.factor(n)) |>
    dplyr::group_by(n) |>
    dplyr::summarise(
      mu_hat_var = var(mu_hat, na.rm=TRUE),
      mu_hat = mean(mu_hat, na.rm=TRUE),

      naive_mu_hat_var = var(naive_mu_hat, na.rm = TRUE),
      naive_mu_hat = mean(naive_mu_hat, na.rm = TRUE),

      true_mu_var = var(true_mu, na.rm = TRUE),
      true_mu = mean(true_mu, na.rm = TRUE),

      tau2_hat_var = var(tau2_hat, na.rm = TRUE),
      tau2_hat = mean(tau2_hat, na.rm = TRUE),

      true_tau2_var = var(true_tau2, na.rm = TRUE),
      true_tau2 = mean(true_tau2, na.rm = TRUE),

      naive_tau2_hat_var = var(naive_tau2_hat, na.rm = TRUE),
      naive_tau2_hat = mean(naive_tau2_hat, na.rm = TRUE))

  mu_hat_rmses <- c()
  naive_mu_hat_rmses <- c()
  tau2_hat_rmses <- c()
  naive_tau2_hat_rmses <- c()

  for(i in 2:n_inter) {
    foo2 <- dcspomic_df |> dplyr::filter(n %in% 1:i)
    mu_hat_rmses[i] <- get_rmse(truth=foo2$true_mu, estimates=foo2$mu_hat)
    naive_mu_hat_rmses[i] <- get_rmse(truth=foo2$true_mu, estimates=foo2$naive_mu_hat)
    tau2_hat_rmses[i] <- get_rmse(truth=foo2$true_tau2, estimates=foo2$tau2_hat)
    naive_tau2_hat_rmses[i] <- get_rmse(truth=foo2$true_tau2, estimates=foo2$naive_tau2_hat)
  }
  df <- dcspomic_df |>
    dplyr::mutate(mu_hat_rmse = mu_hat_rmses,
                  naive_mu_hat_rmse = naive_mu_hat_rmses,
                  tau2_hat_rmse = tau2_hat_rmses,
                  naive_tau2_hat_rmse = naive_tau2_hat_rmses)

  df_mu_long <- df |>
    dplyr::select(n, mu_hat_rmse, naive_mu_hat_rmse) |>
    dplyr::rename("DC-SPOMIC" = "mu_hat_rmse", "naive" = "naive_mu_hat_rmse") |>
    tidyr::pivot_longer(cols = c("DC-SPOMIC", naive), names_to = "Method", values_to = "RMSE")

  df_tau2_long <- df |>
    dplyr::select(n, tau2_hat_rmse, naive_tau2_hat_rmse) |>
    dplyr::rename("DC-SPOMIC" = "tau2_hat_rmse", "naive" = "naive_tau2_hat_rmse") |>
    tidyr::pivot_longer(cols = c("DC-SPOMIC", naive), names_to = "Method", values_to = "RMSE")

  upper_limit <- max(c(df_mu_long$RMSE, df_tau2_long$RMSE), na.rm=TRUE)

  df_mu_long |>
    dplyr::mutate(n = as.numeric(n)) |>
    tidyplots::tidyplot(x = n, y = RMSE, color = Method) |>
    tidyplots::add_data_points(size = 0.25) |>
    tidyplots::add_line() |>
    tidyplots::adjust_font(fontsize=5) |>
    tidyplots::adjust_title(paste("Comparison of", pair, "colocalization estimates for \u03BC"), fontsize = 7) |>
    tidyplots::adjust_x_axis_title("Number of samples") |>
    tidyplots::adjust_y_axis_title(expression("RMSE of \u03BC")) |>
    tidyplots::remove_legend_title() |>
    tidyplots::adjust_legend_position(position = "bottom") |>
    tidyplots::adjust_size(width = 3, height = 2, unit = "in") |>
    tidyplots::adjust_y_axis(limits = c(0, sqrt(upper_limit))) |>
    tidyplots::adjust_x_axis(cut_short_scale = TRUE) |>
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_mu_rmse_plot.png")) |>
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_mu_rmse_plot.pdf")) |>
    # tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_mu_rmse_plot.svg"))
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_mu_rmse_plot.svg"))


  df_tau2_long |>
    dplyr::mutate(n = as.numeric(n)) |>
    tidyplots::tidyplot(x = n, y = RMSE, color = Method) |>
    tidyplots::add_data_points(size = 0.25) |>
    tidyplots::add_line() |>
    tidyplots::adjust_font(fontsize=5) |>
    tidyplots::adjust_title(paste("Comparison of", pair, "colocalization estimates for \u03C4\u00B2"), fontsize = 7) |>
    tidyplots::adjust_x_axis_title("Number of samples") |>
    tidyplots::adjust_y_axis_title(expression("RMSE of \u03C4\u00B2")) |>
    tidyplots::remove_legend_title() |>
    tidyplots::adjust_legend_position(position = "bottom") |>
    tidyplots::adjust_size(width = 3, height = 2, unit = "in") |>
    tidyplots::adjust_y_axis(limits = c(0, upper_limit)) |>
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_tau2_rmse_plot.png")) |>
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_tau2_rmse_plot.pdf")) |>
    # tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_tau2_rmse_plot.svg"))
    tidyplots::save_plot(filename = paste0("output/simulation/intersample_variance/", pair, "_tau2_rmse_plot.svg"))

}

## ---------------------------

## Generate RMSE plots for all cell pairs
get_rmse_plots(results, pair = "(A)_(A)", seed = 123)
get_rmse_plots(results, pair = "(A)_(B)", seed = 123)
get_rmse_plots(results, pair = "(B)_(A)", seed = 123)
get_rmse_plots(results, pair = "(B)_(B)", seed = 123)

## ---------------------------

library(tidyplots)

## Visualize intra and inter- sample examples
spomic_list <- readRDS("spomics_first10x10.rds")
length(spomic_list)

simulation_colors <-
  new_color_scheme(c("A" = "#56B4E9",
                     "B" = "#E69F00"),
                   name = "simulation_color_scheme")


## intra-sample

spomic_list[[1]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intrasample_examples/example1.svg")

spomic_list[[2]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intrasample_examples/example2.svg")

spomic_list[[3]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intrasample_examples/example3.svg")

spomic_list[[4]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intrasample_examples/example4.svg")

## inter-sample

spomic_list[[20]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intersample_examples/example1.svg")

spomic_list[[60]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intersample_examples/example2.svg")

spomic_list[[40]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intersample_examples/example3.svg")

spomic_list[[90]] %>%
  plot_spomic() %>%
  adjust_colors(new_colors = simulation_colors) %>%
  save_plot("output/simulation/intersample_variance/intersample_examples/example4.svg")


