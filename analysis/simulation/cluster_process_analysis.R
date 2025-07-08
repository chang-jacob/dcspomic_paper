
# read cluster process outputs generated on Sherlock
# simulations were split into four runs to improve processing time
# data_path <- "output/simulation/cluster_process"
data_path <- "output/simulation/intersample_variance_v2"

df1 <- readRDS(file.path(data_path, "cluster_process_df1.rds"))
df2 <- readRDS(file.path(data_path, "cluster_process_df2.rds"))
df3 <- readRDS(file.path(data_path, "cluster_process_df3.rds"))
df4 <- readRDS(file.path(data_path, "cluster_process_df4.rds"))

df <- rbind(df1, df2, df3, df4)

# calculate metrics for ground truth and intra-sample variance
devtools::load_all("/Users/jacobchang/Lab/spomic")
library(spatstat)
library(dplyr)
library(tidyr)
library(purrr)
library(pbmcapply)
library(metafor)
library(pbapply)
library(concaveman)


# SIMULATION PARAMETERS
N_SIM <- 100
WINDOW_SIZE <- 2000
R <- 200

BASE_DIR <- "/Users/jacobchang/Lab/dcspomic_paper"


results <- df
intermediate_level <- results |>
  dplyr::group_by(patient, i_j) |>
  dplyr::summarise(patient_mean = mean(colocalization_stat),
                   patient_sigma2 = var(colocalization_stat))

population_colocalization <- intermediate_level |>
  dplyr::group_by(i_j) |>
  dplyr::summarise(mu = mean(patient_mean),
                   tau2 = var(patient_mean))

total_level2 <- results |>
  dplyr::group_by(i_j) |>
  dplyr::summarise(mu = mean(colocalization_stat),
                   tau2 = var(colocalization_stat))



################################################################################
###

get_rmse <- function(truth, estimates) {
  return(sqrt(sum((estimates-truth)^2, na.rm=TRUE) / length(estimates)))
}


n_inter <- 100
n_intra <- 100

# Iterate through each "patient"
# Starting at two because there is no inter-patient variability with a single patient
pair <- "B_A"
results_ab <- results |> dplyr::filter(i_j == pair) # is this going to break since it is not "(A)_(B)"?
head(results_ab)
vec <- 1:100
half <- length(vec) / 2
idxs <- as.vector(rbind(vec[1:half], vec[length(vec):(half + 1)]))
blah = results_ab |>
  group_by(patient) |>
  summarise(v = var(colocalization_stat)) |>
  ungroup() |>
  arrange(desc(v)) |>
  mutate(id = row_number()) |>
  select(patient, id)
head(blah)
results_ab <- inner_join(results_ab, blah)
head(results_ab)

population_colocalization_ab <- population_colocalization |> dplyr::filter(i_j == pair)
set.seed(123)
mu_hats <- c(); tau2_hats <- c();
naive_mu_hats <- c(); naive_tau2_hats <- c();
random_intra_sample <- sample(1:100, replace = FALSE)
df_list <- list()


# for(k in random_intra_sample){
counter = 1
for(k in random_intra_sample){

  print(counter)
  counter = counter + 1
for(i in 2:n_inter) {
  reduced_result <- results_ab |>
    # arrange(desc(colocalization_var)) |>
    # group_by(patient) |>
    # summarise(id = row_number()) |>
    # ungroup()
    # mutate(id = row_number())
    # dplyr::filter(patient %in% 1:i) |>
    dplyr::filter(id %in% idxs[1:i]) |>

    # dplyr::filter(intra_sample == sample(1:n_intra, size = 1))
    dplyr::filter(intra_sample == k)

  model <- metafor::rma.uni(yi = colocalization_stat, vi = colocalization_var, method = "SJ", data = reduced_result)
  # model <- metafor::rma.uni(yi = colocalization_stat, vi = colocalization_var, method = "DL", data = reduced_result)

  # forest(model)

  mu_hats[i] <- coef(summary(model))$estimate
  tau2_hats[i] <- summary(model)$tau2
  naive_mu_hats[i] <- mean(reduced_result$colocalization_stat)
  naive_tau2_hats[i] <- var(reduced_result$colocalization_stat)

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
  mutate(n = as.factor(n)) |>
  group_by(n) |>
  summarise(
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
    naive_tau2_hat = mean(naive_tau2_hat, na.rm = TRUE)
    )

dcspomic_df |>
  mutate(n = factor(n)) |>
  tidyplot(x = n, y = tau2_hat) |>
  add_boxplot()

dcspomic_df |>
  mutate(n = factor(n)) |>
  tidyplot(x = n, y = naive_tau2_hat) |>
  add_boxplot()
  # add_data_points_beeswarm()
  # add_data_points() |>
  # add

mu_hat_rmses <- c(); naive_mu_hat_rmses <- c();
tau2_hat_rmses <- c(); naive_tau2_hat_rmses <- c();
for(i in 2:n_inter) {
  # foo2 <- foo |> filter(n %in% 1:i)
  foo2 <- dcspomic_df |> filter(n %in% 1:i)

  mu_hat_rmses[i] <- get_rmse(truth=foo2$true_mu, estimates=foo2$mu_hat)
  naive_mu_hat_rmses[i] <- get_rmse(truth=foo2$true_mu, estimates=foo2$naive_mu_hat)
  tau2_hat_rmses[i] <- get_rmse(truth=foo2$true_tau2, estimates=foo2$tau2_hat)
  naive_tau2_hat_rmses[i] <- get_rmse(truth=foo2$true_tau2, estimates=foo2$naive_tau2_hat)

}
df <- dcspomic_df |>
# df <- foo |>
  mutate(mu_hat_rmse = mu_hat_rmses,
         naive_mu_hat_rmse = naive_mu_hat_rmses,
         tau2_hat_rmse = tau2_hat_rmses,
         naive_tau2_hat_rmse = naive_tau2_hat_rmses)



# Pivot for ggplot
df_mu_long <- df %>%
  select(n, mu_hat_rmse, naive_mu_hat_rmse) %>%
  rename("DC-SPOMIC" = "mu_hat_rmse",
         "naive" = "naive_mu_hat_rmse") |>
  pivot_longer(cols = c("DC-SPOMIC", naive), names_to = "Method", values_to = "RMSE")

df_tau2_long <- df %>%
  select(n, tau2_hat_rmse, naive_tau2_hat_rmse) %>%
  rename("DC-SPOMIC" = "tau2_hat_rmse",
         "naive" = "naive_tau2_hat_rmse") |>
  pivot_longer(cols = c("DC-SPOMIC", naive), names_to = "Method", values_to = "RMSE")


# Assuming df_mu_long is as before
library(tidyplots)

upper_limit <- max(c(df_mu_long$RMSE, df_tau2_long$RMSE), na.rm=TRUE)

df_mu_long |>
  mutate(n = as.numeric(n)) |>
  tidyplot(x = n, y = RMSE, color = Method) |>
  add_data_points(size = 0.25) |>
  add_line() |>
  adjust_font(fontsize=5) |>
  adjust_title("Comparison of estimates for \u03BC", fontsize = 7) |>
  adjust_x_axis_title("Number of samples") |>
  adjust_y_axis_title(expression("RMSE of \u03BC")) |>
  remove_legend_title() |>
  adjust_legend_position(position = "bottom") |>
  adjust_size(width = 3, height = 2, unit = "in") |>
  adjust_y_axis(limits = c(0, sqrt(upper_limit))) |>
  adjust_x_axis(cut_short_scale = TRUE)

tidyplots::save_plot(
  plot = ggplot2::last_plot(),
  filename = paste0("output/simulation/intersample_variance/mu_rmse_plot.png")
  )
tidyplots::save_plot(
  plot = ggplot2::last_plot(),
  filename = paste0("output/simulation/intersample_variance/mu_rmse_plot.pdf")
)

df_tau2_long |>
  mutate(n = as.numeric(n)) |>
  tidyplot(x = n, y = RMSE, color = Method) |>
  add_data_points(size = 0.25) |>
  add_line() |>
  adjust_font(fontsize=5) |>
  adjust_title("Comparison of estimates for \u03C4\u00B2", fontsize = 7) |>
  adjust_x_axis_title("Number of samples") |>
  adjust_y_axis_title(expression("RMSE of \u03C4\u00B2")) |>
  remove_legend_title() |>
  adjust_legend_position(position = "bottom") |>
  adjust_size(width = 3, height = 2, unit = "in") |>
  adjust_y_axis(limits = c(0, upper_limit))

tidyplots::save_plot(
  plot = ggplot2::last_plot(),
  filename = paste0("output/simulation/intersample_variance/tau2_rmse_plot.png")
  )
tidyplots::save_plot(
  plot = ggplot2::last_plot(),
  filename = paste0("output/simulation/intersample_variance/tau2_rmse_plot.pdf")
)

# ==============================================================================


