library(spomic)
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

intra_sample_simulation <- function(intra_sample_n,
                                    parent_lambda0,
                                    offspring_n0,
                                    offspring_dispersion0,
                                    patient = "",
                                    window_size = WINDOW_SIZE,
                                    plot = FALSE) {
  intra_sample_results_list <- list()
  for(n in 1:intra_sample_n) {
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
                     sd = offspring_dispersion)
      new_y <- rnorm(sum(num_offspring),
                     mean = rep(parents$y, num_offspring),
                     sd = offspring_dispersion)
      
      valid_idx <- (new_x >= 0) & (new_x <= window_size) & (new_y >= 0) & (new_y <= window_size)
      offspring_x <- new_x[valid_idx]
      offspring_y <- new_y[valid_idx]
    }
    
    parent_ppp <- ppp(parents$x, parents$y, window = square(window_size), marks = rep("A", parents$n))
    offspring_ppp <- ppp(offspring_x, offspring_y, window = square(window_size), marks = rep("B", length(offspring_x)))
    marked_process <- superimpose(parent_ppp, offspring_ppp)
    marks(marked_process) <- as.factor(marks(marked_process))
    
    spomic <- create_spomic(ppp_to_df(marked_process, sample_name = ""))
    spomic <- set_spomic_hyperparameters(spomic = spomic,
                                         colocalization_type = "Lcross",
                                         fixed_distance = FALSE,
                                         r = R)
    spomic <- get_spatial_stats(spomic)
    
    intra_sample_results_list[[n]] <- get_spatial_summary(spomic) |>
      mutate(patient = patient,
             intra_sample = n)
  }
  return(bind_rows(intra_sample_results_list))
}

inter_sample_simulation <- function(intra_sample_n) {
  parent_lambda0 = 0.00001 + runif(100, 0, 0.0001)
  offspring_n0 = 3 + rnbinom(100, size=10, mu=2)
  offspring_dispersion0 = abs(200 + rnorm(100, mean=0, sd = 100))
  
  results_list <- pblapply(1:100, function(i){
    print(i)
    intra_sample_simulation(
      intra_sample_n = intra_sample_n,
      parent_lambda0 = parent_lambda0[i],
      offspring_n0 = offspring_n0[i],
      offspring_dispersion0 = offspring_dispersion0[i],
      patient = i
    )
  })
  return(dplyr::bind_rows(results_list))
}

### WARNING: This next step takes a VERY long time!
# Commenting out to avoid running on accident.

# set.seed(123)
# cluster_process_simulation <- inter_sample_simulation(intra_sample_n=100)
# saveRDS(cluster_process_simulation, "cluster_process_df.rds")