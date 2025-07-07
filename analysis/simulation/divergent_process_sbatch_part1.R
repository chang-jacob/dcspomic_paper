devtools::load_all("/home/groups/plevriti/jachang4/spomic")
library(spatstat)
library(dplyr)
library(progress)
library(tidyplots)
# library(ggplot2)
# library(ggpubr)
# library(ggrepel)
# library(cowplot)
# library(patchwork)

# Define fixed intensities for example point process models
square_size <- 2000
window <- owin(square(square_size))

cell_types <- c("(A)", "(B)", "(C)", "(Rare 1)", "(Rare 2)")


# Define lambda functions
LAMBDA_A <- function(x, y, rate) {
  dexp(x, rate = rate)

}
LAMBDA_B <- function(x, y, rate) {
  max_value <- 2000
  dexp(max_value - x, rate = rate)
}

# Intensity constants
LAMBDA_C <- 0.1/1000
LAMBDA_RARE1 <- 0.01/1000
LAMBDA_RARE2 <- 0.005/1000

LAMBDA_B(x=0:2000, y=0:2000, rate=0.004) |> plot()

divergent_spomic <- list()
divergent_intensity_plots <- list()
divergent_spomic_colocalization_results <- list()
set.seed(345)
# rate <- runif(100, min=1e-10, max = 0.005)
rate <- abs(rnorm(1, 0.001, 0.001))
hist(rate)

rate = 0.0025
LAMBDA_A(x=0:2000, y=0:2000, rate) |> plot()
LAMBDA_B(x=0:2000, y=0:2000, rate) |> plot()

for(i in 1:100) {
  print(i)


  # Generate each component point process
  pp1 <- rpoispp(lambda = function(x,y) LAMBDA_A(x, y, rate), win = window)
  pp1$marks <- factor(rep("(A)", pp1$n), levels = cell_types)
  # plot(pp1)

  pp2 <- rpoispp(lambda = function(x,y) LAMBDA_B(x, y, rate), win = window)
  pp2$marks <- factor(rep("(B)", pp2$n), levels = cell_types)
  # plot(pp2)

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
  # plot_spomic(spomic) |> print()
  # plot_cell_proportions(spomic)
  spomic <- set_spomic_hyperparameters(spomic=spomic, r = 200, fixed_distance = TRUE, colocalization_type="Lcross")
  spomic <- get_spatial_stats(spomic)
  divergent_spomic_colocalization_results[[i]] <- get_spatial_summary(spomic)

  if(i == 1) divergent_spomic[[i]] <- spomic

  divergent_intensity_plots[[i]] <- plot_spomic(spomic)
}


