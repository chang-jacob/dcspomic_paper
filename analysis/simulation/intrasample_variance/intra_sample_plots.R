intra_sample_n <- 4
WINDOW_SIZE <- 2000
R <- 200

BASE_DIR <- "/Users/jacobchang/Lab/dcspomic_paper"
set.seed(123)

parent_lambda = 0.00001 + runif(1, 0, 0.0001)
offspring_n = 3 + rnbinom(1, size=10, mu=2)
offspring_dispersion = abs(200 + rnorm(1, mean=0, sd = 100))

  intra_sample_list <- list()
  for(n in 1:intra_sample_n) {
    # parent_lambda <- parent_lambda0 + runif(1, 0, 1/(window_size * window_size))
    # offspring_n <- rpois(n = 1, offspring_n0) + 1
    # offspring_dispersion <- abs(offspring_dispersion0 + rnorm(1, 0, 15))

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

    parent_ppp <- ppp(parents$x, parents$y, window = square(window_size), marks = rep("(A)", parents$n))
    offspring_ppp <- ppp(offspring_x, offspring_y, window = square(window_size), marks = rep("(B)", length(offspring_x)))
    marked_process <- superimpose(parent_ppp, offspring_ppp)
    marks(marked_process) <- as.factor(marks(marked_process))

    spomic <- create_spomic(ppp_to_df(marked_process, sample_name = ""))
    spomic <- set_spomic_hyperparameters(spomic = spomic,
                                         colocalization_type = "Lcross",
                                         fixed_distance = FALSE,
                                         r = R)
    intra_sample_list[[n]] <- spomic
    simulation_colors <-
      new_color_scheme(c("(A)" = "#56B4E9",
                         "(B)" = "#E69F00",
                         "(C)" = "#D3D3D3",
                         "(Rare 1)" = "#009E73",
                         "(Rare 2)" = "#D55E00"),
                       name = "simulation_color_scheme")

    plot_spomic(intra_sample_list[[n]]) |>
      adjust_colors(new_colors = simulation_colors) |>
      tidyplots::save_plot(paste0("output/simulation/intrasample_example", n, ".pdf"),
                           unit = "in", height = 3, width = 3)

  }




