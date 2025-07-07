simulation_colors <-
  new_color_scheme(c("A" = "#56B4E9",
                     "B" = "#E69F00",
                     "(C)" = "#D3D3D3",
                     "(Rare 1)" = "#009E73",
                     "(Rare 2)" = "#D55E00"),
                   name = "simulation_color_scheme")

for(i in 1:10) {
  spomic <- readRDS(paste0("output/simulation/intersample_variance/intersample_spomic_", i, ".rds"))
  plot_spomic(spomic) |>
    adjust_colors(new_colors = simulation_colors)



  save_plot(plot = ggplot2::last_plot(),
            paste0("output/simulation/intersample_variance/intersample_spomic_", i, ".png"),
            bg = "transparent",unit = "in", height = 3, width = 3)
  save_plot(plot = ggplot2::last_plot(),
            paste0("output/simulation/intersample_variance/intersample_spomic_", i, ".pdf"),
            bg = "transparent",unit = "in", height = 3, width = 3)
}

