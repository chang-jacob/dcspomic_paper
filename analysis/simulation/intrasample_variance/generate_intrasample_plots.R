# devtools::load_all("/home/groups/plevriti/jachang4/spomic")
devtools::load_all("/Users/jacobchang/Lab/spomic")

library(tidyplots)

oak_path <- "/oak/stanford/groups/plevriti"
simulation_colors <- new_color_scheme(x = c("(A)" = "#56B4E9",
                                            "(B)" = "#E69F00",
                                            "(C)" = "#D3D3D3",
                                            "(Rare 1)" = "#009E73",
                                            "(Rare 2)" = "#D55E00"),
                                      name = "simulation_color_scheme")

homogeneous_spomics <- readRDS(file.path(oak_path,
                  "jachang4",
                  "dcspomic_output",
                  "intrasample_simulation",
                  "homogeneous_pattern_spomics.rds"))

p <- plot_spomic(homogeneous_spomics[[1]]) %>%
  adjust_colors(new_colors = simulation_colors) %>%
  remove_legend() %>%
  remove_caption() %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_y_axis_labels() %>%
  remove_y_axis_title() %>%
  remove_y_axis_ticks()
  # save_plot("output/simulation/intrasample_variance/homogeneous_process/example1.svg")

ggplot2::ggsave(filename="output/simulation/intrasample_variance/homogeneous_process/example1.svg", plot=p)
