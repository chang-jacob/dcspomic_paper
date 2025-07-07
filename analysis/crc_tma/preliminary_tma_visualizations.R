# Eventually change the directories and imports...

devtools::load_all("/Users/jacobchang/Lab/spomic")
source("analysis/plotting_utils.R")

group1_spomics <- readRDS("output/crc_tma/processed_spomics/lcross/group1_spomics.rds")
group2_spomics <- readRDS("output/crc_tma/processed_spomics/lcross/group2_spomics.rds")


for(i in 1:length(group1_spomics)){
  spomic_scatter <- plot_spomic(group1_spomics[[i]])|>
    tidyplots::adjust_colors(new_colors = crc_tma_palette)
  spomic_proportions <- plot_cell_proportions(group1_spomics[[i]])|>
    tidyplots::adjust_colors(new_colors = crc_tma_palette)

  tidyplots::save_plot(
    plot = spomic_scatter,
    filename = paste0("output/crc_tma/preliminary_tma_visualizations/clr/", group1_spomics[[i]]@details$sample, "_scatter_plot.png")
    )
  tidyplots::save_plot(
    plot = spomic_proportions,
    filename = paste0("output/crc_tma/preliminary_tma_visualizations/clr/", group1_spomics[[i]]@details$sample, "_proportion_plot.png")
  )
}

for(i in 1:length(group2_spomics)){
  spomic_scatter <- plot_spomic(group2_spomics[[i]]) |>
    tidyplots::adjust_colors(new_colors = crc_tma_palette)

  spomic_proportions <- plot_cell_proportions(group2_spomics[[i]])|>
    tidyplots::adjust_colors(new_colors = crc_tma_palette)

  tidyplots::save_plot(
    plot = spomic_scatter,
    filename = paste0("output/crc_tma/preliminary_tma_visualizations/dii/", group2_spomics[[i]]@details$sample, "_scatter_plot.png")
  )
  tidyplots::save_plot(
    plot = spomic_proportions,
    filename = paste0("output/crc_tma/preliminary_tma_visualizations/dii/", group2_spomics[[i]]@details$sample, "_proportion_plot.png")
  )
}
