## ---------------------------
##
## Script name: network_analysis.R
##
## Author: Jake Chang
##
## Date Modified: 2025-08-19
##
## ---------------------------
##
#' Description:
#' This script performs network analysis of colocalization scores
##
## ---------------------------

## load libraries
library(spomic)
library(dcspomic)

library(metafor)
library(tidyverse)
library(tidyplots)
library(pbapply)
library(beepr)
library(janitor)
library(igraph)
library(tidyr)
library(dplyr)

library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)

scratch <- config::get("scratch")
oak <- config::get("oak")
date <- format(Sys.Date(), "%Y%m%d")

source("analysis/plotting_utils.R")

# Group 1: Crohn's like reaction (CLR)
# Group 2: diffuse inflammatory infiltration (DII)

## ---------------------------

## run DC-SPOMIC

group1_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group1_spomics_Lcross_global_envelope.rds"))
group2_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group2_spomics_Lcross_global_envelope.rds"))



## Colocalization network analysis

dcspomic <- create_dcspomic(
  group1_name = "CLR",
  group1_spomics = group1_spomics,
  group2_name = "DII",
  group2_spomics = group2_spomics)
dcspomic <- set_dcspomic_hyperparameters(
  dcspomic,
  r = 75, # 75um
  colocalization_type = "Lcross",
  tau_estimator = "SJ")

dcspomic <- run_dcspomic(dcspomic)

# Remove the cell pairs containing undefined cells and recompute adjusted p-vals
results_without_unknowns <- dcspomic@results$differential_testing[
  !grepl("\\(undefined\\)", dcspomic@results$differential_testing$i_j),
]
results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
results_without_unknowns$holm <- p.adjust(results_without_unknowns$pval, method="holm")
dcspomic@results$differential_testing <- results_without_unknowns

results_sep <- results_without_unknowns %>%
  separate(i_j, into = c("row", "col"), sep = "\\)_\\(", remove = FALSE) %>%
  mutate(
    row = gsub("\\(", "", row),
    col = gsub("\\)", "", col)
  )

## ---------------------------
# Build network 

all_labels <- sort(unique(c(results_sep$row, results_sep$col)))
adj_matrix <- matrix(NA, nrow = length(all_labels), ncol = length(all_labels),
                     dimnames = list(all_labels, all_labels))

# Create matrix with Z-scores from differential testing
for(i in seq_len(nrow(results_sep))) {
  r <- results_sep$row[i]
  c <- results_sep$col[i]
  adj_matrix[r, c] <- results_sep$z_score[i]
}

# Only take the Zscores related to higher colocalization in CLR samples
adj_matrix_df <- as.data.frame(adj_matrix)
adj_matrix_neg <- adj_matrix
adj_matrix_neg[adj_matrix_neg > 0 | is.na(adj_matrix_neg)] <- 0
adj_matrix_neg <- abs(adj_matrix_neg)

adj_matrix_df <- as.data.frame(adj_matrix_neg)

# Make the matrix symmetric to build an indirected graph
sym_adj_matrix <- (adj_matrix_neg + t(adj_matrix_neg)) / 2

colocalization_matrix <- results_without_unknowns

# Convert to graph from adjacency matrix
g <- graph_from_adjacency_matrix(
  sym_adj_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Optional: customize edge width by weight
E(g)$width <- E(g)$weight * 2  # scale for visibility

# # Optional: set node size by degree or other metric
V(g)$size <- degree(g, mode = "in") + 5

# # Optional: add color
E(g)$color <- ifelse(E(g)$weight > 1.96, "red", "gray")

# Run Louvain clustering
cl <- cluster_louvain(g, resolution = 1.25)

# Assign cluster to each node
V(g)$group <- membership(cl)

# Convert igraph to tidygraph
tg <- as_tbl_graph(g)

# Create ggraph plot with label repulsion
set.seed(123)
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "#88888833") + 
  scale_edge_width(range = c(0.2, 2.5)) + 
  
  # Nodes and labels
  geom_node_point(aes(color = as.factor(group)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  labs(color = "Leiden cluster (res=1.25)") + 
  theme_void() +
  theme(legend.position="bottom") + 
  guides(edge_width = "none") 

ggsave("output/crc_tma/lcross_75um/clr_network.png", plot = last_plot(), width = 6, height = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/clr_network.pdf", plot = last_plot(), width = 6, height = 4, units = "in")
ggsave("output/crc_tma/lcross_75um/clr_network.svg", plot = last_plot())

## ---------------------------

# test immune cluster
# Define your 4 main cell types
cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")

# 5-color scheme: replace/add your preferred hex codes
five_colors <- tidyplots::new_color_scheme(
  x = c("#d7d2cb", "#0371b2", "#d36027", "#17b872", "#bb27d3")
)

# Map cell types in the data frame
spomic <- dcspomic@group1_spomics[[8]]
plot_quad <- function(spomic){
  spomic@df |>
    dplyr::mutate(
      cell_type2 = dplyr::case_when(
        cell_type == "(b_cells)" ~ "(b_cells)",
        cell_type == "(tregs)" ~ "(tregs)",
        cell_type == "(cd8_t_cells)" ~ "(cd8_t_cells)",
        cell_type == "(cd4_t_cells_cd45ro)" ~ "(cd4_t_cells_cd45ro)",
        TRUE ~ "(other)"
      )
    ) |>
    dplyr::mutate(
      cell_type2 = factor(cell_type2, levels = c("(other)", "(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)"))
    ) |>
    dplyr::arrange(cell_type2) |>
    tidyplots::tidyplot(x = x, y = y, color = cell_type2) |>
    tidyplots::add_data_points(size = 0.1) |>
    tidyplots::adjust_colors(new_colors = five_colors) |>
    tidyplots::remove_padding() |>
    tidyplots::adjust_title(spomic@details$sample) |>
    tidyplots::adjust_legend_title("Cell type") |>

    spomic_style()
}

plot_quad(dcspomic@group1_spomics[[8]]) |>
  tidyplots::save_plot(filename = "output/crc_tma/tls_cluster_group1.pdf")

plot_quad(dcspomic@group2_spomics[[7]]) |>
  tidyplots::save_plot(filename = "output/crc_tma/tls_cluster_group2.pdf")

## ---------------------------

# Test for differences in average distance in immune cluster

library(dplyr)
library(RANN)

cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")
# df_sub <- dcspomic@group1_spomics[[1]]@df %>% filter(cell_type %in% cell_types)

# Initialize result dataframe
result <- data.frame(
  cell_id = df_sub$cell_id,
  cell_type = df_sub$cell_type,
  nn_b_cells = NA,
  d_b_cells = NA,
  nn_tregs = NA,
  d_tregs = NA,
  nn_cd8 = NA,
  d_cd8 = NA,
  nn_cd4 = NA,
  d_cd4 = NA,
  dist_sum = NA
)

for (i in 1:nrow(df_sub)) {
  this_cell <- df_sub[i, ]
  this_xy <- as.numeric(this_cell[, c("x", "y")])
  dists <- c()
  ids <- c()
  
  for (target_type in cell_types) {
    if (this_cell$cell_type == target_type) {
      # Find nearest neighbor of *same* type, excluding self
      others <- df_sub %>% filter(cell_type == target_type, cell_id != this_cell$cell_id)
    } else {
      others <- df_sub %>% filter(cell_type == target_type)
    }
    # KD-tree is overkill for 1 query, use base distance
    if (nrow(others) == 0) {
      ids <- c(ids, NA)
      dists <- c(dists, NA)
    } else {
      xy_mat <- as.matrix(others[, c("x", "y")])
      diff <- sweep(xy_mat, 2, this_xy, "-")
      euclid <- sqrt(rowSums(diff^2))
      min_idx <- which.min(euclid)
      ids <- c(ids, others$cell_id[min_idx])
      dists <- c(dists, euclid[min_idx])
    }
  }
  # Assign results
  result$nn_b_cells[i] <- ids[1]
  result$d_b_cells[i] <- dists[1]
  result$nn_tregs[i] <- ids[2]
  result$d_tregs[i] <- dists[2]
  result$nn_cd8[i] <- ids[3]
  result$d_cd8[i] <- dists[3]
  result$nn_cd4[i] <- ids[4]
  result$d_cd4[i] <- dists[4]
  result$dist_sum[i] <- sum(dists, na.rm = TRUE)
}
distance_cols <- grep("d_", names(result), value = TRUE)
result[distance_cols] <- lapply(result[distance_cols], function(x) pmin(x, 75))
result$dist_sum <- result$d_b_cells +  result$d_tregs + result$d_cd8 + result$d_cd4

result %>% group_by(cell_type) %>% summarise(mean(dist_sum, na.rm = TRUE))
mean(result$dist_sum, na.rm=TRUE)
hist(result$dist_sum, na.rm=TRUE)

grab_immune_cluster_dist <- function(spomic) {
  cell_types <- c("(b_cells)", "(tregs)", "(cd8_t_cells)", "(cd4_t_cells_cd45ro)")
  df_sub <- spomic@df %>% filter(cell_type %in% cell_types)
  
  # Initialize result dataframe
  result <- data.frame(
    cell_id = df_sub$cell_id,
    cell_type = df_sub$cell_type,
    nn_b_cells = NA,
    d_b_cells = NA,
    nn_tregs = NA,
    d_tregs = NA,
    nn_cd8 = NA,
    d_cd8 = NA,
    nn_cd4 = NA,
    d_cd4 = NA,
    dist_sum = NA
  )
  
  for (i in 1:nrow(df_sub)) {
    this_cell <- df_sub[i, ]
    this_xy <- as.numeric(this_cell[, c("x", "y")])
    dists <- c()
    ids <- c()
    
    for (target_type in cell_types) {
      if (this_cell$cell_type == target_type) {
        # Find nearest neighbor of *same* type, excluding self
        others <- df_sub %>% filter(cell_type == target_type, cell_id != this_cell$cell_id)
      } else {
        others <- df_sub %>% filter(cell_type == target_type)
      }
      # KD-tree is overkill for 1 query, use base distance
      if (nrow(others) == 0) {
        ids <- c(ids, NA)
        dists <- c(dists, NA)
      } else {
        xy_mat <- as.matrix(others[, c("x", "y")])
        diff <- sweep(xy_mat, 2, this_xy, "-")
        euclid <- sqrt(rowSums(diff^2))
        min_idx <- which.min(euclid)
        ids <- c(ids, others$cell_id[min_idx])
        dists <- c(dists, euclid[min_idx])
      }
    }
    # Assign results
    result$nn_b_cells[i] <- ids[1]
    result$d_b_cells[i] <- dists[1]
    result$nn_tregs[i] <- ids[2]
    result$d_tregs[i] <- dists[2]
    result$nn_cd8[i] <- ids[3]
    result$d_cd8[i] <- dists[3]
    result$nn_cd4[i] <- ids[4]
    result$d_cd4[i] <- dists[4]
    result$dist_sum[i] <- sum(dists, na.rm = TRUE)
    result$dist_avg[i] <- mean(dists, na.rm = TRUE)
    
  }
  
  distance_cols <- grep("d_", names(result), value = TRUE)
  result[distance_cols] <- lapply(result[distance_cols], function(x) pmin(x, 75))
  result$dist_sum <- result$d_b_cells +  result$d_tregs + result$d_cd8 + result$d_cd4
  result <- result %>% rowwise() %>% mutate(dist_avg = mean(d_b_cells, d_tregs, d_cd8, d_cd4, na.rm=TRUE))

  result$dist_avg
  
}

g1_immune_dist <- c()
for(i in 1:length(dcspomic@group1_spomics)) {
  print(i)
  g1_immune_dist[i] <- grab_immune_cluster_dist(dcspomic@group1_spomics[[i]])
}

g2_immune_dist <- c()
for(i in 1:length(dcspomic@group2_spomics)) {
  print(i)
  g2_immune_dist[i] <- grab_immune_cluster_dist(dcspomic@group2_spomics[[i]])
}

t.test(g1_immune_dist, g2_immune_dist, alternative = "less")

df <- data.frame(
  value = c(g1_immune_dist, g2_immune_dist),
  group = rep(c("CLR", "DII"),
              c(length(g1_immune_dist), length(g2_immune_dist)))
)


mean(g1_immune_dist) # sample index 8
mean(g2_immune_dist) #sample index 7

t.test(g1_immune_dist, g2_immune_dist, alternative = "less")

library(ggplot2)

# Ensure your groups are ordered so the test is "0" > "1" (change as needed)
df$group <- factor(df$group, levels = c("CLR","DII"))

# One-sided Welch t-test (H1: mean(group "0") > mean(group "1"))
tt <- t.test(value ~ group, data = df, alternative = "less")
p  <- tt$p.value
lab <- if (p < 1e-4) "p < 1e-4 (one-sided t-test)" else sprintf("p = %.4f (one-sided t-test)", p)

# Y positions for the bracket/label
yr <- range(df$value, na.rm = TRUE)
yspan <- diff(yr)
y_bracket <- yr[2] + 0.06 * yspan   # horizontal line height
y_bracket <- 125
tick     <- 7.5          # little vertical tick height

p = ggplot(df, aes(x = group, y = value, color = group)) +
  geom_violin(alpha = 0.15, trim = FALSE) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
  labs(y = "Average distance between cells in immune subnetwork", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  # bracket (1..2 on the x-axis since there are two discrete groups)
  geom_segment(aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE) +
  geom_segment(aes(x = 1, xend = 1, y = y_bracket - tick, yend = y_bracket), inherit.aes = FALSE) +
  geom_segment(aes(x = 2, xend = 2, y = y_bracket - tick, yend = y_bracket), inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = y_bracket + 0.03 * yspan, label = lab) + 
  ggpubr::theme_pubr() + 
  theme(legend.position="none")
p
ggsave(file.path("output/crc_tma/network_comparison.pdf"), height = 5, width = 4, units = "in")
ggsave(file.path("output/crc_tma/network_comparison.png"), height = 5, width = 4, units = "in")
