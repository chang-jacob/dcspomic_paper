## ---------------------------
##
## Script name: crc_dcspomic_multiscale.R
##
## Author: Jake Chang
##
## Date Modified: 2025-07-09
##
## ---------------------------
##
#' Description:
#' This script assesses the stability of DC-SPOMIC results over various length 
#' scales. 
##
## ---------------------------
devtools::load_all("/home/groups/plevriti/jachang4/spomic")
devtools::load_all("/home/groups/plevriti/jachang4/dcspomic")
library(metafor)
library(tidyverse)
library(tidyplots)
library(pbapply)
library(beepr)
library(janitor)

group1_spomics <- readRDS("output/crc_tma/group1_spomics_Lcross_global_envelope.rds")
group2_spomics <- readRDS("output/crc_tma/group2_spomics_Lcross_global_envelope.rds")

volcano_list <- list()
for(r in c(25, 50, 75, 100, 125, 150, 200)) {
  print(r)
  dcspomic <- create_dcspomic(
    group1_name = "CLR",
    group1_spomics = group1_spomics,
    group2_name = "DII",
    group2_spomics = group2_spomics)
  
  dcspomic <- set_dcspomic_hyperparameters(
    dcspomic,
    r = r,
    colocalization_type = "Lcross",
    tau_estimator = "SJ")
  
  dcspomic <- run_dcspomic(dcspomic)
  results_without_unknowns <- dcspomic@results$differential_testing[
    !grepl("\\(undefined\\)", dcspomic@results$differential_testing$i_j),
  ]
  results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
  results_without_unknowns$holm <- p.adjust(results_without_unknowns$pval, method="holm")
  dcspomic@results$differential_testing <- results_without_unknowns
  
  alpha <- 0.05
  p <- plot_volcano(dcspomic, alpha = alpha)
  volcano_list[[as.character(r)]] <- p
  
  p %>% 
    save_plot(filename = paste0("output/crc_tma/lcross_multiscale/volcano_lcross_", r, "um.png")) %>% 
    save_plot(filename = paste0("output/crc_tma/lcross_multiscale/volcano_lcross_", r, "um.pdf")) 
  # %>% 
    # save_plot(filename = paste0("output/crc_tma/lcross_multiscale/volcano_lcross_", r, "um.svg"))
  
  
  
}


