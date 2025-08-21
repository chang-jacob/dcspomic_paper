## Does DC-SPOMIC outperform t-test with smaller sample sizes?
set.seed(123)
shuffled_group1_spomics <- group1_spomics[sample(length(group1_spomics))]
shuffled_group2_spomics <- group2_spomics[sample(length(group2_spomics))]

results_list <- list()
t_subset_results <- list()

dcspomic_pval_res <- list()
t_pval_res <- list()

for(i in 5:17){
  dcspomic_subset <- create_dcspomic(
    group1_name = "CLR",
    group1_spomics = shuffled_group1_spomics[1:i],
    group2_name = "DII",
    group2_spomics = shuffled_group2_spomics[1:i])
  dcspomic_subset <- set_dcspomic_hyperparameters(
    dcspomic_subset,
    r = 75, # 75um
    colocalization_type = "Lcross",
    tau_estimator = "SJ")
  
  dcspomic_subset <- run_dcspomic(dcspomic_subset)
  
  results_without_unknowns <- dcspomic_subset@results$differential_testing[
    !grepl("\\(undefined\\)", dcspomic_subset@results$differential_testing$i_j),
  ]
  results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
  results_without_unknowns$holm <- p.adjust(results_without_unknowns$pval, method="holm")
  dcspomic_subset@results$differential_testing <- results_without_unknowns
  
  results_list[[i]] <- results_without_unknowns %>% filter(FDR < 0.05) %>% pull(i_j)
  
  dcspomic_pval_res[[i]] <- results_without_unknowns %>% mutate(n_samples= i, type = "DCSPOMIC") 
  
  group1_summaries <- lapply(dcspomic_subset@group1_spomics, get_spatial_summary)
  group1_stats <- bind_rows(group1_summaries) %>% filter(!grepl("undefined", i_j))
  
  group2_summaries <- lapply(dcspomic_subset@group2_spomics, get_spatial_summary)
  group2_stats <- bind_rows(group2_summaries) %>% filter(!grepl("undefined", i_j))
  
  cell_pairs <- unique(c(group1_stats$i_j, group2_stats$i_j))
  
  t_results <- bind_rows(lapply(cell_pairs, function(pair) {
    x <- group1_stats %>% filter(i_j == pair) %>% pull(colocalization_stat)
    y <- group2_stats %>% filter(i_j == pair) %>% pull(colocalization_stat)
    t <- t.test(x, y)
    tibble(
      i_j = pair,
      p_value = t$p.value,
      statistic = unname(t$statistic), 
      log2fc = log2(t$estimate[2]+1) - log2(t$estimate[1]+1)
    )
  }))
  t_results$FDR <- p.adjust(t_results$p_value, method = "fdr")
  t_results$holm <- p.adjust(t_results$p_value, method = "holm")
  t_subset_results[[i]] <- t_results %>% filter(FDR < 0.05) %>% pull(i_j)
  t_pval_res[[i]] <- t_results %>% mutate(n_samples = i, type = "t") 
  
  print(i)
}

# A little extra work to get the 18th sample from group2 involved.
dcspomic_subset <- create_dcspomic(
  group1_name = "CLR",
  group1_spomics = shuffled_group1_spomics[1:17],
  group2_name = "DII",
  group2_spomics = shuffled_group2_spomics[1:18])
dcspomic_subset <- set_dcspomic_hyperparameters(
  dcspomic_subset,
  r = 75, # 75um
  colocalization_type = "Lcross",
  tau_estimator = "SJ")

dcspomic_subset <- run_dcspomic(dcspomic_subset)

results_without_unknowns <- dcspomic_subset@results$differential_testing[
  !grepl("\\(undefined\\)", dcspomic_subset@results$differential_testing$i_j),
]
results_without_unknowns$FDR <- p.adjust(results_without_unknowns$pval, method="fdr")
results_without_unknowns$holm <- p.adjust(results_without_unknowns$pval, method="holm")
dcspomic_subset@results$differential_testing <- results_without_unknowns

results_list[[18]] <- results_without_unknowns %>% filter(FDR < 0.05) %>% pull(i_j)

group1_summaries <- lapply(dcspomic_subset@group1_spomics, get_spatial_summary)
group1_stats <- bind_rows(group1_summaries) %>% filter(!grepl("undefined", i_j))

group2_summaries <- lapply(dcspomic_subset@group2_spomics, get_spatial_summary)
group2_stats <- bind_rows(group2_summaries) %>% filter(!grepl("undefined", i_j))

cell_pairs <- unique(c(group1_stats$i_j, group2_stats$i_j))
t_results <- bind_rows(lapply(cell_pairs, function(pair) {
  x <- group1_stats %>% filter(i_j == pair) %>% pull(colocalization_stat)
  y <- group2_stats %>% filter(i_j == pair) %>% pull(colocalization_stat)
  t <- t.test(x, y)
  tibble(
    i_j = pair,
    p_value = t$p.value,
    statistic = unname(t$statistic), 
    log2fc = log2(t$estimate[2]+1) - log2(t$estimate[1]+1)
  )
}))
t_results$FDR <- p.adjust(t_results$p_value, method = "fdr")
t_results$holm <- p.adjust(t_results$p_value, method = "holm")
t_subset_results[[18]] <- t_results %>% filter(FDR < 0.05) %>% pull(i_j)

# foo1 <- bind_rows(dcspomic_pval_res)
# foo2 <- bind_rows(t_pval_res)