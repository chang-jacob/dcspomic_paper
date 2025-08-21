# Explore the potential (likely) relationship between variance estimated by 
# spatial bootstrap and the population size of cell types 

# Hypothesis: 
# I think that most of the time, these will be inversely related. I think for 
# instances where they do not align, it could be because of a dense and  
# isolated group of a relatively rare cell type

library(spomic)
library(dcspomic)

oak <- config::get("oak")

group1_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group1_spomics_Lcross_global_envelope.rds"))
group2_spomics <- readRDS(file.path(oak, "dcspomic_output/crc_tma/group2_spomics_Lcross_global_envelope.rds"))

get_colocalization_and_population <- function(spomic_list, r){
  lapply(spomic_list, function(x) {
    # x = spomic_list[[1]]
    count_table <- table(x@pp$marks)
    count_df <- data.frame(
      cell_type = names(count_table), 
      cell_count = as.vector(count_table))
    count_df$prop <- count_df$cell_count / sum(count_table)
    
    colocalization_df <- x %>%
      set_spomic_hyperparameters(r = r) %>%
      get_spatial_summary()
    
    spl <- strsplit(colocalization_df$i_j, ")_(", fixed = TRUE)
    colocalization_df$i <- sub("^\\(", "", sapply(spl, `[`, 1))
    colocalization_df$j <- sub("\\)$", "", sapply(spl, `[`, 2))
    
    colocalization_df$i <- paste0(sapply(spl, `[`, 1), ")")
    colocalization_df$j <- paste0("(", sapply(spl, `[`, 2))
    
    full_colocalization_df <- colocalization_df %>% 
      left_join(count_df, by=c("i"="cell_type")) %>% 
      rename("N_i" = "cell_count", "prop_i" = "prop") %>% 
      left_join(count_df, by=c("j"="cell_type")) %>% 
      rename("N_j" = "cell_count", "prop_j" = "prop") 
  })
}

g1_df <- bind_rows(
  get_colocalization_and_population(group1_spomics, r=75)
)

g1_df %>% 
  ggplot(aes(x=N_i, y = colocalization_var, color=sample)) + 
  geom_point() + 
  theme_minimal()

g1_df %>% 
  ggplot(aes(x=N_j, y = colocalization_var, color=sample)) + 
  geom_point() + 
  theme_minimal()

g1_df %>% 
  ggplot(aes(x=prop_i, y = colocalization_var, color=sample)) + 
  geom_point() + 
  theme_minimal()

g1_df %>% 
  ggplot(aes(x=prop_j, y = colocalization_var, color=sample)) + 
  geom_point() + 
  theme_minimal()

g1_df %>% 
  ggplot(aes(x=prop_i, y = prop_j, color=colocalization_var)) + 
  geom_point() + 
  theme_minimal()

g1_df %>% 
  ggplot(aes(x=sqrt(N_i), y = sqrt(N_j), color=colocalization_var)) + 
  geom_point() + 
  theme_minimal()

all_g2_spomics_stats <- bind_rows(
  lapply(group2_spomics, function(x) {
    x %>%
      set_spomic_hyperparameters(r = 75) %>%
      get_spatial_summary()
  })
)