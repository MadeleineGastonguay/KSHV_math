# This script runs the full pipeline for data from fixed images of cells 
# with full KSHV (Figure 5)

### Setup ###################################################################### 

# load libraries 
library(tidyverse)
library(here)
library(patchwork)

theme_set(theme_bw())

# Check wd
here()
setwd(here())

# Read in helper functions 
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

### Script inputs ##############################################################

# Folder for results
results_folder <- here("results", "fixed_KSHV")

# Daughter cell intensity data
daughter_cell_file <- "fixed_KSHV_dividing_cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "fixed_KSHV_non_dividing_cells.xlsx"

### Run analysis ###############################################################

# Load in the data
fixed_KSHV_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- fixed_KSHV_data$daughter_cell_data
mother_cell_data <- fixed_KSHV_data$mother_cell_data

cat(length(unique(daughter_cell_data$mother_cell_id)), "daughter cell pairs with", nrow(daughter_cell_data), "clusters\n")
cat(length(unique(mother_cell_data$cell_id)), "non-dividing cells with", nrow(mother_cell_data), "clusters")

# Estimate number of episomes per cell and Replication and Segregation Efficiency
fixed_KSHV_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder,
                                n_prior = list("geom", 0.5), parallel = T)

### Estimates of Replication and Segregation Efficiency ########################

# The ranges in this table are for joint 95% confidence intervals:
estimates <- fixed_KSHV_results$MLE_grid$estimates %>%  
  pivot_longer(everything(), names_sep = "_", names_to = c("metric", "parameter"))

# Calculate marginal confidence intervals (reported in paper):
# Pr
Pr_marginal <- fixed_KSHV_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

# Ps
Ps_marginal <- fixed_KSHV_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range

# Write out estimates
estimates_full <- estimates %>% 
  filter(metric == "MLE") %>% 
  select(parameter, MLE = value) %>%
  mutate(
    joint_CI = estimates %>% filter(metric != "MLE") %>%  
      group_by(parameter) %>% 
      summarise(joint_CI = paste0("(",paste(value, collapse = ", "), ")")) %>% 
      pull(joint_CI),
    
    marginal_CI = c(paste0("(", paste(Pr_marginal, collapse = ","), ")"),
                    paste0("(", paste(Ps_marginal, collapse = ","), ")"))
  ) 

write_csv(estimates_full, here(results_folder, "MLE_parameter_estimates.csv"))

### Make exploratory plots and paper figures ###################################
make_plots(fixed_KSHV_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- fixed_KSHV_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## Show inference results for example images:
# 12, 9, #10

example_inference <- fixed_KSHV_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("Image 12_") | starts_with("Image 9_") | starts_with("Image #10_") ) %>%  
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  rbind(data.frame(cell_id = "Image 9_2", n_epi = 0)) %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("Image 12_1", "Image 12_2", "Image 9_2", "Image 9_1", "Image #10_1", "Image #10_2"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(0,6)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)

### Post-hoc descriptive statistics ############################################

# Estimates of mu and sigma from MCMC
MCMC_summary <- fixed_KSHV_results$all_chains %>% 
  filter(chain == "chain1") %>% 
  summarise(mu_mean = DescTools::Mode(round(mu)), mu_sd = sd(mu),
            sigma_mean = DescTools::Mode(round(sqrt(1/tau))), sigma_sd = sd(sqrt(1/tau))) %>% 
  mutate(n_daughter_pairs = length(unique(daughter_cell_data$mother_cell_id)),
         n_LANA_dots = nrow(mother_cell_data)+nrow(daughter_cell_data))

write_csv(MCMC_summary, here(results_folder, "MCMC_summary.csv"))

# 1 case when we inferred fewer episomes in a cell than observed visually
fixed_KSHV_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu, tau), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(fixed_KSHV_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(fixed_KSHV_data$mother_cell_data, cluster_id, min_episome_in_cluster))
            ) %>% 
  count(n_epi < min_episome_in_cluster)




## Mean values of intensity per cell and number of episomes per cell
mean_intensity_inference <- tibble(
  
  mean_inferred_epi.daughter = fixed_KSHV_results$daughter_cell_samples %>% 
    filter(chain == "chain1") %>% 
    select(-c(chain, iteration)) %>% 
    as.matrix %>% 
    mean,
  
  mean_inferred_epi.mother = fixed_KSHV_results$mother_cell_samples %>% 
    filter(chain == "chain1") %>% 
    select(-c(chain, iteration)) %>% 
    as.matrix %>% 
    mean,
  
  mean_intensity.daughter = fixed_KSHV_data$daughter_cell_data %>% 
    group_by(cell_id) %>% 
    summarise(intensity = sum(total_cluster_intensity)) %>% 
    pull(intensity) %>% mean,
  
  mean_intensity.mother = fixed_KSHV_data$mother_cell_data %>% 
    group_by(cell_id) %>% 
    summarise(intensity = sum(total_cluster_intensity)) %>% 
    pull(intensity) %>% mean
  
) %>% 
  pivot_longer(everything(), names_sep = "[.]", names_to = c("observation", "set")) %>% 
  pivot_wider(names_from = observation)

write_csv(mean_intensity_inference, here(results_folder, "mean_epi_per_cell.csv") )


