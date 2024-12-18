# This script runs the full pipeline for data from images of live cells with 
# full KSHV (figure 6)

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
results_folder <- here("results", "live_KSHV")

# Daughter cell intensity data
daughter_cell_file <- "live_KSHV_dividing_cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "live_KSHV_non_dividing_cells.xlsx"

### Run analysis ###############################################################

# Load in the data
live_KSHV_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- live_KSHV_data$daughter_cell_data
mother_cell_data <- live_KSHV_data$mother_cell_data

cat(length(unique(daughter_cell_data$mother_cell_id)), "daughter cell pairs with", nrow(daughter_cell_data), "clusters\n")
cat(length(unique(mother_cell_data$cell_id)), "non-dividing cells with", nrow(mother_cell_data), "clusters")

# Estimate number of episomes per cell and Replication and Segregation Efficiency
live_KSHV_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, same_mu = F,
                                n_prior = list("geom", 0.5), parallel = T)


### Estimates of Replication and Segregation Efficiency ########################

# The ranges in this table are for joint 95% confidence intervals:
estimates <- live_KSHV_results$MLE_grid$estimates %>%  
  pivot_longer(everything(), names_sep = "_", names_to = c("metric", "parameter"))

# Calculate marginal confidence intervals (reported in paper):
# Pr
Pr_marginal <- live_KSHV_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

# Ps
Ps_marginal <- live_KSHV_results$MLE_grid$grid_search %>%
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
make_plots(live_KSHV_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- live_KSHV_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## Show inference results for example images:
# 16, 14
example_inference <- live_KSHV_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("Cell 16") | starts_with("Cell 14")) %>% 
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  add_row(cell_id = "Cell 14_2", n_epi = 0) %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = paste("Cell", c("16_1", "16_2", "14_2", "14_1")))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, ncol = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,8, by = 2), limits = c(0,8)) +
  # ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank(), plot.margin = margin(5, 12, 5, 5)) +
  coord_flip()

ggsave(here(results_folder, "example_inference.png"), width = 1.5, height = 5)


### Post-hoc descriptive statistics ############################################

# Estimates of mu and sigma from MCMC
MCMC_summary <- live_KSHV_results$all_chains %>% 
  filter(chain == "chain1") %>% 
  summarise(mu_d_mean = DescTools::Mode(round(mu_d)), mu_d_sd = sd(mu_d),
            sigma_d_mean = DescTools::Mode(round(sqrt(1/tau_d))), sigma_d_sd = sd(sqrt(1/tau_d)),
            mu_m_mean = DescTools::Mode(round(mu_m)), mu_m_sd = sd(mu_m),
            sigma_m_mean = DescTools::Mode(round(sqrt(1/tau_m))), sigma_m_sd = sd(sqrt(1/tau_m))) %>% 
  mutate(n_daughter_pairs = length(unique(daughter_cell_data$mother_cell_id)),
         n_LANA_dots_d = nrow(daughter_cell_data),
         n_LANA_dots_m = nrow(mother_cell_data))

write_csv(MCMC_summary, here(results_folder, "MCMC_summary.csv"))


# We never infer fewer episomes than observed visually
live_KSHV_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu_d, tau_d, mu_m, tau_m), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(live_KSHV_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(live_KSHV_data$mother_cell_data, cluster_id, min_episome_in_cluster))
  ) %>% 
  count(n_epi < min_episome_in_cluster)


# Mean values of intensity per cell and number of episomes per cell
mean_intensity_inference <- tibble(
  
  mean_inferred_epi.daughter = live_KSHV_results$daughter_cell_samples %>% 
    filter(chain == "chain1") %>% 
    select(-c(chain, iteration)) %>% 
    as.matrix %>% 
    mean,
  
  mean_inferred_epi.mother = live_KSHV_results$mother_cell_samples %>% 
    filter(chain == "chain1") %>% 
    select(-c(chain, iteration)) %>% 
    as.matrix %>% 
    mean,
  
  mean_intensity.daughter = live_KSHV_data$daughter_cell_data %>% 
    group_by(cell_id) %>% 
    summarise(intensity = sum(total_cluster_intensity)) %>% 
    pull(intensity) %>% mean,
  
  mean_intensity.mother = live_KSHV_data$mother_cell_data %>% 
    group_by(cell_id) %>% 
    summarise(intensity = sum(total_cluster_intensity)) %>% 
    pull(intensity) %>% mean
  
) %>% 
  pivot_longer(everything(), names_sep = "[.]", names_to = c("observation", "set")) %>% 
  pivot_wider(names_from = observation)

write_csv(mean_intensity_inference, here(results_folder, "mean_epi_per_cell.csv") )



