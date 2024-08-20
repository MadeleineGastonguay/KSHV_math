# This script runs the full pipeline for data from fixed images of cells with only KHSV terminal repeats (Figure 2)

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

set.seed(1234)


#####
# Script inputs

# Folder for results
results_folder <- here("results", "live_KSHV")

# Daughter cell intensity data
daughter_cell_file <- "live_KSHV_dividing_cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "live_KSHV_non_dividing_cells.xlsx"

#####

live_KSHV_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- live_KSHV_data$daughter_cell_data
mother_cell_data <- live_KSHV_data$mother_cell_data

live_KSHV_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, same_mu = F,
                                n_prior = list("geom", 0.5), parallel = T)


live_KSHV_results$MLE_grid$estimates


## Adjust confidence intervals to be marginal:

live_KSHV_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

live_KSHV_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range

make_plots(live_KSHV_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- live_KSHV_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## How many cases did we infer fewer episomes than observed visually?
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



## Mean values of intensity per cell and number of episomes per cell
live_KSHV_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

live_KSHV_results$mother_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

live_KSHV_data$daughter_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

live_KSHV_data$mother_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

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


## Get standard deviation of inferred mean and sigma:
live_KSHV_results$all_chains %>% 
  filter(chain == "chain1") %>% 
  select(mu_d, tau_d, mu_m, tau_m) %>% 
  mutate(sd_d = sqrt(1/tau_d), sd_m = sqrt(1/tau_m)) %>% 
  apply(2, sd)

