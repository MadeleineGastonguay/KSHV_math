# This script runs the full pipeline for data from fixed images of cells with only KHSV terminal repeats (Figure 2)

#####
#  Setup 
#####

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

set.seed(1234)

#####
# Script inputs

# Folder for results
results_folder <- here("results", "fixed_8TR")

# Daughter cell intensity data
daughter_cell_file <- "fixed_8TR_dividing_cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "fixed_8TR_non_dividing_cells.xlsx"

#####

fixed_8TR_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- fixed_8TR_data$daughter_cell_data
mother_cell_data <- fixed_8TR_data$mother_cell_data

cat(length(unique(daughter_cell_data$mother_cell_id)), "daughter cell pairs with", nrow(daughter_cell_data), "clusters\n")
cat(length(unique(mother_cell_data$cell_id)), "non-dividing cells with", nrow(mother_cell_data), "clusters")

fixed_8TR_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, 
                                n_prior = list("geom", 0.5), parallel = T, just_Pr = F)

fixed_8TR_results$MLE_grid$estimates

## Adjust confidence intervals to be marginal:

fixed_8TR_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

fixed_8TR_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range


make_plots(fixed_8TR_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- fixed_8TR_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## How many cases did we infer fewer episomes than observed visually?
fixed_8TR_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu, tau), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(fixed_8TR_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(fixed_8TR_data$mother_cell_data, cluster_id, min_episome_in_cluster))
  ) %>% 
  count(n_epi < min_episome_in_cluster)

## Show inference results for example images:
# 30, 31, 3, 4
example_inference <- fixed_8TR_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("30") | starts_with("31") | starts_with("3_") | starts_with("4_")) %>% 
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("30_2", "30_1", "31_1", "31_2", "3_2", "3_1", "4_1"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,8, by = 2), limits = c(0,8)) +
  # ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)

## Mean values of intensity per cell and number of episomes per cell
fixed_8TR_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

fixed_8TR_results$mother_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

fixed_8TR_data$daughter_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

fixed_8TR_data$mother_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

## Get standard deviation of inferred mean and sigma:

fixed_8TR_results$all_chains %>% 
  filter(chain == "chain1") %>% 
  select(mu, tau) %>% 
  mutate(sd = sqrt(1/tau)) %>% 
  apply(2, sd)


