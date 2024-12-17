#####
# Simulating an exponentially growing cell population under selection
# Generate Figure 8
#####

### Setup ###################################################################### 
library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
library(Hmisc)
library(ggrepel)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))

### Script inputs ##############################################################
## Run 100 trials of simulations with varying parameters
ntrials <- 100

## Number of cells to start with
n_cells <- 1

## Number of episomes per cell at the start
n_epi <- 3

## Indicator if simulations should be rerun (TRUE) or loaded from previous run (FALSE)
rerun <- FALSE

## create out folder
out_folder <- here("results","simulations_expo_selection")
if(!file.exists(out_folder)) dir.create(out_folder, showWarnings = F)

if(!rerun){
  # Load results from previous run if they exist
  load(here(out_folder, "expo_pop_3epi.RData"))
}else{
  # Run simulations
  
  ## set up a grid of parameters to simulate with and collect all results:
  pReps = 0.8
  pSegs = c(0.8, 0.9)
  param_grid <- expand_grid(pRep = pReps, pSeg = pSegs, n_cells_start = n_cells, n_epi_start = n_epi)
  
  expo_3epi_df <- param_grid %>% 
    pmap_df(function(pRep, pSeg, n_cells_start, n_epi_start) {
      result <- exponential_growth(pRep, pSeg, nIts = 150000, nRuns = 100, n_cells_start, n_epi_start, selection = TRUE, max_epi = 9)
      cbind(data.frame(pRep = pRep, pSeg = pSeg, n_cells_start = n_cells_start), result)
    }) %>% 
    rename(Pr = pRep, Ps = pSeg, episomes_per_cell = episomes, trial = run, n_cells = n_cells_start) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps)),
           number_of_cells = frac*total)
  
  # Calculate average number of episomes per cell as population size increases
  averages <- expo_3epi_df %>% 
    group_by(total, Pr, Ps, trial) %>% 
    filter(episomes_per_cell != -1) %>% 
    summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
    ungroup()
  
  save(expo_3epi_df, param_grid, averages, 
       file = here(out_folder, "expo_pop_3epi.RData"))
}

### Figure 8 ###################################################################
# Plot average number of episomes per cell and fraction of population with k-episomes
# as the cell population increases

figure_plot <- expo_3epi_df %>% 
  filter(episomes_per_cell == -1, Ps == "80%") %>% 
  ggplot(aes(total, frac, group = interaction(Pr, Ps))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", color = "black") + 
  labs(x = "Population Size", y = "Average number of episomes per cell", color = "Segregation Efficiency",
       title = "80% Replication, 80% Segregation")  +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) + 
  
  expo_3epi_df %>% 
  filter(episomes_per_cell == -1, Ps == "90%") %>% 
  ggplot(aes(total, frac, group = interaction(Pr, Ps))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", color = "black") + 
  labs(x = "Population Size", y = "Average number of episomes per cell", color = "Segregation Efficiency",
       title = "80% Replication, 90% Segregation")  +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) +
  
  expo_3epi_df %>% 
  filter(episomes_per_cell != -1, Ps == "80%") %>% 
  ggplot(aes(total, frac, color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  # Add labels directly to lines:
  geom_text_repel(data = ~filter(., total == max(total), episomes_per_cell <= 6, episomes_per_cell > 0) %>% 
                    group_by(episomes_per_cell)  %>% summarise(frac = mean(frac)),
                  aes(1e5, frac, label = episomes_per_cell), show.legend = F,
                  nudge_x = 0.25, nudge_y = 0.01, min.segment.length = 0, box.padding = 0.25, bg.color = "white") +
  labs(x = "Population Size", y = "Relative Abundance", color = "Episomes\nper cell") + 

  
  expo_3epi_df %>% 
  filter(episomes_per_cell != -1, Ps == "90%") %>% 
  ggplot(aes(total, frac, color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  # Add labels directly to lines:
  geom_text_repel(data = ~filter(., total == max(total), episomes_per_cell <= 5, episomes_per_cell > 0) %>% 
                    group_by(episomes_per_cell)  %>% summarise(frac = mean(frac)),
                  aes(1e5, frac, label = episomes_per_cell), show.legend = F,
                  nudge_x = 0.25, nudge_y = 0.01, min.segment.length = 0, box.padding = 0.25, bg.color = "white") +
  labs(x = "Population Size", y = "Relative Abundance", color = "Episomes\nper cell") + 
  
  plot_layout(guides = "collect") &
  scale_color_manual(values = safe_colorblind_palette) & 
  theme(legend.justification = c(0.5, -0.05)) &
  scale_x_log10(breaks = 10^(0:5), labels = expression(1, 10, 10^2, 10^3, 10^4, 10^5))

ggsave(here(out_folder, "exponential_selection_figure.pdf"), figure_plot, width = 10, height = 7)


