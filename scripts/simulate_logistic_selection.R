## 
# Logistic growth with selection
# Simulating a tumor with 1000 cells that require KSHV epsiomes for viability
##

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))

out_folder <- here("results", "logistic_simulations_with_selection")
if(!dir.exists(out_folder)) dir.create(out_folder)

# values chosen based on preliminary simulations that show segregation does not play a roll
pReps <- c(0.1, 0.3, 0.4, 0.5, 0.8, 0.9)
pSegs <- c(0.1, 0.8, 0.9)
param_grid <- expand_grid(pRep = pReps, pSeg = pSegs)

# Use 1000 cells: preliminary simulations show that threshold does not depend on this number,
# but time until tumor is eliminated does
n_cells <- 1000

PEL_selection <- param_grid %>% 
  pmap(extinction, nTrials = 100, n_epi = 3, selectAgainstZero = T, n_cells = n_cells)

PEL_selection_df <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], PEL_selection[[i]]$Totals)) 

PEL_selection_with_total <- PEL_selection_df %>% 
  mutate(tumor_size = PEL_selection_df %>% select(6:15) %>% apply(1, sum)) %>% 
  pivot_longer(6:ncol(.)) %>% 
  mutate(name = factor(name, levels = c("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "tumor_size")))


supplemental_figure1 <- PEL_selection_with_total %>%
  ggplot(aes(time, value, group = interaction(trial,name), color = name)) +
  geom_line(alpha = 0.5) +
  scale_color_manual(values = c(safe_colorblind_palette[1:10], "darkgray")) +
  ggh4x::facet_grid2(pRep ~ pSeg, scales = "free_x", independent = "x", labeller = "label_both")

ggsave(here(out_folder, "sup_fig1.pdf"), supplemental_figure1, width = 20, height = 20)

supplemental_figure2 <- PEL_selection_with_total %>%
  filter(name == "tumor_size") %>% 
  group_by(time = round(time, 1), Pr = pRep*100, Ps = pSeg*100, trial, name) %>% 
  summarise(value = mean(value)) %>% 
  ggplot(aes(time, value)) +
  stat_summary(fun.data = "mean_cl", geom = "smooth") +
  ggh4x::facet_grid2(pRep ~ pSeg, scales = "free_x", independent = "x", labeller = "label_both")

ggsave(here(out_folder, "sup_fig2.pdf"), supplemental_figure2, width = 20, height = 20)

#####
# Make these plots:
# Fraction of populations that go extinct at each combination

# Size of the population as a fucntion of model parameters for those that don't go extinct
PEL_selection_with_total %>%
  filter(name == "tumor_size") %>% 
  group_by(pRep = pRep*100, pSeg = pSeg*100, trial) %>%
  filter(time >= ifelse(min(value) == 0, max(time), 500)) %>%
  summarise(value = mean(value)) %>% 
  ggplot(aes(pRep, value)) +
  geom_point(aes(color = factor(pSeg))) +
  labs(x = "Replication Efficiency (%)", y = "Steady State Population Size", color = "Ps (%)")
