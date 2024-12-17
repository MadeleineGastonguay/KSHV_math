#####
# Simulating a tumor that requires KSHV episomes for viability, such as PEL,
# under hypothetical therapies that reduce either replication or segregation efficiency
# Generate Figure 9, Figure S11, and Figure S12
### Details of stochastic simulations:
# - For each scenario, 5 tumors are simulated stochastically
# - Tumors grow exponentially with a doubling time of 6 days until they reach a size of 1e5 cells. 
# - The tumor starts from one cell latently infected with KSHV, assuming 80% replication 
#   efficiency and 90% segregation efficiency
# - Cells that are created without any episomes die immediately
# - The average lifespan of cells is varied from 1 day to immortal, and the birth rate
#   is calculated so that the tumor doubles in size every 6 days. This birth rate must be faster than 
#   the rate required to observe a doubling time of 6 days in the absence of selection to overcome
#   the additional cell death.
# - Once tumors reach 1e5 cells, either replication efficiency or segregation efficiency is reduced
# - Tumors are simulated for 250 days or until they reach 1e6 cells
#####

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))


### Tumor growth simulations without treatment #################################
## For each d, find b so that observed doubling time is 6 days in the absence of treatment

## Use MLE estimates of replication and segregation efficiency 
pRep <- 0.8
pSeg <- 0.9

## cells with episomes live an average of 1 day
d_1 <- 1/1 # death rate (1/life span; per day)
# Without selection, we expect cells to divide every 1/(log(2)/6 + 1) = 0.9 days for a doubling time of 6 days
b_1 <- log(2)/1.77 + d_1 # birth rate (1/time to division; per day)
1/b_1 # with selection, cells divide after an average of 0.7 days 

# Simulate tumor with selection in the absence of treatment
one_day_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_1, d = d_1, nIts = 1e8,
  nRuns = 30, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes live an average of 2 days
d_2 <- 1/2 
# Without selection, we expect cells to divide every 1.6 days for a doubling time of 6 days
b_2 <- log(2)/2.585 + d_2
1/b_2 # with selection, cells divide after an average of 1.3 days

# Simulate tumor with selection in the absence of treatment
two_days_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_2, d = d_2, nIts = 1e8,
  nRuns = 15, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes live an average of 5 days
d_5 <- 1/5 
# Without selection, we expect cells to divide every 3.2 days for a doubling time of 6 days
b_5 <- log(2)/3.58 + d_5 
1/b_5 # with selection, cells divide after an average of 2.5 days

# Simulate tumor with selection in the absence of treatment
d5_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_5, d = d_5, nIts = 1e8,
  nRuns = 10, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes live an average of 10 days
d_10 <- 1/10 
# Without selection, we expect cells to divide every 4.6 days for a doubling time of 6 days
b_10 <- log(2)/4.1 + d_10 
1/b_10 # with selection, cells divide after an average of 3.7 days

# Simulate tumor with selection in the absence of treatment
d10_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_10, d = d_10, nIts = 1e8,
  nRuns = 5, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes live an average of 20 days
d_20 <- 1/20 
# Without selection, we expect cells to divide every 6 days for a doubling time of 6 days
b_20 <- log(2)/4.41 + d_20 
1/b_20 # with selection, cells divide after an average of 4.8 days

# Simulate tumor with selection in the absence of treatment
d20_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_20, d = d_20, nIts = 1e8,
  nRuns = 5, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes live an average of 40 days
d_40 <- 1/40 
# Without selection, we expect cells to divide every 7.1 days for a doubling time of 6 days
b_40 <- log(2)/4.6 + d_40 
1/b_40 # with selection, cells divide after an average of 5.7 days

# Simulate tumor with selection in the absence of treatment
d40_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_40, d = d_40, nIts = 1e8,
  nRuns = 5, stop_size = 1e5, selection = T, max_epi = 9
)

## cells with episomes never die
d_0 <- 0 
# Without selection, we expect cells to divide every 8.7 days for a doubling time of 6 days
b_0 <- log(2)/4.815 + d_0
1/b_0 # with selection, cells divide after an average of 6.9 days

# Simulate tumor with selection in the absence of treatment
no_death_no_treat_test <- exponential_growth(
  pRep, pSeg, b = b_0, d = d_0, nIts = 1e8,
  nRuns = 5, stop_size = 1e5, selection = T, max_epi = 9
)

## Confirm that all scenarios have a doubling time of ~6 days
no_treatment <- rbind(
  one_day_no_treat_test %>% mutate(d = d_1),
  two_days_no_treat_test %>% mutate(d = d_2),
  d5_no_treat_test %>% mutate(d = d_5),
  d10_no_treat_test %>% mutate(d = d_10),
  d20_no_treat_test %>% mutate(d = d_20),
  d40_no_treat_test %>% mutate(d = d_40),
  no_death_no_treat_test %>% mutate(d = d_0)
)

# write_csv(no_treatment, here("results", "PEL_simulations_with_selection", "tumor_simulations_without_treatment.csv"))
no_treatment <- read_csv(here("results", "PEL_simulations_with_selection", "tumor_simulations_without_treatment.csv"))

no_treatment %>%  
  filter(total == 1e5 | total == 1e5/2) %>% 
  group_by(d, total, run) %>% 
  filter(time == min(time), episomes == 0, pRep == 0.8) %>% 
  distinct %>% 
  group_by(d, run) %>% 
  summarise(dt = diff(time)) %>% 
  ggplot(aes(factor(1/d), dt))  + 
  geom_hline(yintercept = 6, lty = "dashed", alpha = 0.5) + 
  geom_point() + 
  ylim(c(0,7)) + 
  labs(x = "Cell lifespan", y = "Tumor doubling time (days)")


## Compare stochastic simulations to deterministic case
no_treatment %>%   
  group_by(interaction(run, d)) %>%
  filter(max(total) >= 1e4) %>% 
  # Standardize simulations to start at time = 0 when tumor size is 1e3
  mutate(t = time[min(which(total == 1e4))], time = time - t) %>%
  ggplot(aes(time, total, color = factor(d), group = interaction(run,d))) + 
  # Plot stochastic simulations
  geom_line() + 
  # Calculate and plot deterministic simulation with a doubling time of 6 days
  geom_line(data = tibble(time = seq(0, 100, by = 0.1), total = exp(log(2)/6*time)), 
            aes(time - log(1e4)*6/log(2), total), inherit.aes = F, alpha = 0.75) +
  scale_y_log10(limits = c(1e4, 1e5)) +
  labs(x = "Time from 1000 cells", y = "Tumor size", color = "Death rate",
       title = "PEL growth before treatment",
       subtitle = "Goal doubling time = 6 days") +
  xlim(c(0, 20))   + 
  facet_wrap(~1/d)


### Simulate a theoretical intervention reducing replication efficiency ########

## Reduced replication efficiencies to simulate:
pRep_reduced <- seq(0.1, 0.5, by = 0.1)

## wrapper function to simulate for a tumor given reduced pRep values, birth, and death rate
sim_reduced_rep <- function(pRep, pSeg, pRep_reduced,  d, b, combine = NULL){
  PEL_simulations(
    pRep, pSeg, # Replication and segregation efficiency before treatment is applied
    c(pRep_reduced, pRep), pSeg, # Replication and segregation efficiency after treatment is applied
    b = b, d = d, # Birth and death rate
    selection = T, # Simulate with selection
    treatment_size = 1e5,  # apply treatment when tumor has 1e5 cells
    max_epi = 9,  # Cap at 9 episomes per cell for book-keeping
    nRuns = 5, # simulate 5 tumors
    stop_time = 250, stop_size = 1e6, # Stop simulation when tumor reaches 1e6 cells or 250 days
    add_to = combine # Results to combine with if running additional scenarios
  )
}

## Cell lifespan of 1 day:
d1_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_1, b_1) 

# Additional replication thresholds to test based on first set of simulations
d1_pRep_reduction <- sim_reduced_rep(pRep, pSeg, c(0.6, 0.7, 0.75), d_1, b_1, combine = d1_pRep_reduction) 

## Cell lifespan of 2 days:
d2_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_2, b_2) 

# Additional replication thresholds to test based on first set of simulations
d2_pRep_reduction <- sim_reduced_rep(pRep, pSeg, c(0.6, 0.65, 0.7), d_2, b_2, combine = d2_pRep_reduction) 
# d2_pRep_reduction <- sim_reduced_rep(pRep, pSeg, 0.65, d_2, b_2, combine = d2_pRep_reduction) 

## Cell lifespan of 5 days:
d5_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_5, b_5) 

# Additional replication thresholds to test based on first set of simulations
d5_pRep_reduction <- sim_reduced_rep(pRep, pSeg, c(0.55, 0.6), d_5, b_5, combine = d5_pRep_reduction) 

## Cell lifespan of 10 days:
d10_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_10, b_10) 

# Additional replication thresholds to test based on first set of simulations
d10_pRep_reduction <- sim_reduced_rep(pRep, pSeg, 0.35, d_10, b_10, combine = d10_pRep_reduction) 

## Cell lifespan of 20 days:
d20_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_20, b_20) 

# Additional replication thresholds to test based on first set of simulations
d20_pRep_reduction <- sim_reduced_rep(pRep, pSeg, 0.25, d_20, b_20, combine = d20_pRep_reduction) 

## Cell lifespan of 40 days:
d40_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_40, b_40) 

# Additional replication thresholds to test based on first set of simulations
d40_pRep_reduction <- sim_reduced_rep(pRep, pSeg, 0.15, d_40, b_40, combine = d40_pRep_reduction) 

## No death of cells with episomes:
no_death_pRep_reduction <- sim_reduced_rep(pRep, pSeg, pRep_reduced, d_0, b_0) 

# Additional replication thresholds to test based on first set of simulations
no_death_pRep_reduction <- sim_reduced_rep(pRep, pSeg, 0, d_0, b_0, combine = no_death_pRep_reduction) 


### Plot simulation results for replication reduction ##########################
# Compare results
all_results <- rbind(
  d1_pRep_reduction %>% mutate(d = d_1),
  d2_pRep_reduction %>% mutate(d = d_2),
  d5_pRep_reduction %>% mutate(d = d_5),
  d10_pRep_reduction %>% mutate(d = d_10),
  d20_pRep_reduction %>% mutate(d = d_20),
  d40_pRep_reduction %>% mutate(d = d_40),
  no_death_pRep_reduction %>% mutate(d = d_0)
  ) %>% 
  mutate(lifespan = factor(1/d, levels = sort(unique(1/.$d)), 
                           labels = paste("Baseline cell lifespan:", sort(unique(1/.$d)), "days")),
         pRep = factor(pRep, levels = rev(sort(unique(.$pRep))), 
                       labels = c("Baseline", paste0(rev(sort(unique(.$pRep)))[-1]*100, "%")))) 

write_csv(all_results, file = here("results", "PEL_simulations_with_selection", "PEL_replication_threshold_simulations.csv"))

## Read in results from previous runs
# all_results <- read_csv(here("results", "PEL_simulations_with_selection", "PEL_replication_threshold_simulations.csv")) 
# all_results <- all_results %>% mutate(
#   pRep = factor(pRep, levels = rev(sort(unique(.$pRep)))),
#   lifespan = fct_inorder(lifespan)
# )


## Figure S11
PEL_threshold_simulations <- all_results %>% 
  group_by(run, d) %>% 
  mutate(t = min(time[which(total == 1e5)]),
         time = time - t) %>% 
  ggplot(aes(time, total, group = interaction(run, pRep, d), color = factor(pRep))) + 
  geom_line(alpha = 0.5) +
  ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pRep) %>%
                             # filter(time <= 200, total >= 1e3, pRep != "Baseline") %>%
                             filter(time <= 125, total >= 1e3, pRep != "Baseline") %>%
                             filter(time == max(time), episomes == 0) %>% distinct,
                           aes(x = time, y = total, label = pRep),
                           show.legend = FALSE, fontface = 2, nudge_x = 0.25, nudge_y = 0.25, #hjust = 0,  , 
                           min.segment.length = 0, box.padding = 0.25, bg.color = "white") +
  scale_colour_scico_d(palette = "acton", end = 0.9) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  scale_y_log10(breaks = c(10^(3:6)), labels = expression(10^3, 10^4, 10^5, 10^6)) +
  coord_cartesian(ylim = c(1e3, 1.3e6), xlim = c(-50, 150)) + 
  facet_wrap(~lifespan) +
  labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
       color = "Replication\nEfficiency (%)")  + 
  theme(legend.position = "none")

ggsave(here("results", "PEL_simulations_with_selection", "PEL_replication_threshold_simulations_full.pdf"), 
       PEL_threshold_simulations,
       width = 12, height = 10)


### Simulate a theoretical intervention reducing segregation efficiency ########

## Reduced segregation efficiencies to simulate:
pSeg_reduced <- seq(0, 0.5, by = 0.1)

## wrapper function to simulate for a tumor given reduced pRep values, birth, and death rate
sim_reduced_seg <- function(pRep, pSeg, pSeg_reduced,  d, b, combine = NULL){
  PEL_simulations(
    pRep, pSeg, # Replication and segregation efficiency before treatment is applied
    pRep, c(pSeg_reduced, pSeg), # Replication and segregation efficiency after treatment is applied
    b = b, d = d, # Birth and death rate
    selection = T, # Simulate with selection
    treatment_size = 1e5,  # apply treatment when tumor has 1e5 cells
    max_epi = 9,  # Cap at 9 episomes per cell for book-keeping
    nRuns = 5, # simulate 5 tumors
    stop_time = 250, stop_size = 1e6, # Stop simulation when tumor reaches 1e6 cells or 250 days
    add_to = combine # Results to combine with if running additional scenarios
  )
}

## Cell lifespan of 1 day:
d1_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_1, b_1) 

## Cell lifespan of 2 days:
d2_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_2, b_2) 

## Cell lifespan of 5 days:
d5_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_5, b_5) 

## Cell lifespan of 10 days:
d10_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_10, b_10) 

## Cell lifespan of 20 days:
d20_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_20, b_20) 

## Cell lifespan of 40 days:
d40_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_40, b_40) 

## No death of cells with episomes:
no_death_pSeg_reduction <- sim_reduced_seg(pRep, pSeg, pSeg_reduced, d_0, b_0) 

### Plot results ###############################################################

all_results_seg <- rbind(
  d1_pSeg_reduction %>% mutate(d = d_1),
  d2_pSeg_reduction %>% mutate(d = d_2),
  d5_pSeg_reduction %>% mutate(d = d_5),
  d10_pSeg_reduction %>% mutate(d = d_10),
  d20_pSeg_reduction %>% mutate(d = d_20),
  d40_pSeg_reduction %>% mutate(d = d_40),
  no_death_pSeg_reduction %>% mutate(d = d_0)
) %>% 
  mutate(lifespan = factor(1/d, levels = sort(unique(1/.$d)), 
                           labels = paste("Baseline cell lifespan:", sort(unique(1/.$d)), "days")),
         pSeg = factor(pSeg, levels = rev(sort(unique(.$pSeg))), 
                       labels = c("Baseline", paste0(rev(sort(unique(.$pSeg)))[-1]*100, "%")))) 

write_csv(all_results_seg, file = here("results", "PEL_simulations_with_selection", "PEL_segregation_threshold_simulations.csv"))
# all_results_seg <- read_csv(here("results", "PEL_simulations_with_selection", "PEL_segregation_threshold_simulations.csv"))
# all_results_seg <- all_results_seg %>% mutate(
#   pRep = factor(pSeg, levels = rev(sort(unique(.$pSeg)))),
#   lifespan = fct_inorder(lifespan)
# )

## Figure S12
PEL_threshold_simulations_seg <- all_results_seg %>%
  # filter(run == 1) %>%
  group_by(run, d) %>% 
  mutate(t = min(time[which(total == 1e5)]),
         time = time - t) %>% 
  ggplot(aes(time, total, group = interaction(run, pSeg, d), color = pSeg)) + 
  geom_line(alpha = 0.5) +
  ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pSeg) %>%
                             filter(pSeg != "Baseline") %>% 
                             # filter(time <= 200, total >= 1e3) %>%
                             mutate(t = abs(time - 15)) %>% 
                             filter(t == min(t), episomes == 0) %>% 
                             arrange(pSeg),
                           aes(x = time, y = total, label = pSeg),
                           direction = "y", nudge_x = 30,
                           show.legend = FALSE, fontface = 2, box.padding = 0.25,
                           segment.size = 0.25, size = 3, bg.color = "white") +
  scale_colour_scico_d(palette = "oslo", end = 0.9) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  scale_y_log10(breaks = c(10^(3:6)), labels = expression(10^3, 10^4, 10^5, 10^6)) +
  coord_cartesian(ylim = c(1e3, 1e6), xlim = c(-50, 50)) +
  facet_wrap(~lifespan) +
  labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
       color = "Segregation\nEfficiency (%)") +
  theme(legend.position = "none")

ggsave(here("results", "PEL_simulations_with_selection", "PEL_segregation_threshold_simulations_full.pdf"), 
       PEL_threshold_simulations_seg,
       width = 12, height = 10)
       #width = 8, height = 7)

### Figure for paper with average cell life span of 1 day ######################

## Figure 9
PEL_paper_figure <- all_results %>% filter(d == d_1) %>% 
  group_by(run, d) %>% 
  mutate(t = min(time[which(total == 1e5)]),
         time = time - t,
         pRep = fct_recode(pRep, `Baseline (80%)` = "Baseline")) %>% 
  ggplot(aes(time, total,  color = factor(pRep))) + 
  geom_line(alpha = 0.5, aes(group = interaction(run, pRep, d))) +
  ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pRep) %>%
                             filter(time <= 60, total >= 1e3) %>%
                             summarise(t = min(max(time), 60),
                                       total = max(total[time == t]), time = t) %>% 
                             # filter(time == max(time), episomes == 0) %>% distinct %>% 
                             arrange(pRep),
                           aes(x = time, y = total, label = pRep),
                           show.legend = FALSE, fontface = 2, nudge_x = c(-10, 10, 10, 10, 15, 10, -15, -10, -15), nudge_y = 0.1, #hjust = 0,  ,
                           min.segment.length = 0, box.padding = 0.25, size = 3,  segment.size = 0.25, bg.color = "white") +
  scale_colour_scico_d(palette = "acton", end = 0.9) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  coord_cartesian(ylim = c(1e3, 1e6), xlim = c(-50, 100)) + 
  labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
       color = "Replication\nEfficiency (%)", title = "Reduction in Replication Efficiency") +
  
  plot_spacer() +
  
  all_results_seg %>%
  filter(d == d_1) %>% 
  group_by(run, d) %>% 
  mutate(t = min(time[which(total == 1e5)]),
         time = time - t,
         pSeg = fct_recode(pSeg, `Baseline (90%)` = "Baseline")) %>% 
  ggplot(aes(time, total, group = interaction(run, pSeg, d), color = pSeg)) + 
  geom_line(alpha = 0.5) +
  ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pSeg) %>%
                             # filter(time <= 200, total >= 1e3) %>%
                             mutate(t = abs(time - 15)) %>% 
                             filter(t == min(t), episomes == 0) %>% 
                             arrange(pSeg),
                           aes(x = time, y = total, label = pSeg),
                           direction = "y", nudge_x = c(-20, rep(20, 6)),
                           show.legend = FALSE, fontface = 2, box.padding = 0.25,
                           segment.size = 0.25, size = 3, bg.color = "white") +
  scale_colour_scico_d(palette = "oslo", end = 0.9) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  coord_cartesian(ylim = c(1e3, 1e6), xlim = c(-50, 40)) +
  labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
       color = "Segregation\nEfficiency (%)", title = "Reduction in Segregation Efficiency")  + 
  
  plot_layout(nrow = 1, widths = c(1, 0.1, 1)) & 
  scale_y_log10(labels = expression(10^3, 10^4, 10^5, 10^6), breaks = 10^c(3:6)) &
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(here("results", "PEL_simulations_with_selection", "PEL_simulations.pdf"),
       PEL_paper_figure, width = 9, height = 4.5)

### Figure for MIDAS poster ####################################################

rep_poster <- all_results %>% filter(d == 1/5) %>% 
   group_by(run, d) %>% 
   mutate(t = min(time[which(total == 1e5)]),
          time = time - t) %>% 
   mutate(pRep = fct_recode(pRep, `Baseline (80%)` = "Baseline")) %>% 
   ggplot(aes(time, total, group = interaction(run, pRep, d), color = factor(pRep))) + 
   geom_line(alpha = 0.5) +
   ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pRep) %>%
                              filter(time <= 200, total >= 1e3) %>%
                              filter(time == max(time), episomes == 0) %>% distinct,
                            aes(x = time, y = total, label = pRep),
                            show.legend = FALSE, fontface = 2, nudge_x = 0.25, #hjust = 0,  ,
                            min.segment.length = 0, box.padding = 0.25, size = 7) +
   scale_colour_scico_d(palette = "acton", end = 0.9) +
   guides(color = guide_legend(override.aes = list(alpha = 1))) + 
   # scale_y_log10() +
   coord_cartesian(ylim = c(1e3, 1e6), xlim = c(-50, 250)) + 
   # facet_wrap(~lifespan) +
   labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
        color = "Replication\nEfficiency (%)", title = "Reduction in Replication Efficiency") + 
   theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank(),
         plot.background = element_blank(), legend.box.background = element_blank()) + 
   theme(legend.position = "none") +
   scale_y_log10(labels = expression(10^3, 10^4, 10^5, 10^6), breaks = 10^c(3:6)) 

ggsave(here(out_folder, "MIDAS_poster_rep.png"), rep_poster, bg = "transparent",
       width = 8, height = 7)

seg_poster <- all_results_seg %>%
  filter(d == 1/5) %>% 
  group_by(run, d) %>% 
  mutate(t = min(time[which(total == 1e5)]),
         time = time - t) %>% 
  mutate(pSeg = fct_recode(pSeg, `Baseline (90%)` = "Baseline")) %>% 
  ggplot(aes(time, total, group = interaction(run, pSeg, d), color = pSeg)) + 
  geom_line(alpha = 0.5) +
  # ggrepel::geom_text_repel(data = . %>% group_by(lifespan, pSeg) %>%
  #                            # filter(time <= 200, total >= 1e3) %>%
  #                            filter(time == max(time), episomes == 0) %>% distinct,
  #                          aes(x = time, y = total, label = pSeg),
  #                          show.legend = FALSE, fontface = 2, nudge_y = 0.3, #hjust = 0,  ,
  #                          min.segment.length = 0, box.padding = 0.25, size = 7) +
  scale_colour_scico_d(palette = "oslo", end = 0.9) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  # scale_y_log10() +
  coord_cartesian(ylim = c(1e3, 1e6), xlim = c(-50, 40)) +
  # facet_wrap(~lifespan) +
  labs(x = "Time from treatment (days)", y = "Tumor Size (cells)", 
       color = "Segregation\nEfficiency (%)", title = "Reduction in Segregation Efficiency")  + 
  theme(plot.title = element_text(hjust = 0.5), legend.background = element_rect(),
        plot.background = element_blank(), legend.box.background = element_blank()) + 
  theme(legend.position = c(0.95,0.05), legend.justification = c(1,0)) +
  scale_y_log10(labels = expression(10^3, 10^4, 10^5, 10^6), breaks = 10^c(3:6)) 

ggsave(here(out_folder, "MIDAS_poster_seg.png"), seg_poster, bg = "transparent",
       width = 8, height = 7)


