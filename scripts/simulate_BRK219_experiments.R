#####
# Simulating an exponentially growing high copy number KSHV-infected cell population to emulate BRK219 experiments
### Details of stochastic simulations:
# For each passage, we fit a birth-death model to cell growth data to estimate the birth rate (b) assuming no cell death (d = 0)
# We initialize the number of episomes per cell based on the distribution of LANA dots per cell at day 0
# We simulate cell growth using the estimated birth rate for each time period using the 
# SUM159-informed estimates of 80% replication efficiency and 90% segregation efficiency
# We calculate the percent of cells with episomes, the average number of episomes per cell, and teh distribution of episome copy number per cell over time
#####


### Setup ######################################################################
library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
library(ggrepel)
library(ggthemes)
library(scales)
library(fitdistrplus)
library(ggdist)
library(ggnewscale)
theme_set(theme_minimal())

source(here("scripts", "functions_simulations.R"))
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

out_folder <- here("results", "brk219")
if(!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

### Inputs #####################################################################
LANA_dots <- read_csv(here("data", "derived", "brk219_full_LANA_dots.csv"))
cell_growth <- read_csv(here("data", "derived", "brk219_cell_growth.csv"))

# summarise mean LANA dots
LANA_summary <- LANA_dots %>% 
  group_by(day) %>% 
  summarise(mean = mean(LANA_dots),
            median = median(LANA_dots),
            var = var(LANA_dots),
            sd = sd(LANA_dots),
            n_zero = sum(LANA_dots == 0),
            total_dots = sum(LANA_dots),
            sem = sd/sqrt(n()))

### Plot distribution of LANA dots per day #####################################
LANA_boxplots <- LANA_dots %>% 
  ggplot(aes(factor(day), LANA_dots)) + 
  geom_boxplot(outlier.shape = NA, fill = "gray", alpha = 0.7) + 
  ggbeeswarm::geom_beeswarm(size = 1, alpha = 0.7, method = "center") +
  labs(x = "Day", y = "LANA dots per cell") +
  theme(panel.grid.major.x = element_blank())

ggsave(here(out_folder, "LANA_boxplots.pdf"), LANA_boxplots, width = 5, height = 4)

### Define variation in initial episomes based on distribution of LANA dots at day 0 ####
day0_LANA <- LANA_dots %>% filter(day == 0) %>% pull(LANA_dots)

# Negative binomial fits better than poisson because variance > mean
fit_nb_day0 <- fitdist(day0_LANA, "nbinom")
fit_nb_day0$estimate
# size       mu 
# 11.34859 24.78032 
fit_confint <- confint(fit_nb_day0)
#         2.5 %   97.5 %
# size  3.948792 18.74839
# mu   22.061602 27.49903

### Fit variable rates of cell growth to data ##################################
cut_times <- c(4.5,9.5,16,21.5,26.5,31.5,38.5,43.5) # Passaging times
no_growth <- c(14,25) # Times without cell grwoth

# Function to simulate cell growth:
simulate_cell_growth_variable <- function(rs, cut_times, no_growth){
  # cut_times and no-growth times are fixed, fitting r for each growth period
  # No-growth indicates the start of no growth, it is assumed that growth stops until the next cut time
  sim_growth_full <- NULL
  for(i in 1:length(cut_times)){
    t_start <- ifelse(i == 1, 0, cut_times[i-1])
    # t_end   <- ifelse(i == 1, cut_times[i], cut_times[i] - cut_times[i-1])
    t_end <- cut_times[i] - t_start - 0.1
    
    if(length(no_growth)>0 & no_growth[1] < cut_times[i]){
      stop_grow <- no_growth[1] - t_start
      no_growth <- no_growth[-1]
      day2 <- seq(stop_grow + 0.1, t_end, by = 0.1)
    }else{
      stop_grow <- Inf
      day2 <- NULL
    }
    
    r <- rs[i]
    day <- seq(0, min(t_end, stop_grow), by = 0.1)
    
    sim_growth1 <- 1e5*exp(r*day)
    sim_growth2 <- sim_growth1[length(sim_growth1)]*exp(0*day2)
    sim_growth_full <- rbind(sim_growth_full, data.frame(day = c(day,day2) + t_start, cells = c(sim_growth1, sim_growth2)))
  }
  
  return(sim_growth_full)
}


# Cost function to fit to data:
cost <- function(pars){
  no_growth <- c(14,25)
  cut_times <- c(4.5,9.5,16,21.5,26.5,31.5,38.5,43.5)
  sim <- simulate_cell_growth_variable(pars, cut_times, no_growth) 
  data_days <- cell_growth$day
  idx <- data_days < 43.5
  data_days <- data_days[idx]
  sim_at_time <- approx(sim$day, sim$cells, data_days)$y/1e5
  diff <- cell_growth$live_cells[idx] - sim_at_time
  return(sum(diff^2))
}

# Optimize growth rates:
best_rs <- optim(rep(0.6,8),cost, method = "L-BFGS-B",lower = 0)

# Plot simulated growth compared to data
slow_growth_plots <- simulate_cell_growth_variable(best_rs$par, cut_times, no_growth) %>% 
  ggplot(aes(day, cells/1e5)) + 
  geom_line(aes(color = "Updated cell growth simulations"), linewidth = 1, alpha = 0.7) + 
  geom_point(data=cell_growth, aes(day, live_cells)) + 
  theme(legend.position = "none") + 
  labs(color = "", x= "Time (days)", y = "Number of cells (.10^5)") + 
  scale_color_manual(values = c("darkred")) + 
  geom_label(data = data.frame(time = c(0,cut_times), dt = c(NA, round(log(2)/best_rs$par,1))) %>% 
               mutate(x = time - c(0,diff(time))/2) %>% filter(!is.na(dt)),
             aes(x, Inf, label = paste0("Td: ", dt), color = "Updated cell growth simulations"), 
             show.legend = F, vjust = 1.1) + 
  coord_cartesian(clip = "off")

ggsave(here(out_folder, "slow_cell_growth.png"), slow_growth_plots, width = 8, height = 3, bg = "white")


### Simulate KSHV dynamics with MLE parameters (start with 1e4 cells instead of 1e5 for speed)  ####
# Wrapper function to account for serial passaging and periods of no growth
sim_passage_wrapper_vary_b <- function(
    pRep, pSeg, bs, d, passage_times, initial_episomes, no_growth, i = 0
){
  
  if(i %% 10 == 0) cat("\n=== parameter set", i, "===\n")
  
  previous <- exponential_growth(
    pRep, pSeg, b = bs[1], d = d, nIts = 1e8,
    nRuns = 1, stop_size = 10^6.8, selection = F, max_epi = length(initial_episomes)-1,
    n_cells_start = sum(initial_episomes), stop_time = passage_times[1],
    initial_conditions = t(as.matrix(initial_episomes))
  )
  
  out <- previous
  
  for(i in 2:length(passage_times)){
    initial_conditions <- previous %>% group_by(run) %>% 
      filter(time == max(time), episomes != -1) %>% 
      # cut back to initial number of episomes
      mutate(num = frac*sum(initial_episomes)) %>% 
      ungroup %>% 
      distinct(run, episomes, num) %>% 
      pivot_wider(names_from = episomes, values_from = num) %>% 
      select(-run) %>% 
      as.matrix()
    
    if(length(no_growth) > 0 & no_growth[1] < passage_times[i]){
      stop_time <- no_growth[1] - passage_times[i-1]
      no_growth <- no_growth[-1]
    }else{
      stop_time <- passage_times[i]-passage_times[i-1]
      # start_times <- previous %>% group_by(run) %>% 
      #   filter(time == max(time)) %>% ungroup %>% distinct(run, time) %>% pull(time)
    }
    
    previous <- exponential_growth(
      pRep, pSeg, b = bs[i], d = d, nIts = 1e8, stop_size = 10^(6.8),
      nRuns = 1, selection = F, max_epi = length(initial_episomes)-1,
      n_cells_start = sum(initial_episomes), initial_conditions = initial_conditions,
      start_times = passage_times[i-1], stop_time = stop_time
    )
    
    out <- rbind(out, previous %>% filter(time != max(time)))
  }
  
  return(out)
}


# Simulate with MLE from SUM159 cells:
n_cell_start = 1e4
max_epi = 100
n_epi_start_fewer <- rnbinom(n_cell_start, mu = fit_nb_day0$estimate[2], size = fit_nb_day0$estimate[1])
initial_conditions_fewer <- rep(0,max_epi + 1)
count_epi <- table(n_epi_start_fewer)
initial_conditions_fewer[as.numeric(names(count_epi))+1] <- unname(count_epi)  

SUM159_MLE_vary_b <- sim_passage_wrapper_vary_b(pRep = 0.8, pSeg = 0.9, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer, no_growth)

# Simulate upper and lower bounds of confidence interval from fixed_KSHV images to capture uncertainty:
SUM159_confidence_intervals <- read_csv(here("results", "fixed_KSHV", "MLE_parameter_estimates.csv"))

Pr_CI <- as.numeric(str_extract_all(SUM159_confidence_intervals[[1,"marginal_CI"]], "[0-9.]+")[[1]])
Ps_CI <- as.numeric(str_extract_all(SUM159_confidence_intervals[[2,"marginal_CI"]], "[0-9.]+")[[1]])
mean_CI <- fit_confint[2,]

# Use the best estimate of size parameter as simulations indicate that results are not sensitive to 
# initial variation in episome counts per cell 
sample_initial_epi <- function(mu, size, n_cell_start = 1e4, max_epi = 150){
  n_epi_start_fewer <- rnbinom(n_cell_start, mu = mu, size = size)
  initial_conditions_fewer <- rep(0,max_epi + 1)
  count_epi <- table(n_epi_start_fewer)
  initial_conditions_fewer[as.numeric(names(count_epi))+1] <- unname(count_epi)
  return(initial_conditions_fewer)
}
big_epi <- sample_initial_epi(max(mean_CI), fit_nb_day0$estimate[1])
small_epi <- sample_initial_epi(min(mean_CI), fit_nb_day0$estimate[1])

brk219_CIs_upper <- expand_grid(pRep = Pr_CI[2], pSeg = Ps_CI) %>% mutate(i = 1:n()) %>% 
  pmap(sim_passage_wrapper_vary_b, bs = best_rs$par, d = 0, passage_times = cut_times, 
       initial_episomes = big_epi, no_growth = no_growth)

brk219_CIs_lower <- expand_grid(pRep = Pr_CI[1], pSeg = Ps_CI) %>% mutate(i = 1:n()) %>% 
  pmap(sim_passage_wrapper_vary_b, bs = best_rs$par, d = 0, passage_times = cut_times, 
       initial_episomes = small_epi, no_growth = no_growth)

brk219_CIs_df <- rbind(
  1:2 %>% map_df(function(i) cbind(expand_grid(pRep = Pr_CI[2], pSeg = Ps_CI)[i,], distinct(brk219_CIs_upper[[i]])) %>% mutate(trial = i)),
  1:2 %>% map_df(function(i) cbind(expand_grid(pRep = Pr_CI[1], pSeg = Ps_CI)[i,], distinct(brk219_CIs_lower[[i]])) %>% mutate(trial = i+2))  
)


save(brk219_CIs_df, file = here(out_folder, "CI_simulations.rds"))
# load(here(out_folder, "CI_simulations.rds"))


### Compare average number of episomes to rate of LANA dots decay ##############

## Calculate rate of LANA dot decay: 

# Convert day to generation for each time period: 
# match time period to growth rate (b)
temp <- cut_times
ts_full <- c()
for(i in 1:nrow(LANA_dots)){
  if(LANA_dots[i,"day"] < temp[1]){
    ts_full <- c(ts_full, temp[1])
  }else{
    temp <- temp[-1]
    ts_full <- c(ts_full, temp[1])
  }
}
bs <- data.frame(cut_time = cut_times, b = best_rs$par)
# calculate generataion = day*b
temp_df_full <- LANA_dots %>% mutate(cut_time = ts_full) %>% left_join(bs) %>% mutate(generation = day*b)

## Fit decay model with negative binomial and poisson distribution in MASS:
# Rate of decay = Pr-1 (since we converted from day to generation)

fit_nb_regression <- MASS::glm.nb(LANA_dots ~ generation, data = temp_df_full)
confint(fit_nb_regression)

fit_pois_regression <- MASS::glm(LANA_dots ~ generation, data = temp_df_full, family = "poisson")
confint(fit_pois_regression)

BIC(fit_nb_regression)
BIC(fit_pois_regression)


### Simulate with Pr fit to Brk.219 data: ######################################

# Pull the average number of LANA dots per cell at day 0 from the negative binomial regression
# (this is the fitted intercept, which needs to be converted from the log scale)
mu = exp(fit_nb_regression$coefficients[1])
initial_conditions_fewer_brk <- sample_initial_epi(
  mu, size = fit_nb_day0$estimate[1], n_cell_start = 1e4, max_epi = 125
  )

# Use pRep estimated from fitting exponential decay to  LANA dots in BRK.219
brk219_pRep <- round(1+fit_nb_regression$coefficients["generation"], 2)

BRK219_MLE_vary_b <- sim_passage_wrapper_vary_b(pRep = brk219_pRep, pSeg = 0.9, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer_brk, no_growth)


### Compare rate of LANA dot loss to predictions from model ###################
# Create common time grid
t_common <- seq(0, max(brk219_CIs_df$time), length.out = 200)

interpolated_sims <- brk219_CIs_df %>% 
  filter(pSeg == 1, episomes == -1) %>% 
  group_by(trial, pRep) %>% 
  reframe(
    frac = approx(x = time, y = frac, xout = t_common, rule = 2)$y,
    time = t_common
  ) %>% 
  filter(!is.na(frac))

ribbon_df <- interpolated_sims %>% 
  select(time, frac, pRep) %>% 
  mutate(pRep = ifelse(pRep == 0.9, "max", "min")) %>% 
  pivot_wider(names_from = pRep, values_from = frac) 


summary_plot <- brk219_CIs_df %>% 
  filter(episomes == -1) %>% 
  arrange(desc(pRep), desc(pSeg)) %>% 
  mutate(parameters = fct_inorder(interaction(pRep, pSeg)),
         pRep = ifelse(pRep == max(pRep), "Upper", "lower"),
         pSeg = ifelse(pSeg == max(pSeg), "Upper", "lower"))  %>% 
  ggplot(aes(time, frac)) +
  geom_ribbon(data = ribbon_df , 
              aes(time, ymin = min, ymax = max, fill = "Best estimates from\nfixed KSHV images and 95% PI"),
              alpha = 0.2, inherit.aes = F, color = NA) +
  geom_line(data = SUM159_MLE_vary_b %>% filter(episomes == -1), linewidth = 1.5, linetype = "solid", alpha = 0.8,
            aes(color = "Best estimates from\nfixed KSHV images and 95% PI")) + #  color = "gray39",
  geom_line(data = BRK219_MLE_vary_b %>% filter(episomes == -1), linewidth = 1.5, alpha = 0.8, aes(color = "Best fit to LANA dots")) + 
  geom_point(data = LANA_summary, aes(day, mean, color = "LANA dot mean and 95% CI"), size = 2.5) + 
  geom_errorbar(data = LANA_summary, aes(day, ymin = mean - 1.96*sem, ymax = mean + 1.96*sem, color = "LANA dot mean and 95% CI"), inherit.aes = F, width = 0.7, key_glyph = "vpath") + 
  theme(legend.position = "inside", legend.position.inside = c(1,1), legend.justification = c(1.1, 1.1), legend.background = element_rect(color = NA)) +
  labs(x = "Time (days)", y = "Average number of episomes or LANA dots per cell", color = NULL, fill = NULL) + 
  scale_color_manual(values = c("gray39", "black", "black"), aesthetics = c("color", "fill")) +
  guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))


ggsave(here(out_folder, "LANA_comparison_summary.pdf"), summary_plot, width = 5, height = 4)

