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
theme_set(theme_minimal())

source(here("scripts", "functions_simulations.R"))
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

out_folder <- here("results", "brk219")

### Inputs #####################################################################
LANA_dots <- read_csv(here("data", "derived", "brk219_full_LANA_dots.csv"))
GFP <- read_csv(here("data", "derived", "brk219_longitudinal.csv"))
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

### Define variation in initial episomes based on distribution of LANA dots at day 0 ####
day0_LANA <- LANA_dots %>% filter(day == 0) %>% pull(LANA_dots)

# Negative binomial fits better than poisson because variance > mean
fit_nb <- fitdist(day0_LANA, "nbinom")
fit_nb$estimate

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

### Compare cell growth data and percent GFP data ##############################

# Gray boxes denote periods where the two datasets don't align:
compare_cell_growth_GFP <- ggplot() +
  annotate(geom = "rect", xmin = 10, xmax = 18, ymin = -Inf, ymax = Inf, color = "lightgray", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 21, xmax = 25, ymin = -Inf, ymax = Inf, color = "lightgray", alpha = 0.1) +
  annotate(geom = "rect", xmin = 26, xmax = 29, ymin = -Inf, ymax = Inf, color = "lightgray", alpha = 0.1) +
  annotate(geom = "rect", xmin = 31, xmax = 35, ymin = -Inf, ymax = Inf, color = "lightgray", alpha = 0.1) +
  geom_point(data=cell_growth, aes(day, live_cells/20*100, color = "Live cells"), size = 3) + 
  geom_line(data=cell_growth, aes(day, live_cells/20*100, color = "Live cells"), linewidth = 1) + 
  geom_point(data=GFP, aes(day, percent_GFP_live, color = "GFP+ cells"), size = 3) + 
  geom_line(data=GFP, aes(day, percent_GFP_live, color = "GFP+ cells"), linewidth = 1) +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Percent GFP+ cells",
    limits = c(1,100),
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./100*20, name="Number of cells (.10^5)")
  ) + 
  
  theme(
    axis.title.y = element_text(color = "darkblue", size=13),
    axis.title.y.right = element_text(color = "gold3", size=13)
  ) +
  scale_color_manual(values = c("darkblue", "gold3"))  + 
  labs(x = "Time (days)", color = "")


ggsave(here(out_folder, "compare_GFP_cell_count.png"), compare_cell_growth_GFP, width = 8, height= 3.5)


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
n_epi_start_fewer <- rnbinom(n_cell_start, mu = fit_nb$estimate[2], size = fit_nb$estimate[1])
initial_conditions_fewer <- rep(0,max_epi + 1)
count_epi <- table(n_epi_start_fewer)
initial_conditions_fewer[as.numeric(names(count_epi))+1] <- unname(count_epi)  

brk219_MLE_vary_b <- sim_passage_wrapper_vary_b(pRep = 0.8, pSeg = 0.9, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer, no_growth)

# Calculate percent of cells with episomes over time
percent0_vary_growth <- brk219_MLE_vary_b %>% 
  filter(episomes == 0) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(linewidth = 1) + 
  geom_point(data = GFP, aes(day, percent_GFP_live)) + 
  labs(x = "Time (days)", y = "Percent of cells with episomes") + 
  ylim(c(0,100))

ggsave(here(out_folder, "percent0_vary_growth.png"), percent0_vary_growth, width = 4, height = 4)

# Compare observed to simulated value at day 30
brk219_MLE_vary_b %>% filter(episomes == 0) %>% filter(time == time[which.min(abs(time-30))])

### Propogate uncertainty in parameters ########################################

## Define grid of probable Pr and Ps values based on the posteriors from SUM159
fixed_KSHV_data <- load_data("fixed_KSHV_non_dividing_cells.xlsx", "fixed_KSHV_dividing_cells.xlsx")

daughter_cell_data <- fixed_KSHV_data$daughter_cell_data
mother_cell_data <- fixed_KSHV_data$mother_cell_data

fixed_KSHV_results <- run_pipeline(
  daughter_cell_data, mother_cell_data, 
  here("results", "fixed_KSHV"),
  n_prior = list("geom", 0.5), parallel = T, just_Pr = F,
  overwrite = F
)

grid_prob <- fixed_KSHV_results$MLE_grid$grid_search
set.seed(300)
idx <- sample(1:nrow(grid_prob), prob = grid_prob$probability, size = 100, replace = T)

param_grid <- grid_prob[idx, c("Pr", "Ps")] %>% 
  rename(pRep = Pr, pSeg = Ps) %>% mutate(i = 1:n())


# Sample combinations of Pr and Ps according to their likelihood in previous analysis
brk219_uncertainty_slower_growth <- param_grid %>% 
  pmap(sim_passage_wrapper_vary_b, b = best_rs$par, d = 0, no_growth = no_growth,
       initial_episomes = initial_conditions_fewer, passage_times = cut_times)

# save(brk219_uncertainty_slower_growth, param_grid, file =  here(out_folder, "brk219_sims_slower_growth.rds"))
load(here(out_folder, "brk219_sims_slower_growth.rds"))

brk219_uncertainty_slow_growth_df <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], brk219_uncertainty_slower_growth[[i]]) %>% mutate(trial = i))  %>% 
  distinct()

## Plot range of simulated trajectories

traj_by_Pr_slow <- brk219_uncertainty_slow_growth_df %>% 
  filter(episomes == 0) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = pRep, group = trial), alpha = 0.6) + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.74) + 
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == 0), linewidth = 1) + 
  scale_color_viridis_c("Pr") + 
  labs(x = "Time (days)", y = "Percent of cells with episomes")

traj_by_Ps_slow <- brk219_uncertainty_slow_growth_df %>% 
  filter(episomes == 0) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = pSeg, group = trial), alpha = 0.6) + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.74) + 
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == 0), linewidth = 1) + 
  scale_color_viridis_c("Ps") + 
  labs(x = "Time (days)", y = "Percent of cells with episomes")

## Function to Calculate 95% prediction interval for either the average number of episomes or the percent of cells with episomes across simulations:
calculate_percentiles <- function(df, outcome){
  
  if(outcome == "mean" | outcome == "avg"){
    df <- df %>% filter(episomes == -1)
  }else if(outcome == "perc_epi"){
    df <- df %>% filter(episomes == 0) %>% mutate(frac = 100*(1-frac))
  }
  
  # Create common time grid
  t_common <- seq(0, max(df$time), length.out = 200)
  
  interpolated_sims <- df %>% 
    group_by(trial) %>% 
    reframe(
      perc = approx(x = time, y = frac, xout = t_common, rule = 2)$y,
      time = t_common
    ) %>% 
    filter(!is.na(perc))
  
  # Calculate mean and 95% CI at each time point
  percentiles <- interpolated_sims %>%
    group_by(time) %>%
    summarise(
      mean = mean(perc),
      median = median(perc),
      lower = quantile(perc, 0.025),
      upper = quantile(perc, 0.975),
      .groups = 'drop'
    )
  
  return(percentiles)
}


percentiles_0_slow_growth <- calculate_percentiles(brk219_uncertainty_slow_growth_df, "perc_epi")

traj_percentiles_slow <- ggplot(percentiles_0_slow_growth, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.3) +
  # geom_line(aes(y = mean, color = "mean"),  linewidth = 1) +
  # geom_line(aes(y = median, color = "median"),  linewidth = 1) +
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == 0),
            aes(time, 100*(1-frac), color = "MLE"), linewidth = 1 ) + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.75) +
  labs(x = "Time (days)", y = "Percent of cells with episomes") + 
  # scale_color_manual(values = c("blue", "darkgreen", "black"))
  scale_color_manual(values = c("black")) + 
  scale_fill_manual(values = c("blue")) + 
  theme(legend.title = element_blank())

summary_percent_slow <- traj_by_Pr_slow + traj_percentiles_slow & plot_annotation(tag_levels = "a") &
  theme(legend.position = "inside", legend.position.inside = c(0,0), legend.justification = c(-0.1,-0.1),
        legend.background = element_rect(color = NA))

ggsave(here(out_folder, "percent_epi_remaining_slow_growth.png"), summary_percent_slow, width = 8, height = 4, bg = "white")



### Compare simulated decay in average episome per cell to decay in average LANA dots and mean GFP ###############################

average_decay_by_Pr_slow_growth <- brk219_uncertainty_slow_growth_df %>% 
  filter(episomes == -1) %>% 
  ggplot(aes(time, frac/max(frac)*100)) + 
  geom_line(aes(color = pRep, group = trial), alpha = 0.6) + 
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == -1), linewidth  = 1) +
  geom_point(data = GFP, aes(day, mean_GFP_live/max(mean_GFP_live, na.rm = T)*100, shape = "mean GFP"), inherit.aes = F, size = 2.5) +
  geom_point(data = LANA_summary, aes(day, mean/max(mean, na.rm = T)*100, shape = "mean LANA dots"), inherit.aes = F, size = 2.5) +
  geom_errorbar(data = LANA_summary %>% mutate(min = mean-1.96*sem, max = mean+1.96*sem),
                aes(day, ymin = min/max(mean, na.rm = T)*100, ymax = max/max(mean, na.rm = T)*100), inherit.aes = F) +
  labs(x = "Time (days)", y = "Percent of signal remaining", color = "Replication\nEfficiency", shape = "") + 
  scale_color_viridis_c(limits = c(0.5,0.95)) +
  # scale_y_log10(breaks = c(1,10,100)) + 
  coord_cartesian(ylim = c(0,100))


percentiles_mean_slow_growth <- calculate_percentiles(brk219_uncertainty_slow_growth_df, "mean")

mean_percentiles_slow_growth <- percentiles_mean_slow_growth %>% 
  mutate(across(-time, ~.x/max(.x)*100)) %>% 
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% PI"), alpha = 0.3) +
  # geom_line(aes(y = mean, color = "mean"),  linewidth = 1) +
  # geom_line(aes(y = median, color = "median"),  linewidth = 1) +
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == -1),
            aes(time, frac/max(frac)*100, color = "MLE"), linewidth = 1 ) + 
  geom_point(data = GFP, aes(day, mean_GFP_live/max(mean_GFP_live, na.rm = T)*100, shape = "mean GFP"), inherit.aes = F, size = 2.5) + 
  geom_point(data = LANA_summary, aes(day, mean/max(mean, na.rm = T)*100, shape = "mean LANA dots"), inherit.aes = F, size = 2.5) + 
  geom_errorbar(data = LANA_summary %>% mutate(min = mean-1.96*sem, max = mean+1.96*sem),
                aes(day, ymin = min/max(mean, na.rm = T)*100, ymax = max/max(mean, na.rm = T)*100), inherit.aes = F) +
  labs(x = "Time (days)", y = "Percent of signal remaining", color = "", shape = "", fill = "") + 
  # scale_color_manual(values = c("blue", "darkgreen", "black")) + 
  scale_color_manual(values = c("black")) + 
  scale_fill_manual(values = c("blue")) + 
  # scale_y_log10(breaks = c(1,10,100)) + 
  coord_cartesian(ylim = c(0,100))

average_decay_plot <- average_decay_by_Pr_slow_growth  + mean_percentiles_slow_growth &
  plot_annotation(tag_levels = "a",title = "Average episomes per cell across all live cells") & 
  theme(plot.title = element_text(hjust = 0.5, face = 2))

ggsave(here("results", "brk219", "average_decay_slow_growth_LANA_plus_GFP.png"), average_decay_plot, width = 12, height = 5, bg = "white")

### Compare average number of episomes to LANA dots alone #######################

average_decay_by_Pr_slow_growth2 <- brk219_uncertainty_slow_growth_df %>% 
  filter(episomes == -1) %>% 
  ggplot(aes(time, frac)) + 
  geom_line(aes(color = pRep, group = trial), alpha = 0.6) + 
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == -1), linewidth  = 1) +
  geom_point(data = LANA_summary, aes(day, mean, shape = "mean LANA dots"), inherit.aes = F, size = 2.5) +
  geom_errorbar(data = LANA_summary %>% mutate(min = mean-1.96*sem, max = mean+1.96*sem),
                aes(day, ymin = min, ymax = max), inherit.aes = F) +
  labs(x = "Time (days)", y = "Average episomes per cell", color = "Replication\nEfficiency", shape = "") + 
  scale_color_viridis_c(limits = c(0.5,0.95)) 


mean_percentiles_slow_growth2 <- percentiles_mean_slow_growth %>% 
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% PI"), alpha = 0.3) +
  geom_line(data = brk219_MLE_vary_b %>% filter(episomes == -1),
            aes(time, frac, color = "MLE"), linewidth = 1 ) + 
  geom_point(data = LANA_summary, aes(day, mean, shape = "mean LANA dots"), inherit.aes = F, size = 2.5) +
  geom_errorbar(data = LANA_summary %>% mutate(min = mean-1.96*sem, max = mean+1.96*sem),
                aes(day, ymin = min, ymax = max), inherit.aes = F) +
  labs(x = "Time (days)", y = "Average episomes per cell", color = "", shape = "", fill = "") + 
  # scale_color_manual(values = c("blue", "darkgreen", "black")) + 
  scale_color_manual(values = c("black")) + 
  scale_fill_manual(values = c("blue"))  


average_decay_plot2 <- average_decay_by_Pr_slow_growth2  + 
  guides(
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)
  ) +
  mean_percentiles_slow_growth2 +
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  ) &
  plot_annotation(tag_levels = "a") & 
  theme(plot.title = element_text(hjust = 0.5, face = 2)) &
  theme(legend.position = "inside", legend.position.inside = c(1,1), legend.justification = c(1,1), 
        legend.background = element_rect(fill = "white", color= "white"), legend.spacing.y = unit(0, "npc")) 

ggsave(here("results", "brk219", "average_decay_slow_growth.png"), average_decay_plot2, width = 9, height = 5, bg = "white")

ggsave(here("results", "brk219", "average_decay_PI.png"), mean_percentiles_slow_growth2 + theme(legend.position = "none"), 
       width = 6, height = 4, bg = "white")


### Compare distribution of episomes per cell to distribution of LANA dots per cell over time  ####

## Sample the distribution of episomes per cell at times when experimental data were recorded according to simulations: 
sample_episomes <- function(simualtion){
  times <- approx(unique(simualtion %>% pull(time)), xout = unique(LANA_dots$day))$y
  
  # Create data frame with episomes:
  temp <- simualtion %>% 
    filter(round(time,1) %in% unique(LANA_dots$day), episomes != -1) %>% 
    group_by(day = round(time)) %>% 
    filter(time == min(time))
  
  epi_for_dist3 <- NULL
  for(Day in unique(temp$day)){
    t <- temp %>% filter(day == Day)
    total <- unique(t$total)
    episomes <- c()
    for(i in 1:nrow(t)){
      episomes <- c(episomes, rep(t$episomes[i], t$frac[i]*total))
    }
    epi_for_dist3 <- rbind(epi_for_dist3, data.frame(day = Day, episomes = episomes))
  }
  
  return(epi_for_dist3)
}

# Maximum likelihood estimates
epi_for_dist <- sample_episomes(brk219_MLE_vary_b)

## Add in simulations with Pr = 87%:
brk219_Pr87 <- sim_passage_wrapper_vary_b(pRep = 0.87, pSeg = 0.9, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer, no_growth)
epi_for_dist2 <- sample_episomes(brk219_Pr87)

## Add in simulations with Pr = 87% and random segregation
brk219_Pr87_randomSeg <- sim_passage_wrapper_vary_b(pRep = 0.87, pSeg = 0.5, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer, no_growth)
epi_for_dist3 <- sample_episomes(brk219_Pr87_randomSeg)

## Add in simulations with Pr = 90% and random segregation
brk219_Pr8_randomSeg <- sim_passage_wrapper_vary_b(pRep = 0.8, pSeg = 0.5, bs = best_rs$par, d = 0, cut_times, initial_conditions_fewer, no_growth)
epi_for_dist4 <- sample_episomes(brk219_Pr8_randomSeg)

## Plot results: s
simulation_data_frame <- rbind(
  epi_for_dist %>% mutate(sim = "Simulations with 80% Replication Efficiency and 90% Segregation Efficiency"),
  epi_for_dist4 %>% mutate(sim = "Simulations with 80% Replication Efficiency and 50% Segregation Efficiency"),
  epi_for_dist2 %>% mutate(sim = "Simulations with 87% Replication Efficiency and 90% Segregation Efficiency"),
  epi_for_dist3 %>% mutate(sim = "Simulations with 87% Replication Efficiency and Random Segregation")
) %>% 
  mutate(sim = fct_inorder(sim))

# calculate ks test for each day: 
ks_summary <- NULL
for(Day in seq(0,20,by=2)){
  sim1 <- epi_for_dist %>% filter(day == Day) %>% pull(episomes)
  sim2 <- epi_for_dist2 %>% filter(day == Day) %>% pull(episomes)
  sim3 <- epi_for_dist3 %>% filter(day == Day) %>% pull(episomes)
  sim4 <- epi_for_dist4 %>% filter(day == Day) %>% pull(episomes)
  data <- LANA_dots %>% filter(day == Day) %>% pull(LANA_dots)
  ks_summary<- rbind(ks_summary, data.frame(day = Day, V1 = ks.test(sim1, data)$p.value, 
                                            V2 = ks.test(sim4, data)$p.value,
                                            V3= ks.test(sim2, data)$p.value,
                     V4 = ks.test(sim3, data)$p.value))
}

ks_summary <- ks_summary %>% 
  pivot_longer(!day, names_to = "sim", values_to = "p.value")  %>% 
  mutate(sim = fct_inorder(sim), sim = factor(sim, labels = levels(simulation_data_frame$sim)))

plot_distribution_comparison <- function(simulation){
  p <- simulation_data_frame %>%
     filter(as.numeric(sim) == simulation) %>% 
     ggplot(aes(factor(day), episomes)) + 
     geom_violin(alpha = 0.3, aes(fill = sim,  color = sim), adjust = 2, scale = "width") + 
     stat_dots(
       data = LANA_dots %>% rename(episomes = LANA_dots),
       layout = "swarm",      # Replicates beeswarm spacing
       side = "both",         # Distributes points evenly on both sides of the axis
       dotsize = 1,           # Adjust dot size
       binwidth = NA,         # Automatically calculates best bin width,
       color = "black",
       fill = "black",
       alpha = 0.5
     ) +
     geom_text(data = ks_summary %>% filter(p.value < 0.01, as.numeric(sim) == simulation), aes(factor(day), y = Inf, label = "*"), vjust = 1.1, size = 7) + 
     geom_line(data = . %>% group_by(day, sim) %>% summarise(mean = mean(episomes)), aes(x =as.numeric(as.factor(day)), y = mean, color = sim), linewidth = 1) + 
     geom_line(data = LANA_summary, aes(x =as.numeric(as.factor(day)), y = mean), linewidth = 1, alpha = 0.7) + 
     labs(x = "Day", y = "Number of episomes or\nLANA dots per cell") + 
     theme(legend.position = "none") +
     facet_wrap(~sim, ncol = 1) + 
     scale_color_viridis_d(end = 0.8, drop = F) + 
     scale_fill_viridis_d(end = 0.8, drop = F) 
  
  return(p)
  
}

compare_dists_over_time3 <- plot_distribution_comparison(1) + plot_distribution_comparison(3) + plot_distribution_comparison(4) &
  plot_layout(ncol = 1) & theme(strip.text = element_text(size = 12))

ggsave(here(out_folder, "distributions_of_episomes_over_time.png"), compare_dists_over_time3, width = 9, height = 10)

### Find values of replication and segregation efficiency that might better explain the observed data #####

## Simulate alternative pRep and pSeg: 
parameter_sensitivity_grid <- expand_grid(
  pRep = c(0.8, 0.85, 0.9, 0.92, 0.94, 0.95, 0.96, 1),
  pSeg = c(0.1, 0.4, 0.5, 0.7, 0.9, 1)
)

parameter_sensitivity_sims <- parameter_sensitivity_grid %>% 
  mutate(i = 1:nrow(.)) %>% 
  pmap(sim_passage_wrapper_vary_b,
       bs = best_rs$par, d = d_0, passage_times = cut_times, 
       initial_episomes = initial_conditions_fewer, no_growth = no_growth)

brk219_slow_growth_sensitivity <- 1:nrow(parameter_sensitivity_grid) %>% 
  map_df(function(i) cbind(parameter_sensitivity_grid[i,], parameter_sensitivity_sims[[i]]) %>% mutate(trial = i))  %>% 
  distinct() %>% 
  arrange(pRep, pSeg) %>% 
  mutate(param = fct_inorder(paste0("Pr: ",pRep*100, "%, Ps: ", pSeg*100, "%"))) 

### Plot the number of cells without episomes for each simulation and compare to data:

# 1) Is there a segregation efficiency for 80% replication efficiency that better explains the data?
pSeg_sensitivity <- brk219_slow_growth_sensitivity %>% 
  filter(episomes == 0, pRep == 0.8) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = param), linewidth = 1, alpha = 0.7)  + 
  geom_line(data = ~filter(.x, pSeg == 0.9), lty = "dashed")  + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.75) + 
  labs(x = "Time (days)", y = "Percent of cells with episomes", color = "Model Parameters") +
  scale_color_manual(values = scales::brewer_pal(palette = "Blues")(7)[2:7]) + 
  guides(color = guide_legend(reverse = TRUE))

# 2) Is there a replication efficiency for 90% segregation efficiency that better explains the data?
pRep_sensitivity <- brk219_slow_growth_sensitivity %>% 
  filter(episomes == 0, pSeg == 0.9) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = param), linewidth = 1, alpha = 0.7)  + 
  geom_line(data = ~filter(.x, pRep == 0.8), lty = "dashed")  + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.75) + 
  labs(x = "Time (days)", y = "Percent of cells with episomes", color = "Model Parameters") +
  scale_color_manual(values = scales::brewer_pal(palette = "RdPu")(8)[2:8]) + 
  guides(color = guide_legend(reverse = TRUE))

individual_sensitivity <- pSeg_sensitivity + ggtitle("Replication Efficiency: 80%") + 
  pRep_sensitivity + ggtitle("Segregation Efficiency: 90%") & 
  plot_annotation(tag_levels = "a") & ylim(c(0,100)) & theme(plot.title = element_text(hjust = 0.5))
ggsave(here(out_folder, "individual_parameter_sensitivies.png"), individual_sensitivity, width = 12, height = 4, bg = "white")

# 3) Is there a different combination that explains the data best?
pRep_pSeg_sensitivity <- brk219_slow_growth_sensitivity %>% 
  filter(episomes == 0, pRep != 0.8) %>%
  arrange(pRep) %>% 
  mutate(pRep = fct_inorder(paste0("Replication Efficiency: ", pRep*100, "%"))) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = factor(pSeg*100)), linewidth = 1, alpha = 0.7)  + 
  geom_line(data = filter(brk219_slow_growth_sensitivity, pRep == 0.8, pSeg == 0.9, episomes == 0) %>% select(-pRep, -pSeg), lty = "dashed")  + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.75) + 
  ggrepel::geom_text_repel(data = ~group_by(.x, pRep, pSeg) %>% filter(time == max(time)), 
                           aes(label = pSeg*100, color = factor(pSeg*100)), xlim = c(43,48), 
                           min.segment.length = 0, show.legend = F, size = 3) +
  labs(x = "Time (days)", y = "Percent of cells with episomes", color = "Segregation\nEfficiency (%)") +
  scale_color_manual(values = scales::brewer_pal(palette = "Blues")(7)[2:7]) +
  facet_wrap(~pRep) + ylim(c(0,100)) + theme(strip.text = element_text(size = 12)) + 
  guides(color = guide_legend(reverse = TRUE))

ggsave(here(out_folder, "full_parameter_sensitivities.png"), pRep_pSeg_sensitivity, width = 12, height = 6, bg = "white")

pRep_pSeg_sensitivity_subset <- rbind(brk219_slow_growth_sensitivity, brk219_slow_growth_sensitivity2) %>% 
  filter(episomes == 0, pRep %in% c(0.8, 0.85, 0.9, 0.95)) %>%
  arrange(pRep) %>% 
  mutate(pRep = fct_inorder(paste0("Replication Efficiency: ", pRep*100, "%"))) %>% 
  ggplot(aes(time, 100*(1-frac))) + 
  geom_line(aes(color = factor(pSeg*100)), linewidth = 1, alpha = 0.7)  + 
  geom_line(data = filter(brk219_slow_growth_sensitivity, pRep == 0.8, pSeg == 0.9, episomes == 0) %>% select(-pRep, -pSeg), lty = "dashed")  + 
  geom_point(data = GFP, aes(day, percent_GFP_live), size = 2, alpha = 0.75) + 
  ggrepel::geom_text_repel(data = ~group_by(.x, pRep, pSeg) %>% filter(time == max(time)), 
                           aes(label = pSeg*100, color = factor(pSeg*100)), xlim = c(43,48), 
                           min.segment.length = 0, show.legend = F, size = 3) +
  labs(x = "Time (days)", y = "Percent of cells with episomes", color = "Segregation\nEfficiency (%)") +
  scale_color_manual(values = scales::brewer_pal(palette = "Blues")(7)[2:7]) +
  facet_wrap(~pRep) + ylim(c(0,100)) + theme(strip.text = element_text(size = 12)) + 
  guides(color = guide_legend(reverse = TRUE))

ggsave(here(out_folder, "subset_parameter_sensitivities.png"), pRep_pSeg_sensitivity_subset, width = 8, height = 6, bg = "white")


### Plot the average number of episomes per cell for each simulation and compare to data:

pRep_pSeg_sensitivity_average <- brk219_slow_growth_sensitivity %>% 
  filter(episomes == -1, pRep != 0.8) %>%
  arrange(pRep) %>% 
  mutate(pRep = fct_inorder(paste0("Replication Efficiency: ", pRep*100, "%"))) %>% 
  ggplot(aes(time, frac)) + 
  geom_line(aes(color = factor(pSeg*100)), linewidth = 1, alpha = 0.7)  + 
  geom_line(data = filter(brk219_slow_growth_sensitivity, pRep == 0.8, pSeg == 0.9, episomes == -1) %>% select(-pRep, -pSeg), lty = "dashed")  + 
  geom_point(data = LANA_summary, aes(day, mean, shape = "mean LANA dots"), inherit.aes = F, size = 2.5) +
  geom_errorbar(data = LANA_summary %>% mutate(min = mean-sd, max = mean+sd),
                aes(day, ymin = min, ymax = max), inherit.aes = F) +
  labs(x = "Time (days)", y = "Average number of episomes", color = "Segregation\nEfficiency (%)", shape = "") +
  scale_color_manual(values = scales::brewer_pal(palette = "Blues")(7)[2:7]) +
  facet_wrap(~pRep) + theme(strip.text = element_text(size = 12)) + 
  guides(color = guide_legend(reverse = TRUE))

ggsave(here(out_folder, "full_parameter_sensitivities_average_epi.png"), pRep_pSeg_sensitivity_average, width = 12, height = 6, bg = "white")

