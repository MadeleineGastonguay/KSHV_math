# This script benchmarks the likelihood in synthetic data

### Setup ###################################################################### 

# load libraries 
library(tidyverse)
library(here)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

### Script inputs ##############################################################

set.seed(400)

# test a variety of Pr and Ps values and numbers of daughter cell pairs:
values <- c(0.1, 0.3, 0.5, 0.8, 0.9)
n <- c(10, 40, 100)

sim_parameters <- expand_grid(Pr = values, Ps = values, n = n)

### Functions to simulate cell division ########################################

# Function to simulate episome data
sim_one_cell <- function(X0, Pr, Ps, id){
  # How many episomes replicate?
  r <- sum(rbinom(X0, 1, Pr))
  # How many replicated episomes segregate?
  s <- sum(rbinom(r, 1, Ps))
  # Assign k replicated pairs to cell 1
  k <- sum(rbinom(r-s, 1, 0.5))
  # Assign j singletons to cell 1
  j <- sum(rbinom(X0-r, 1, 0.5))
  
  # Count episomes in cell 1 and cell 2
  X1 <- s + 2*k + j
  X2 <- X0 + r - (s + 2*k + j)
  
  data.frame(r, s, k, j, X1 = max(X1, X2), X2 = min(X1, X2))
}

# Function to simulate multiple cells 
# X0s is a vector of length n_cells with the initial number of episomes. Alternatively, it can be a single value if all cells have the same to start.
# Pr and Ps are the probability of replication and division, respectively
# n_cells is the number of cells to simulate
simulate_multiple_cells <- function(X0s, Pr, Ps, n_cells){
  if(length(X0s) != n_cells & length(X0s) != 1) stop("Incorect dimension for X0s")
  
  data.frame(X0 = X0s, Pr, Ps, id = 1:n_cells) %>% 
    bind_cols(pmap_df(., sim_one_cell))
}

### Simulate data ##############################################################

sim_data <- sim_parameters %>% 
  pmap(function(Pr, Ps, n){
    # Sample the number of episomes in the mother cell from a poisson distribution
    X0s <- sample(1:100, n, replace = T, prob = dpois(1:100, 3))
    # Simulate cell division for each mother cell 
    simulate_multiple_cells(Pr = Pr, Ps = Ps, n_cells = n, X0s = X0s )
  } ) %>% 
  list_rbind(names_to = "simulation")


### Apply MLE search to synthetic data #########################################

## Assuming X0 is known (ideal scenario):
MLE_with_known_X0 <- sim_data %>% 
  group_by(simulation) %>% 
  nest() %>% 
  pull(data) %>% 
  map(run_grid_search, viz = F, known_X0 = T, marginal_CI = T)

MLE_with_known_X0_estimates <- MLE_with_known_X0 %>% lapply("[[", "estimates") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% 
  mutate(rep = 1)

MLE_with_known_X0_grid <- MLE_with_known_X0 %>% lapply("[[", "grid_search") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)


## Assuming X0 is unknown (limitation of our experimental data) :
MLE_with_unknown_X0 <- sim_data %>%
  group_by(simulation) %>%
  nest() %>%
  pull(data) %>%
  map(run_grid_search, viz = F, known_X0 = F, lambda = 3, marginal_CI = T)

MLE_with_unknown_X0_estimates <- MLE_with_unknown_X0 %>% lapply("[[", "estimates") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)

MLE_with_unknown_X0_grid <- MLE_with_unknown_X0 %>% lapply("[[", "grid_search") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)

data <- list(sim_data)

## Repeat 9 more times:
for(i in 2:10){
  sim_data.2 <- sim_parameters %>% 
    pmap(function(Pr, Ps, n){
      X0s <- sample(1:100, n, replace = T, prob = dpois(1:100, 3))
      simulate_multiple_cells(Pr = Pr, Ps = Ps, n_cells = n, X0s = X0s )
    } ) %>% 
    list_rbind(names_to = "simulation")
  
  data[[i]] <- sim_data.2
  
  MLE_with_unknown_X0.2 <- sim_data.2 %>% 
    group_by(simulation) %>%
    nest() %>%
    pull(data) %>%
    map(run_grid_search, viz = F, known_X0 = F, lambda = 3)
  
  MLE_with_unknown_X0_estimates <- rbind(MLE_with_unknown_X0_estimates, 
                                         MLE_with_unknown_X0.2 %>% lapply("[[", "estimates") %>% 
                                           list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_unknown_X0_grid <- rbind(MLE_with_unknown_X0_grid, 
                                    MLE_with_unknown_X0.2 %>% lapply("[[", "grid_search") %>% 
                                      list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_known_X0.2 <- sim_data.2 %>% 
    group_by(simulation) %>%
    nest() %>%
    pull(data) %>%
    map(run_grid_search, viz = F, known_X0 = T)
  
  MLE_with_known_X0_estimates <- rbind(MLE_with_known_X0_estimates, 
                                       MLE_with_known_X0.2 %>% lapply("[[", "estimates") %>% 
                                         list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_known_X0_grid <- rbind(MLE_with_known_X0_grid, 
                                  MLE_with_known_X0.2 %>% lapply("[[", "grid_search") %>% 
                                    list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
}


## Gather and save results from multiple simulations
MLE_with_unknown_X0_estimates <- MLE_with_unknown_X0_estimates %>%
  left_join(sim_parameters %>% mutate(simulation = 1:nrow(.)),  by = "simulation")

MLE_with_known_X0_estimates <- MLE_with_known_X0_estimates %>%
  left_join(sim_parameters %>% mutate(simulation = 1:nrow(.)),  by = "simulation")


# save(MLE_with_known_X0_estimates,  data,  MLE_with_known_X0_grid,
#      file = here("results", "benchmarking", "MLE_results_knownX0.RData"))

# save(MLE_with_unknown_X0_estimates, data, MLE_with_unknown_X0_grid,
#      file = here("results", "benchmarking", "MLE_results_unknownX0.RData"))


### Plot results ###############################################################

## Estimate versus true for Pr when X0 is unknown
Pr_estimates_unknown <- MLE_with_unknown_X0_estimates %>% 
  mutate(n = factor(paste(n, "Daughter cell pairs"), levels = paste(c(10, 40, 100), "Daughter cell pairs")),
         Ps = paste0("Segregation Efficiency: ", Ps*100, "%")) %>% 
  ggplot(aes(Pr*100, MLE_Pr*100, color = min_Pr <= Pr & max_Pr >= Pr)) + 
  geom_abline(alpha = 0.5) +  
  # Plot confidence intervals for each simulated data set
  geom_linerange(aes(ymin = min_Pr*100, ymax=  max_Pr*100), 
                 position = position_dodge2(width = 7), alpha = 0.7) + 
  # Plot ML estimate for each simulated data set
  geom_point(position = position_dodge2(width = 7)) + 
  # Plot mean of ML estimates across simulated data sets in red
  geom_point(data = . %>% group_by(Pr, Ps, n) %>% summarise(MLE_Pr = mean(MLE_Pr)), 
             color = "#DF536B", fill = "#DF536B", shape = 24, alpha = 0.7, size= 2.5) + 
  # Format plot
  facet_grid(n ~ Ps) + 
  theme(legend.position = "none") + 
  # Color confidence intervals that exclude true estimate in black
  scale_color_manual(values = c("black", "gray")) + 
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  labs(x = "True Replication Efficiency (%)", y = "Maximum Likleihood Estimate of Replication Efficiency (%)")

ggsave(here("results", "benchmarking", "Pr_estimates.png"), Pr_estimates_unknown, width = 12, height = 8)


## Estimate versus true for Ps when X0 is unknown
Ps_estimates_unknown <- MLE_with_unknown_X0_estimates %>% 
  mutate(n = factor(paste(n, "Daughter cell pairs"), levels = paste(c(10, 40, 100), "Daughter cell pairs")),
         Pr = paste0("Replication Efficiency: ", Pr*100, "%")) %>% 
  ggplot(aes(Ps*100, MLE_Ps*100, color = min_Ps <= Ps & max_Ps >= Ps)) + 
  geom_abline(alpha = 0.5) + 
  geom_linerange(aes(ymin = min_Ps*100, ymax=  max_Ps*100), 
                 position = position_dodge2(width = 7), alpha = 0.7) + 
  geom_point(position = position_dodge2(width = 7)) + 
  geom_point(data = . %>% group_by(Ps, Pr, n) %>% summarise(MLE_Ps = mean(MLE_Ps)), 
             color = "#DF536B", fill = "#DF536B", shape = 24, alpha = 0.7, size= 2.5) + 
  facet_grid(n ~ Pr) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("black", "gray")) + 
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  labs(x = "True Segregation Efficiency (%)", y = "Maximum Likleihood Estimate of Segregation Efficiency (%)")

ggsave(here("results", "benchmarking", "Ps_estimates.png"), Ps_estimates_unknown, width = 12, height = 8)


## Estimate versus true for Pr when X0 is known
Pr_estimates_known <- MLE_with_known_X0_estimates %>% 
  mutate(n = factor(paste(n, "Daughter cell pairs"), levels = paste(c(10, 40, 100), "Daughter cell pairs")),
         Ps = paste0("Segregation Efficiency: ", Ps*100, "%")) %>% 
  ggplot(aes(Pr*100, MLE_Pr*100, color = min_Pr <= Pr & max_Pr >= Pr)) + 
  geom_abline(alpha = 0.5) + 
  geom_linerange(aes(ymin = min_Pr*100, ymax=  max_Pr*100), 
                 position = position_dodge2(width = 7), alpha = 0.7) + 
  geom_point(position = position_dodge2(width = 7)) + 
  geom_point(data = . %>% group_by(Pr, Ps, n) %>% summarise(MLE_Pr = mean(MLE_Pr)), 
             color = "#DF536B", fill = "#DF536B", shape = 24, alpha = 0.7, size= 2.5) + 
  facet_grid(n ~ Ps) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("black", "gray")) + 
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  labs(x = "True Replication Efficiency (%)", y = "Maximum Likleihood Estimate of Replication Efficiency (%)")

ggsave(here("results", "benchmarking", "Pr_estimates_knownX0.png"), Pr_estimates_known, width = 12, height = 8)


## Estimate versus true for Ps when X0 is known
Ps_estimates_known <- MLE_with_known_X0_estimates %>% 
  mutate(n = factor(paste(n, "Daughter cell pairs"), levels = paste(c(10, 40, 100), "Daughter cell pairs")),
         Pr = paste0("Replication Efficiency: ", Pr*100, "%")) %>% 
  ggplot(aes(Ps*100, MLE_Ps*100, color = min_Ps <= Ps & max_Ps >= Ps)) + 
  geom_abline(alpha = 0.5) + 
  geom_linerange(aes(ymin = min_Ps*100, ymax=  max_Ps*100), 
                 position = position_dodge2(width = 7), alpha = 0.7) + 
  geom_point(position = position_dodge2(width = 7)) + 
  geom_point(data = . %>% group_by(Ps, Pr, n) %>% summarise(MLE_Ps = mean(MLE_Ps)), 
             color = "#DF536B", fill = "#DF536B", shape = 24, alpha = 0.7, size= 2.5) + 
  facet_grid(n ~ Pr) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("black", "gray")) + 
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,100)) + 
  labs(x = "True Segregation Efficiency (%)", y = "Maximum Likleihood Estimate of Segregation Efficiency (%)")

ggsave(here("results", "benchmarking", "Ps_estimates_knownX0.png"), Ps_estimates_known, width = 12, height = 8)
