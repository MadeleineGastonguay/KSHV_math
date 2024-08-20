## This script assess the performance of the MCMC algorithm on simulated data

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)
# library(furrr)

# future::plan(multisession, workers = 4)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

#####
# Define new gibbs sampler 
#####
run_gibbs_new <- function(tau0, mu0, I, n_iterations, ns = NA){
  clusters <- names(I)
  q <- length(I)
  
  # Initial guesses
  tau <- rep(NA, n_iterations)
  tau[1] <- tau0
  mu <- rep( NA, n_iterations)
  mu[1] <- mu0
  
  # initialize matrix of ns
  n <- matrix(NA, ncol = q, nrow = n_iterations)
  if(any(is.na(ns))){
    n[1,] <- round(I/mu[1])
    n[1,][n[1,] == 0] <- 1  
    colnames(n) <- clusters
  }else{
    n[1,] <- ns
    colnames(n) <- names(ns)
  }
  # make sure n is in correct order
  n <- n[,clusters]
  
  nks <- seq(1,100) # define possible values of number of episomes per cluster
  for(j in 2:n_iterations){
    for(k in 1:q){
      # define the probability of each nk given observed data and other parameters
      nk_like <- nks %>% sapply(log_likelihood_n, mu = mu[j-1], sigma2 = 1/tau[j-1], I = I[k]) # log likelihood
      nk_probs <- exp(nk_like - max(nk_like))/sum(exp(nk_like - max(nk_like)), na.rm = T) #probabilities
      # sample nk
      n[j, k] <- sample(nks, 1, prob = nk_probs, replace =  T)
    }
    # sample mu:
    mu[j] <- rnorm(1, sum(I)/sum(n[j,]), sqrt(1/(tau[j-1]*sum(n[j,]))))
    # sample tau:
    tau[j] <- rgamma(1, q/2 + 1, 0.5*(sum(I^2/n[j,]) - 2*mu[j]*sum(I) + mu[j]^2*sum(n[j,])))
  }
  
  return(cbind(iteration = 1:n_iterations, mu, tau, as.data.frame(n)))
}

# function to run Gibbs sampling
run_gibbs_old <- function(tau0, mu0, I, n_iterations, ns = NA){
  clusters <- names(I)
  q <- length(I)
  
  # Initial guesses
  tau <- rep(NA, n_iterations)
  tau[1] <- tau0
  mu <- rep( NA, n_iterations)
  mu[1] <- mu0
  
  # initialize matrix of ns
  n <- matrix(NA, ncol = q, nrow = n_iterations)
  if(any(is.na(ns))){
    n[1,] <- round(I/mu[1])
    n[1,][n[1,] == 0] <- 1  
    colnames(n) <- clusters
  }else{
    n[1,] <- ns
    colnames(n) <- names(ns)
  }
  # make sure n is in correct order
  n <- n[,clusters]
  
  nks <- seq(1,100) # define possible values of number of episomes per cluster
  for(j in 2:n_iterations){
    for(k in 1:q){
      # define the probability of each nk given observed data and other parameters
      nk_like <- nks %>% sapply(log_likelihood_n, mu = mu[j-1], sigma2 = 1/tau[j-1], I = I[k]) # log likelihood
      nk_probs <- exp(nk_like - max(nk_like))/sum(exp(nk_like - max(nk_like)), na.rm = T) #probabilities
      # sample nk
      n[j, k] <- sample(nks, 1, prob = nk_probs, replace =  T)
    }
    # sample mu:
    mu[j] <- rnorm(1, sum(I/n[j,])/q, sqrt(1/(q*tau[j-1])))
    # sample tau:
    tau[j] <- rgamma(1, q/2, 0.5*sum((I/n[j,]-mu[j])^2))
  }
  
  return(cbind(iteration = 1:n_iterations, mu, tau, as.data.frame(n)))
}

#####
# Simulate data
#####

set.seed(400)

# # simulate data:
# real_mu <- 550 # based on observed data
# real_sigma2 <- 100000 # based on observed data
# real_ns <- sample(1:5, 230, replace = T, prob = dpois(1:5, 1)/sum(dpois(1:5, 1))) # based on minimum number of episomes in fixed KSHV data
# q <- length(real_ns) # number of clusters
# I <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q)

# simulate data:
# real_mu <- 925 # based on observed data
# real_sigma2 <- 100000 # based on observed data
# real_ns <- c(rep(1,105), rep(2,33), rep(3,7), rep(4,2)) # based on observed data
# q <- length(real_ns) # number of clusters
# I <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q)

names(I) <- paste0("n", 1:q)

#####
# Apply new and old samplers and compare results
#####
n_iterations <- 10000
tau0 <- 1/real_sigma2
mu0 <- 1000
# really_old_results1 <- run_gibbs_really_old(tau0, mu0, I, n_iterations)
old_results1 <- run_gibbs_old(tau0, mu0, I, n_iterations)
new_results1 <- run_gibbs_new(tau0, mu0, I, n_iterations)

tau0 <- 1e-4
mu0 <- 300
# really_old_results2 <- run_gibbs_really_old(tau0, mu0, I, n_iterations)
old_results2 <- run_gibbs_old(tau0, mu0, I, n_iterations)
new_results2 <- run_gibbs_new(tau0, mu0, I, n_iterations)


old_results <- rbind(old_results1 %>% mutate(chain = "chain1"),
                     old_results2 %>% mutate(chain = "chain2")) %>% 
  filter(iteration > 500)

new_results <- rbind(new_results1 %>% mutate(chain = "chain1"),
                     new_results2 %>% mutate(chain = "chain2")) %>% 
  filter(iteration > 500)

save(old_results, new_results, real_mu, real_sigma2, I, real_ns, file = "results/figure_out_MCMC_change_prior.RData")

p <- rbind(old_results %>% mutate(alg = "old"),
      new_results %>% mutate(alg = "new")) %>% 
      # really_old_results %>% mutate(alg = "really old")) %>% 
  ggplot(aes(iteration, mu, color= interaction(chain, alg))) + 
  geom_point(shape = NA) + 
  geom_line() + 
  geom_hline(yintercept = real_mu) + 
  ylim(c(0, max(c(old_results$mu, new_results$mu, real_mu)))) +
  theme(legend.position = "bottom")
  # facet_wrap(~alg)

ggMarginal(p, margins = "y", groupColour = T)

p2 <- rbind(old_results %>% mutate(alg = "old"),
      new_results %>% mutate(alg = "new")) %>% 
      # really_old_results %>% mutate(alg = "really old")) %>%
  ggplot(aes(iteration, 1/tau, color= interaction(chain, alg))) + 
  geom_line() + 
  geom_point(shape = NA) +
  geom_hline(yintercept = real_sigma2) +
  theme(legend.position = "bottom") +
  ylim(c(0, max(c(1/old_results$tau, 1/new_results$tau, real_sigma2)))) + 
  labs(y = "sigma^2")
  # facet_wrap(~alg)

ggMarginal(p2, margins = "y", groupColour = T)

convergence_test_old <- convergence_results(old_results)
convergence_test_new <- convergence_results(new_results)

old_results1 %>% select(starts_with("n")) %>% 
  pivot_longer(everything(), names_to = "cluster", values_to = "episomes") %>% 
  count(cluster, episomes) %>% 
  group_by(cluster) %>% 
  filter(n == max(n)) %>% 
  merge(data.frame(real_ns) %>% mutate(cluster= paste0("n", 1:nrow(.)))) %>% 
  count(real_ns == episomes)

new_results1 %>% select(starts_with("n")) %>% 
  pivot_longer(everything(), names_to = "cluster", values_to = "episomes") %>% 
  count(cluster, episomes) %>% 
  group_by(cluster) %>% 
  filter(n == max(n)) %>% 
  merge(data.frame(real_ns) %>% mutate(cluster= paste0("n", 1:nrow(.)))) %>% 
  count(real_ns == episomes)

# #####
# # Try a different mu 
# #####
# 
# set.seed(400)
# 
# # # simulate data pt 2:
# # real_mu <- 550 # based on observed data
# # real_sigma2 <- 10000 # based on observed data
# # real_ns_case2 <- sample(1:5, 147, replace = T, prob = dpois(1:5, 1)/sum(dpois(1:5, 1))) # based on minimum number of episomes in fixed KSHV data
# # q <- length(real_ns) # number of clusters
# # I_case2 <- matrix(rnorm(q, real_mu*real_ns_case2, sqrt(real_sigma2*real_ns_case2)), ncol = q)
# # 
# # names(I_case2) <- paste0("n", 1:q)
# 
# #####
# # Apply new and old samplers and compare results
# #####
# n_iterations <- 1000
# tau0 <- 1/real_sigma2
# mu0 <- 1000
# # really_old_results1 <- run_gibbs_really_old(tau0, mu0, I, n_iterations)
# old_results1_case2 <- run_gibbs_old(tau0, mu0, I_case2, n_iterations)
# new_results1_case2 <- run_gibbs_new(tau0, mu0, I_case2, n_iterations)
# 
# tau0 <- 1e-4
# mu0 <- 300
# # really_old_results2 <- run_gibbs_really_old(tau0, mu0, I, n_iterations)
# old_results2_case2 <- run_gibbs_old(tau0, mu0, I_case2, n_iterations)
# new_results2_case2 <- run_gibbs_new(tau0, mu0, I_case2, n_iterations)
# 
# 
# old_results_case2 <- rbind(old_results1_case2 %>% mutate(chain = "chain1"),
#                            old_results2_case2 %>% mutate(chain = "chain2")) %>% 
#   filter(iteration > 500)
# 
# new_results_case2 <- rbind(new_results1_case2 %>% mutate(chain = "chain1"),
#                            new_results2_case2 %>% mutate(chain = "chain2")) %>% 
#   filter(iteration > 500)
# 
# save(old_results_case2, new_results_case2, real_mu, real_sigma2, I_case2, real_ns_case2, file = "results/figure_out_MCMC_case2_change_prior.RData")
# 
# p_case2 <- rbind(old_results_case2 %>% mutate(alg = "old"),
#                  new_results_case2 %>% mutate(alg = "new")) %>% 
#   # really_old_results %>% mutate(alg = "really old")) %>% 
#   ggplot(aes(iteration, mu, color= interaction(chain, alg))) + 
#   geom_point(shape = NA) + 
#   geom_line() + 
#   geom_hline(yintercept = real_mu) + 
#   ylim(c(0, max(c(old_results$mu, new_results$mu, real_mu)))) +
#   theme(legend.position = "bottom")
# # facet_wrap(~alg)
# 
# ggMarginal(p_case2, margins = "y", groupColour = T)
# 
# p2_case2 <- rbind(old_results_case2 %>% mutate(alg = "old"),
#                   new_results_case2 %>% mutate(alg = "new")) %>% 
#   # really_old_results %>% mutate(alg = "really old")) %>%
#   ggplot(aes(iteration, 1/tau, color= interaction(chain, alg))) + 
#   geom_line() + 
#   geom_point(shape = NA) +
#   geom_hline(yintercept = real_sigma2) +
#   theme(legend.position = "bottom") +
#   ylim(c(0, max(c(1/old_results$tau, 1/new_results$tau, real_sigma2)))) + 
#   labs(y = "sigma^2")
# # facet_wrap(~alg)
# 
# ggMarginal(p2_case2, margins = "y", groupColour = T)
# 
# old_results1_case2 %>% select(starts_with("n")) %>% 
#   pivot_longer(everything(), names_to = "cluster", values_to = "episomes") %>% 
#   count(cluster, episomes) %>% 
#   group_by(cluster) %>% 
#   filter(n == max(n)) %>% 
#   merge(data.frame(real_ns_case2) %>% mutate(cluster= paste0("n", 1:nrow(.)))) %>% 
#   count(real_ns == episomes)
# 
# new_results1_case2 %>% select(starts_with("n")) %>% 
#   pivot_longer(everything(), names_to = "cluster", values_to = "episomes") %>% 
#   count(cluster, episomes) %>% 
#   group_by(cluster) %>% 
#   filter(n == max(n)) %>% 
#   merge(data.frame(real_ns_case2) %>% mutate(cluster= paste0("n", 1:nrow(.)))) %>% 
#   count(real_ns == episomes)
