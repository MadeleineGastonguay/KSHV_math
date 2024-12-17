# Systematically determine sensitivity to prior 

### Setup ###################################################################### 

library(tidyverse)
library(here)
library(furrr)
library(scico)
library(patchwork)

source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

plan(multisession(workers = 8))

theme_set(theme_classic())

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


### Compare possible prior distributions #######################################

# Use geometric distribution with various values for p
geom_params <- expand_grid(ns = 1:100, param = seq(0.1,0.9,by= 0.1)) %>% 
  group_by(param) %>% 
  mutate(p = dgeom(ns, param)/sum(dgeom(1:100, param)))

# Compare geometric priors to poisson prior
prior_plot <- geom_params %>% 
  filter(param <= 0.6) %>% 
  ggplot(aes(ns, p, color = paste0("geom(", param, ")"), group = param)) + 
  geom_line(size = 1) + 
  # add poisson prior
  geom_line(aes(ns, dpois(ns, 1), color = "poisson(1)"), size = 1) + 
  scale_x_continuous(limits = c(0,10), breaks = 0:10) + 
  labs(color = "prior for n_k", #caption = "poisson(1) in black for reference", 
       x = "number of episomes in a cluster", y= "prior probability")  + 
  scale_color_manual(values = c(safe_colorblind_palette[1:6], "black")) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) 

ggsave(here("results","supplemental_figures_updated_pdf_n_prior", "n_prior_plot.png"), prior_plot, width = 4, height = 3)

### Simulate data similar to fixed 8TR conditions ##############################
## Assume there are likely 2 episomes per cluster
set.seed(400)

## Define parameters of intensity distribution for a single episome
real_mu <- 525 # based on observed data
real_sigma2 <- 140^2 # based on observed data
q <- 230 # number of measured LANA dots

## sample real number of episomes per cluster
real_ns <- sample(1:5, q, replace = T, prob = dpois(1:5, 2)/sum(dpois(1:5, 2))) 
hist(real_ns)

## simulate intensity of each cluster 
I <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q) 
names(I) <- paste0("n", 1:q)


### Apply MCMC #################################################################
n_iterations <- 10000 # number of iterations for MCMC

## Try various initial conditions for MCMC
hyper_parameters <- expand_grid(tau0 = c(2/real_sigma2, 0.5/real_sigma2), 
                                mu0 = c(300, 1000)) %>% 
  mutate(chain = 1:nrow(.))

## Apply MCMC with each prior
cat('starting\n')
# Poisson
baseline <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("pois", 1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done with poisson prior\n')

geom0.1 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.2 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.2), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.3 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.3), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.4 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.4), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.5 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.5), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.6 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.6), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

### Plot results ###############################################################

all_runs <- rbind(
  baseline %>% mutate(prior = "poisson(1)"),
  geom0.1 %>% mutate(prior = "geom(0.1)"),
  geom0.2 %>% mutate(prior = "geom(0.2)"),
  geom0.3 %>% mutate(prior = "geom(0.3)"),
  geom0.4 %>% mutate(prior = "geom(0.4)"),
  geom0.5 %>% mutate(prior = "geom(0.5)"),
  geom0.6 %>% mutate(prior = "geom(0.6)")
) %>% 
  filter(iteration > 1000) 

inferred_param <- all_runs %>% group_by(prior) %>% 
  summarise(mu = DescTools::Mode(round(mu)), sigma = DescTools::Mode(round(sqrt(1/tau))))

## Trace plots of mu with each prior
mu_plot <- all_runs %>% 
  filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)"))) %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior, nrow = 1) + 
  geom_hline(yintercept = real_mu) + 
  geom_text(data = data.frame(NA), aes(1100, real_mu, label = real_mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param %>% filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)"))), 
             aes(yintercept = mu), lty = "dashed") +
  geom_text(data = inferred_param %>% filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)"))), 
            aes(1100, mu, label = mu), 
            hjust = 0, vjust = 1.1, inherit.aes = F) +
  ylim(c(0, max(all_runs$mu))) + 
  labs(y = bquote(mu), x = "Iteration")

## Trace plots of sigma with each prior
sigma_plot <- all_runs %>% 
  filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)")))  %>% 
  ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
  geom_line() + 
  facet_wrap(~prior, nrow = 1) + 
  geom_hline(yintercept = sqrt(real_sigma2)) + 
  geom_text(data = data.frame(NA), aes(1100, sqrt(real_sigma2), label = sqrt(real_sigma2)), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param %>% filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)"))) , aes(yintercept = sigma), lty = "dashed") +
  geom_text(data = inferred_param %>% filter(!(prior %in% c("geom(0.2)", "geom(0.3)", "poisson(1)"))) , aes(1100, sigma, label = sigma), 
            hjust = 0, vjust = 1.1, inherit.aes = F) +
  ylim(c(0, max(sqrt(1/all_runs$tau)))) + 
  labs(y = bquote(sigma), x = "Iteration")


sensitivity_plot <- mu_plot + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  sigma_plot + 
  theme(strip.background = element_blank(), strip.text = element_blank()) + 
  plot_layout(ncol = 1, guides = "collect")

ggsave(here("results", "supplemental_figures", "prior_sensitivity.png"), sensitivity_plot, width = 8, height = 4.5)

### Check convergence heuristics for all cases #################################
convergence <- all_runs %>% 
  group_split(prior) %>% 
  future_map_dfr(function(x) convergence_results(select(x, -prior)) %>% mutate(prior = unique(x$prior)))

convergence %>% 
  filter(prior == "geom(0.1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence %>% 
  filter(prior == "geom(0.2)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "geom(0.3)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "geom(0.4)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence %>% 
  filter(prior == "geom(0.5)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "poisson(1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  plot_layout(ncol = 3, byrow = T, guides = "collect") &
  geom_point(color="gray") &
  geom_point(data = . %>% filter(name %in% c("mu", "tau")) ,
             aes(color = as.character(name))) &
  facet_wrap(~prior, scales = "free_y") & 
  geom_vline(xintercept = 1, lty = "dashed") &
  xlim(c(0.99, 1.05)) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) &
  labs(color = "")

### Compare estimates of episomes per cluster to synthetic data ################
cluster_modes <- all_runs %>%
  filter(chain == "1") %>% 
  select(-mu, -tau, -iteration, -chain) %>% 
  pivot_longer(starts_with("n"), values_to = "n_epi", names_to = "cluster_id") %>% 
  count(prior, cluster_id, n_epi) %>% 
  group_by(prior, cluster_id) %>% 
  filter(n == max(n)) %>% 
  select(-n) %>% 
  ungroup() %>% 
  left_join(data.frame(real_ns) %>% mutate(cluster_id = paste0("n", 1:nrow(.)))) 

cluster_modes %>% 
  count(prior, n_epi, real_ns) %>% 
  ggplot(aes(real_ns, n_epi, fill = n)) + 
  geom_tile() + 
  facet_wrap(~prior) + 
  labs(x = "real nk", y= "estimated nk") + 
  geom_abline(lty = "dashed")

cluster_modes %>% 
  mutate(correct = 
           factor(ifelse(n_epi == real_ns, "correct inference", "incorrect inference"),
                  levels = c("incorrect inference", "correct inference"))) %>% 
  ggplot(aes(n_epi, fill = correct)) + 
  geom_bar(position = "stack") +
  facet_wrap(~prior) + 
  geom_bar(aes(real_ns, color = "simulated data"), fill = NA, inherit.aes = F) + 
  scale_color_manual(values = c("black") )+ 
  scale_fill_manual(values = c("red", "lightgreen")) + 
  labs(color = "", fill = "") + 
  labs(x = "number of episomes per cluster", y = "count") + 
  geom_text(data = .  %>% count(prior, correct) %>% group_by(prior) %>% mutate(p = n/sum(n)) %>% 
              filter(correct == "correct inference"),
            aes(Inf, Inf, label = paste0("accuracy: ", round(p*100), "%")),
            hjust = 1, vjust = 1.1, inherit.aes = F, show.legend = F) 

sds <- all_runs %>% 
  filter(chain == "1") %>% 
  pivot_longer(starts_with("n"), names_to = "cluster_id") %>% 
  group_by(prior, cluster_id) %>% 
  summarise(sd = sd(value))

tibble(intensity = I[1,], cluster_id = names(I)) %>% 
  left_join(cluster_modes) %>% 
  left_join(sds) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~prior, scales = "free_y") +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution")



### Repeat analysis with greater variance in simulated data ####################

## Simulate data
set.seed(400)
real_sigma2 <- 350^2 
I_lv <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q)
names(I_lv) <- paste0("n", 1:q)

## Apply MCMC
cat('starting\n')
baseline_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("pois", 1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done with baseline\n')

geom0.1_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.2_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.2), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.3_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.3), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.4_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.4), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.5_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.5), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.6_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.6), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

## Plot results
all_runs_lv <- rbind(
  baseline_lv %>% mutate(prior = "poisson(1)"),
  geom0.1_lv %>% mutate(prior = "geom(0.1)"),
  geom0.2_lv %>% mutate(prior = "geom(0.2)"),
  geom0.3_lv %>% mutate(prior = "geom(0.3)"),
  geom0.4_lv %>% mutate(prior = "geom(0.4)"),
  geom0.5_lv %>% mutate(prior = "geom(0.5)"),
  geom0.6_lv %>% mutate(prior = "geom(0.6)")
) %>% 
  filter(iteration > 1000) 

inferred_param_lv <- all_runs_lv %>% group_by(prior) %>% 
  summarise(mu = DescTools::Mode(round(mu)), sigma = DescTools::Mode(round(sqrt(1/tau))))

all_runs_lv %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_mu) + 
  geom_text(data = data.frame(NA), aes(1100, real_mu, label = real_mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param_lv, aes(yintercept = mu), lty = "dashed") +
  geom_text(data = inferred_param_lv, aes(1100, mu, label = mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(all_runs_lv$mu))) + 
  labs(y = "mu")

all_runs_lv %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_sigma2) + 
  ylim(c(0, max(1/all_runs$tau))) + 
  labs(y = "sigma^2")

all_runs_lv %>% 
  ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = sqrt(real_sigma2)) + 
  geom_text(data = data.frame(NA), aes(1100, sqrt(real_sigma2), label = sqrt(real_sigma2)), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param_lv, aes(yintercept = sigma), lty = "dashed") +
  geom_text(data = inferred_param_lv, aes(1100, sigma, label = sigma), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(sqrt(1/all_runs_lv$tau))+10)) + 
  labs(y = "sigma")


cluster_modes_lv <- all_runs_lv %>%
  filter(chain == "1") %>% 
  select(-mu, -tau, -iteration, -chain) %>% 
  pivot_longer(starts_with("n"), values_to = "n_epi", names_to = "cluster_id") %>% 
  count(prior, cluster_id, n_epi) %>% 
  group_by(prior, cluster_id) %>% 
  filter(n == max(n)) %>% 
  select(-n) %>% 
  ungroup() %>% 
  left_join(data.frame(real_ns) %>% mutate(cluster_id = paste0("n", 1:nrow(.)))) 

cluster_modes_lv %>% 
  count(prior, n_epi, real_ns) %>% 
  ggplot(aes(real_ns, n_epi, fill = n)) + 
  geom_tile() + 
  facet_wrap(~prior) + 
  labs(x = "real nk", y= "estimated nk") + 
  geom_abline(lty = "dashed")

cluster_modes_lv %>% 
  mutate(correct = 
           factor(ifelse(n_epi == real_ns, "correct inference", "incorrect inference"),
                  levels = c("incorrect inference", "correct inference"))) %>% 
  ggplot(aes(n_epi, fill = correct)) + 
  geom_bar(position = "stack") +
  facet_wrap(~prior) + 
  geom_bar(aes(real_ns, color = "simulated data"), fill = NA, inherit.aes = F) + 
  scale_color_manual(values = c("black") )+ 
  scale_fill_manual(values = c("red", "lightgreen")) + 
  labs(color = "", fill = "") + 
  labs(x = "number of episomes per cluster", y = "count") + 
  geom_text(data = .  %>% count(prior, correct) %>% group_by(prior) %>% mutate(p = n/sum(n)) %>% 
              filter(correct == "correct inference"),
            aes(Inf, Inf, label = paste0("accuracy: ", round(p*100), "%")),
            hjust = 1, vjust = 1.1, inherit.aes = F, show.legend = F) 

convergence_lv <- all_runs_lv %>% 
  group_split(prior) %>% 
  future_map_dfr(function(x) convergence_results(select(x, -prior)) %>% mutate(prior = unique(x$prior)))

convergence_lv %>% 
  filter(prior == "geom(0.1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence_lv %>% 
  filter(prior == "geom(0.2)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "geom(0.3)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "geom(0.4)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence_lv %>% 
  filter(prior == "geom(0.5)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "poisson(1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  plot_layout(ncol = 3, byrow = T, guides = "collect") &
  geom_point(color="gray") &
  geom_point(data = . %>% filter(name %in% c("mu", "tau")) ,
             aes(color = as.character(name))) &
  facet_wrap(~prior, scales = "free_y") & 
  geom_vline(xintercept = 1, lty = "dashed") &
  xlim(c(0.99, 1.05)) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) &
  labs(color = "")


sds_lv <- all_runs_lv %>% 
  filter(chain == "1") %>% 
  pivot_longer(starts_with("n"), names_to = "cluster_id") %>% 
  group_by(prior, cluster_id) %>% 
  summarise(sd = sd(value))

tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
  left_join(cluster_modes_lv) %>% 
  left_join(sds_lv) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~prior, scales = "free_y") +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution")


rbind(
tibble(intensity = I[1,], cluster_id = names(I)) %>%
  left_join(cluster_modes) %>% 
  left_join(sds) %>% mutate(data = "baseline simulation"),

tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
  left_join(cluster_modes_lv) %>% 
  left_join(sds_lv) %>% mutate(data = "simulation with larger variance")
) %>% 
  filter(prior == "geom(0.5)" | prior == "geom0.5") %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~data, ncol = 1) +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution") +
  theme_bw()

rbind(
  tibble(intensity = I[1,], cluster_id = names(I)) %>%
    left_join(cluster_modes) %>% 
    left_join(sds) %>% mutate(data = "baseline simulation") %>% 
    filter(prior == "geom0.5") %>% 
    left_join(inferred_param),
  
  tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
    left_join(cluster_modes_lv) %>% 
    left_join(sds_lv) %>% mutate(data = "simulation with larger variance") %>% 
    filter(prior == "geom(0.5)") %>% 
    left_join(inferred_param_lv)
) %>% 
  mutate(intensity = intensity/mu) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~data, ncol = 1) +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data Normalized by Inferred Mean", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution") +
  scale_x_continuous(breaks = 0:10) + 
  theme_bw()
  



