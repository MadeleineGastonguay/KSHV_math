# This script visualizes experimental data from longitdunal experiments in Brk219

### Setup ###################################################################### 

library(here)
library(tidyverse)
library(fitdistrplus)
theme_set(theme_minimal())

### Script inputs ##############################################################
## Read in data:
LANA_dots <- read_csv(here("data", "derived", "brk219_LANA_dots.csv"))
GFP <- read_csv(here("data", "derived", "brk219_longitudinal.csv"))
cell_growth <- read_csv(here("data", "derived", "brk219_cell_growth.csv"))
percent_dead <- read_csv(here("data", "derived", "brk219_percent_dead.csv"))


### Plot cell growth ###########################################################

cell_growth %>% 
  ggplot(aes(day, log10(live_cells*1e5))) + 
  geom_point() + 
  geom_line()

# Define passaging times
cut_times <- c(4.5,9.5,16,21.5,26.5,31.5,38.5,43.5)

percent_dead_plot <- percent_dead %>% 
  pivot_longer(!day) %>% 
  ggplot(aes(day, value, color = name)) + 
  geom_point(size = 3) + 
  geom_line(linewidth = 1) + 
  scale_color_brewer(palette = "Set2") + 
  geom_vline(data = data.frame(t = cut_times), aes(xintercept = t)) + 
  coord_cartesian(ylim = c(0,100)) + 
  labs(x = "Time (days)", y = "Percent dead")

### Plot percent of cells with episomes remaining ##############################
GFP %>% 
  ggplot(aes(day, percent_GFP_live)) + 
  geom_point() + 
  geom_line() + 
  ylim(c(0,100)) + 
  labs(x = "Time (day)", y= "Percent of cells with episomes")
  

### Plot LANA_dot ##############################################################
LANA_summary <- LANA_dots %>% 
  group_by(day) %>% 
  summarise(mean = mean(LANA_dots),
            median = median(LANA_dots),
            sd = sd(LANA_dots),
            var = var(LANA_dots),
            min = min(LANA_dots), 
            max = max(LANA_dots),
            n = n(),
            total = sum(LANA_dots))

LANA_dots %>% 
  mutate(day_f = factor(as.character(day), levels = as.character(0:12))) %>% 
  ggplot(aes(day_f, LANA_dots)) + 
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm() + 
  geom_line(data = LANA_summary, aes(day+1, mean, color = "mean")) + 
  geom_line(data = LANA_summary, aes(day+1, median, color = "median")) + 
  scale_x_discrete(drop = F) 

### Estimate variation in initial number of episomes from LANA dots at day 0 ###

# Fit multiple distributions
day0_LANA <- LANA_dots %>% filter(day == 0) %>% pull(LANA_dots)
fit_pois <- fitdist(day0_LANA, "pois")
fit_nb <- fitdist(day0_LANA, "nbinom")

# Compare 
par(mfrow = c(2, 2))
denscomp(list(fit_pois, fit_nb), 
         legendtext = c("Poisson", "Negative Binomial"))
qqcomp(list(fit_pois, fit_nb), 
       legendtext = c("Poisson", "Negative Binomial"))
cdfcomp(list(fit_pois, fit_nb), 
        legendtext = c("Poisson", "Negative Binomial"))
ppcomp(list(fit_pois, fit_nb), 
       legendtext = c("Poisson", "Negative Binomial"))

# Negative ninimal has better BIC:
gofstat(list(fit_pois, fit_nb), 
        fitnames = c("Poisson", "Negative Binomial"))

# Compare each distribution to observed
LANA_dots %>% 
  filter(day == 0) %>% 
  ggplot(aes(LANA_dots)) + 
  geom_density(aes(fill = "observed"), alpha = 0.5) + 
  geom_density(data = tibble(LANA_dots = rpois(31, fit_pois$estimate)), aes(fill = "poisson"), alpha = 0.5) + 
  geom_density(data = tibble(LANA_dots = rnbinom(31, size = fit_nb$estimate[1], mu = fit_nb$estimate[2])), 
               aes(fill = "negative binomial"), alpha = 0.5)

# Poisson doesn't capture full variance
LANA_dots %>% 
  filter(day == 0) %>% 
  mutate(dist = "observed") %>% 
  dplyr::select(LANA_dots, dist) %>% 
  rbind(
    tibble(LANA_dots = rpois(31, fit_pois$estimate)) %>% mutate(dist = "poisson"),
    tibble(LANA_dots = rnbinom(31, size = fit_nb$estimate[1], mu = fit_nb$estimate[2])) %>% mutate(dist = "negative binomial")
  ) %>% 
  ggplot(aes(LANA_dots, dist)) + 
  geom_boxplot(aes(fill = dist), outlier.shape = NA, alpha = 0.5) + 
  ggbeeswarm::geom_beeswarm(aes(color = dist), size = 2)

### Plot rate of average GFP decay and LANA dots ###############################

percent_df <- rbind(
  GFP %>% 
    select(day, median_GFP_live, mean_GFP_live) %>% 
    pivot_longer(!day) %>% mutate(data = "GFP"),
  LANA_summary %>% 
    select(day, mean, median) %>% 
    pivot_longer(!day) %>% mutate(data = "LANA dots")
) %>% 
  filter(!is.na(value)) %>% 
  group_by(data, name) %>% 
  mutate(percent = value/max(value)*100,
         name = str_remove(name, "_GFP_live")) %>% 
  ungroup %>% 
  filter(name != "median") 

percent_df %>% 
  ggplot(aes(day, log10(percent)-2, shape = data)) +
  geom_smooth(method = "lm", aes(color = data), se= F, formula = y ~ x + 0) + 
  geom_point(size = 2.5, alpha = 0.75) + 
  geom_line() +
  labs(x = "Time (days)", y = "Percent of signal remaining", title = "Mean LANA dots or GFP intensity across live cells") + 
  scale_y_continuous(limits = c(0, 2)-2, labels = function(x) 10^(x+2), breaks = c(0,1,2)-2)  +
  scale_color_manual(values = palette.colors()[2:3]) + 
  theme(legend.position = "inside", legend.position.inside = c(0.8,0.5), legend.justification = c(1,1))


### Fit rate of GFP loss and LANA dot loss #####################################

GFP_expo_fit <- nls(mean_GFP_live ~ a * exp(-b * day), data = GFP, 
           start = list(a = 4500000, b = 0.1))

GFP_linear_fit <- nls(mean_GFP_live ~ a - b * day, data = GFP, 
                start = list(a = 4500000, b = 0.1))

BIC(GFP_expo_fit)
BIC(GFP_linear_fit)

LANA_expo_fit <- nls(LANA_dots ~ a * exp(-b * day), data = LANA_dots, 
                    start = list(a = 25, b = 0.1))

LANA_linear_fit <- nls(LANA_dots ~ a - b * day, data = LANA_dots, 
                     start = list(a = 25, b = 0.1))

BIC(LANA_expo_fit)
BIC(LANA_linear_fit)

mean_LANA_expo_fit <- nls(mean ~ a * exp(-b * day), data = LANA_summary, 
                     start = list(a = 25, b = 0.1))

mean_LANA_linear_fit <- nls(mean ~ a - b * day, data = LANA_summary, 
                       start = list(a = 25, b = 0.1))

BIC(mean_LANA_expo_fit)
BIC(mean_LANA_linear_fit)
