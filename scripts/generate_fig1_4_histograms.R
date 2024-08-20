## This script re-does histograms from Franceline to match the style of the rest of the functions 

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)

theme_set(theme_bw(base_size = 15))


#####
# Make data
#####

# Figure 1: p8TR LANA dot distributions:

fig1_data <- tibble(
  LANA_dots = 0:30, 
  percent = c(4.8, 36.5, 24, 19.2, 2.9, 0, 1, 1, 0, 0, 1, 0, 0, 1.9, 1, 1.9, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
)

fig1 <- fig1_data %>% 
  # filter(LANA_dots < 29) %>% 
  mutate(LANA_dots = fct_inorder(as.character(LANA_dots))) %>%
  ggplot(aes(LANA_dots, percent)) + 
  geom_col() + 
  geom_text(aes(label = sprintf("%0.1f", percent)), vjust = -0.3) +
  labs(x = "Number of LANA dots per nucleus",
       y = "Percentage of nuclei",
       title = "LANA dot distribution 14 weeks pi") + 
  theme(plot.title = element_text(size = 20, face = 2, hjust = 0.5)) +
        # axis.title = element_text(size = 15)) + 
  ylim(c(0,40))

ggsave(here("results", "figure1_histogram.png"), fig1, width = 10, height = 4)

# Figure 2: KSHV LANA dot distribution:
fig2_data <- tibble(
  LANA_dots = 0:7,
  percent = c(05.1, 80.1, 11.8, 1.5, 1.5, 0, 0, 0)
)


fig2 <- fig2_data %>% 
  # filter(LANA_dots < 5) %>% 
  mutate(LANA_dots = fct_inorder(as.character(LANA_dots))) %>%
  ggplot(aes(LANA_dots, percent)) + 
  geom_col() + 
  geom_text(aes(label = sprintf("%0.1f", percent)), vjust = -0.3) +
  labs(x = "Number of LANA dots per nucleus",
       y = "Percentage of nuclei",
       title = "Episome distribution") + 
  theme(plot.title = element_text(size = 20, face = 2, hjust = 0.5)) +
  ylim(c(0,100))

ggsave(here("results", "figure2_histogram.png"), fig2, width = 7, height = 3)

