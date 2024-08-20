## This script runs and plots simulations of a cell population with constant size

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
library(Hmisc)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))

#####
# Run 100 trials of simulations with varying parameters
#####
ntrials <- 100

# Number of cells to start with
n_cells <- 1

# Number of episomes per cell at the start
n_epi <- 3

# Indicator if simulations should be rerun (TRUE) or loaded from previous run (FALSE)
rerun <- TRUE

out_folder <- here("results","simulations_expo_selection")

if(!file.exists(out_folder)) dir.create(out_folder, showWarnings = F)

if(!rerun){
  load(here(out_folder, "expo_pop_3epi.RData"))

}else{
  
  ## set up a grid of parameters to simulate with and collect all results:
  pReps = 0.8
  pSegs = c(0.8, 0.9)
  # n_cells = seq(5,40, by = 5)
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
  
  averages <- expo_3epi_df %>% 
    group_by(total, Pr, Ps, trial) %>% 
    filter(episomes_per_cell != -1) %>% 
    summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
    ungroup()
  
  save(expo_3epi_df, param_grid, averages, 
       file = here(out_folder, "expo_pop_3epi.RData"))
}

expo_3epi_df %>% 
  filter(episomes_per_cell != -1) %>% 
  ggplot(aes(log10(total), log10(frac + 1/10000), color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  stat_summary(data = expo_3epi_df %>% filter(episomes_per_cell == -1),
               fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(data = expo_3epi_df %>% filter(episomes_per_cell == -1),
               fun.data="mean_cl_normal", geom="smooth", color = "black") + 
  facet_grid(Pr ~ Ps, labeller = "label_both") + 
  labs(x = "Log10(Population Size)", y = "Log10(Fraction of Population)", color = "k-episome cells, k = ...",
       caption = "Black line denotes the average number of episomes per cell")

fraction_of_pop <- expo_3epi_df %>% 
  filter(episomes_per_cell != -1) %>% 
  ggplot(aes(log10(total), frac*100, color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  facet_grid(Pr ~ Ps, labeller = "label_both") + 
  labs(x = "Log10(Population Size)", y = "Percent of Population", color = "k-episome cells, k = ...",
       title = "Percent of population with k episomes per cell",
       subtitle = "80% replication efficiency, varying segregation efficiency") 

ggsave(here(out_folder, "fraction_of_population.pdf"), fraction_of_pop, width = 10, height = 4)

average_epi <- expo_3epi_df %>% 
  filter(episomes_per_cell == -1) %>% 
  ggplot(aes(log10(total), frac, color = Ps, group = interaction(Pr, Ps))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  labs(x = "Log10(Population Size)", y = "Average", color = "Segregation Efficiency",
       title = "Average number of episomes per cell",
       subtitle = "80% replication efficiency, varying segregation efficiency") +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) + 
  scale_color_brewer(palette = "Set1") 

ggsave(here(out_folder, "average_episomes.pdf"), average_epi, width = 4, height = 4)

expo_3epi_df %>% 
  filter(episomes_per_cell == -1) %>% 
  ggplot(aes(total, frac, color = Ps, group = interaction(Pr, Ps, trial))) + 
  geom_line(alpha = 0.25) +
  # stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  # stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  labs(x = "Log10(Population Size)", y = "Average", color = "Segregation Efficiency",
       title = "Average number of episomes per cell",
       subtitle = "80% replication efficiency, varying segregation efficiency") +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) + 
  scale_color_brewer(palette = "Set1")  


######
# Figures for paper

figure_plot <- expo_3epi_df %>% 
  filter(episomes_per_cell == -1, Ps == "80%") %>% 
  ggplot(aes(log10(total), frac, group = interaction(Pr, Ps))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", color = "black") + 
  labs(x = "Log10(Population Size)", y = "Average number of episomes per cell", color = "Segregation Efficiency",
       title = "80% Replication, 80% Segregation")  +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) + 
  
  expo_3epi_df %>% 
  filter(episomes_per_cell == -1, Ps == "90%") %>% 
  ggplot(aes(log10(total), frac, group = interaction(Pr, Ps))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", color = "black") + 
  labs(x = "Log10(Population Size)", y = "Average number of episomes per cell", color = "Segregation Efficiency",
       title = "80% Replication, 90% Segregation")  +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) +
  coord_cartesian(ylim = c(0, 3)) +
  
  expo_3epi_df %>% 
  filter(episomes_per_cell != -1, Ps == "80%") %>% 
  ggplot(aes(log10(total), frac, color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  # facet_grid(Pr ~ Ps, labeller = "label_both") + 
  labs(x = "Log10(Population Size)", y = "Relative Abundance", color = "Episomes\nper cell") +
       # title = "Percent of population with k episomes per cell",
       # subtitle = "80% replication efficiency, varying segregation efficiency") +
  
  expo_3epi_df %>% 
  filter(episomes_per_cell != -1, Ps == "90%") %>% 
  ggplot(aes(log10(total), frac, color = as.factor(episomes_per_cell), group = factor(episomes_per_cell))) + 
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  # facet_grid(Pr ~ Ps, labeller = "label_both") + 
  labs(x = "Log10(Population Size)", y = "Relative Abundance", color = "Episomes\nper cell") +
       # title = "Percent of population with k episomes per cell",
       # subtitle = "80% replication efficiency, varying segregation efficiency") 
  
  plot_layout(guides = "collect") &
  scale_color_manual(values = safe_colorblind_palette) & 
  theme(legend.justification = c(0.5, -0.05)) 

ggsave(here(out_folder, "exponential_selection_figure.pdf"), figure_plot, width = 10, height = 7)


expo_3epi_df %>% 
  filter(episomes_per_cell != -1, Ps == "90%") %>% 
  group_by(total = round(log10(total), 2), Pr, Ps, episomes_per_cell) %>% 
  summarise(abundance = mean(number_of_cells)) %>% 
  ggplot(aes(total, abundance, fill = as.factor(episomes_per_cell))) + 
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Log10(Population Size)", y = "Percent of Population", color = "k-episome cells, k = ...") +
  scale_fill_manual(values = safe_colorblind_palette)

## Plot average episome over time for case with moderate defect
simple_plot <-  averages %>% 
  filter(Ps == "90%") %>% 
  ggplot(aes(log10(total), avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") +
  labs(title = "80% Replication, 90% Segregation") +
  coord_cartesian(c(0,700)) +
  
  averages %>% 
  filter(Pr == "100%", Ps == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "100% Replication, 90% Segregation") +
  coord_cartesian(c(0,700)) +
  
  averages %>% 
  filter(Pr == "90%", Ps == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "90% Replication, 100% Segregation") +
  coord_cartesian(c(0,120)) +
  
  averages %>% 
  filter(Pr == "90%", Ps == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "90% Replication, 90% Segregation") +
  coord_cartesian(c(0,120)) +
  
  plot_layout(ncol = 2, byrow = T) &
  # geom_rug(data = . %>% group_by(trial) %>% filter(time == max(time))) &
  ylim(0,5) & 
  theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", y = "Average number of episomes per cell")


ggsave(here(out_folder, "simple_averages.png"), simple_plot, width = 8, height = 8)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


trajectory_plot <- function(Pr, Ps, xlim){
  PR <- Pr
  PS <- Ps
  
  plot <-  averages %>% 
    filter(Pr == PR, Ps == PS) %>% 
    ggplot(aes(time, avg, group = trial)) + 
    geom_line(alpha = 0.25, color= "gray") 
  
  if(!(PR == "100%" & PS == "100%")){
    temp <- time_to_expo %>% 
      mutate(Pr = paste0(pRep*100, "%"),
             Ps = paste0(pSeg*100, "%")) %>% 
      filter(Pr == PR, Ps == PS) %>% 
      pull(expoTime) 
    
    mean_time <- round(mean(temp))
    sd_time <- round(sd(temp))
    
    label <- str_interp("T[E]: ${mean_time} %+-% ${sd_time}")
    
    plot <- plot + 
      geom_text(data = data.frame(NA), aes(mean_time, 0, 
                                           label = label),
                inherit.aes = F, parse = T, vjust = -0.3) +
      geom_point(data = data.frame(NA), aes(mean_time, 0),
                 inherit.aes = F) 
  }
  
  plot +
    labs(title = str_interp("${PR} Replication, ${PS} Segregation")) +
    labs(x = "Time (generations)",y = "Average number of episomes per cell") +
    ylim(c(0,5)) +
    coord_cartesian(c(0,xlim)) 
}

distribution_plot1 <- function(Pr = "100%", Ps = "100%", xlim, facet = F, Pr_facet = "all"){
  PR <- Pr
  PS <- Ps
  
  if(facet){
    df <- expo_3epi_df_long %>% 
      group_by(Pr, Ps, time = round(time), episomes_per_cell) %>% 
      summarise(number_of_cells = mean(number_of_cells))  
    if(!any(Pr_facet == "all")) df <- df %>% filter(Pr %in% Pr_facet)
  }else{
    df <- expo_3epi_df_long %>% 
      filter(Pr == PR, Ps == PS) %>% 
      group_by(time = round(time), episomes_per_cell) %>% 
      summarise(number_of_cells = mean(number_of_cells)) 
  }
  
  if(facet){
    if(Pr_facet == "100%"){
      df <- df %>% rbind(data.frame(Pr = factor("100%", levels = levels(expo_3epi_df_long$Pr)), 
                                    Ps = factor("100%", levels = levels(expo_3epi_df_long$Ps)), 
                                    time = 1:750, 
                                    episomes_per_cell = factor("3", levels = levels(expo_3epi_df_long$episomes_per_cell)), 
                                    number_of_cells = 1000))
    }
    add_zeros <- function(data, id){
      if(id$Pr == "100%" & id$Ps == "100%") return(data)
      return(rbind(data, 
                   data.frame(time = (max(data$time) + 1):750, 
                              episomes_per_cell = factor("0", levels = levels(expo_3epi_df_long$episomes_per_cell)), 
                              number_of_cells = 1000)))
    } 
    
    df <- df %>% group_by(Pr, Ps) %>% group_modify(add_zeros)
  }else if(PR == "100%" & PS == "100%"){
    df <- df %>% rbind(data.frame(time = 1:750, episomes_per_cell = factor("3", levels = levels(expo_3epi_df_long$episomes_per_cell)), number_of_cells = 1000))
  } 
  
  plot <- df %>% 
    ggplot(aes(time, number_of_cells, fill = episomes_per_cell, color = episomes_per_cell)) + 
    geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + 
    scale_fill_manual(values = safe_colorblind_palette) + 
    scale_color_manual(values = safe_colorblind_palette) + 
    coord_cartesian(c(0,xlim)) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1))  +
    labs(x = "Time (generations)", fill = "Episomes per cell", y = "Percent of\ncell population") +
    guides(color = "none") + 
    theme(panel.grid = element_blank())
  
  if(facet) plot <- plot + facet_grid(Pr ~ Ps, labeller = "label_both")
  
  plot
}

new_plot <- trajectory_plot(Pr = "100%", Ps = "100%", 700) +
  
  trajectory_plot(Pr = "100%", Ps = "90%", 700) +
  
  distribution_plot1(Pr = "100%", Ps = "100%", 700) + 
  
  distribution_plot1(Pr = "100%", Ps = "90%", 700) +
  
  trajectory_plot(Pr = "90%", Ps = "100%", 120) +
  
  trajectory_plot(Pr = "90%", Ps = "90%", 120) +
  
  distribution_plot1(Pr = "90%", Ps = "100%", 120) +
  
  distribution_plot1(Pr = "90%", Ps = "90%", 120) +
  
  plot_layout(ncol = 2, byrow = T, heights = c(3,1, 3, 1), guides = "collect") &
  # theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", fill = "Episomes per cell") &
  guides(color = "none")

ggsave(here(out_folder, "new_simulations_fig.png"), new_plot, width = 10, height = 10)

supplemental_plot <- distribution_plot1(xlim = 700, facet = T, Pr_facet = "100%") + 
  theme(axis.title.x = element_blank()) +
  distribution_plot1(xlim = 250, facet = T, Pr_facet = "95%") + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  distribution_plot1(xlim = 120, facet = T, Pr_facet = "90%") +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  distribution_plot1(xlim = 60, facet = T, Pr_facet = "75%") + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  distribution_plot1(xlim = 30, facet = T, Pr_facet = "50%")  + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  plot_layout(ncol = 1,guides = "collect")

ggsave(here(out_folder, "supplemental_episome_distribution.png"), supplemental_plot , width = 12, height = 8)

#####
# Plot histograms with distribution of episomes in population at varying time points
#####
distribution_plot2 <- function(generation, PR, PS){
  df <- expo_3epi_df_long %>% 
    filter(pRep == PR, pSeg == PS) %>% 
    group_by(trial, pRep, pSeg) %>% 
    filter((abs(time - generation)) == min(abs(time - generation))) %>% 
    mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) 
  
  # average <- df %>% group_by(trial) %>% summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) 
  
  df %>% 
    group_by(episomes_per_cell, Pr, Ps) %>% 
    summarise(percent_of_cells = mean(frac)*100) %>% 
    ggplot(aes(episomes_per_cell, percent_of_cells)) + 
    geom_bar(stat = "identity") + 
    labs(x = "Number of episomes per cell", y = "Percent of cells with a\ngiven number of episomes",
         title = paste("Division", generation)) + 
    ylim(c(0, 100)) + 
    scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(-1,10))
}

plot_epi_distribution <- function(Pr, Ps){
  plot <- distribution_plot2(1, Pr,Ps) + theme(axis.title.y = element_text(size = 10)) + 
    distribution_plot2(10, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    distribution_plot2(20, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    distribution_plot2(100, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    plot_layout(nrow = 1) &
    # plot_annotation(title = str_interp("${Pr*100}% Replication, ${Ps*100}% Segregation")) &
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(), plot.title = element_text(size = 8)) &
    scale_y_continuous(breaks = c(0, 50, 100), limits = c(0,100)) & 
    labs(y = "Percent of\ncell population")
  
  # # Use the tag label as an x-axis label
  # wrap_elements(panel = plot) +
  #   labs(tag = "Number of episomes per cell") +
  #   theme(
  #     plot.tag = element_text(size = rel(1)),
  #     plot.tag.position = "bottom"
  #   )
  plot
}

plot_epi_distribution(1,0.9)
plot_epi_distribution(0.9,1)
plot_epi_distribution(0.9,0.9)

x_label <- ggdraw() + 
  draw_label(
    "Number of episomes per cell",
    x = 0.5,
    hjust = 0.5,
    size = 10,
    y = 0.8
    # vjust = 1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(-10, 0, 0, 0),
    plot.background = element_rect(color = "white")
  )
## add epi distributions to "simple plot"
new_plot2 <- plot_grid(
  trajectory_plot(Pr = "100%", Ps = "100%", 700),
  trajectory_plot(Pr = "100%", Ps = "90%", 700),
  plot_epi_distribution(1,1) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  plot_epi_distribution(1,0.9) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  x_label, x_label,
  trajectory_plot(Pr = "90%", Ps = "100%", 120),
  trajectory_plot(Pr = "90%", Ps = "90%", 120),
  plot_epi_distribution(0.9,1) ,#+ theme(plot.margin = margin(0,1,0,1)),
  plot_epi_distribution(0.9,0.9) ,#+ theme(plot.margin = margin(0,1,0,1)),
  x_label, x_label,
  
  byrow = T, ncol = 2, rel_heights = c(2,1, 0.1, 2,1, 0.1),
  axis = "lr", align = "v"
)

ggsave(here(out_folder, "new_simulations_fig_v2.png"), new_plot2, width = 10, height = 10, bg = "white")


#####
# Supplemental figures of simulation summary
#####

time_to_expo_summary <- time_to_expo %>% 
  mutate(Pr = paste0(pRep*100, "%"),
         Ps = paste0(pSeg*100, "%")) %>% 
  filter(!(pRep == 1 & pSeg == 1)) %>% 
  group_by(n_cells, Pr, Ps, pRep, pSeg) %>% 
  summarise(min = min(expoTime), sd = sd(expoTime), expoTime = mean(expoTime))

sup1 <- time_to_expo_summary %>% 
  ggplot(aes(pSeg*100, expoTime, color = as.factor(n_cells), group = n_cells)) + 
  geom_point(size = 2) + geom_line(lwd = 1) + 
  geom_errorbar(aes(ymin = expoTime - sd, ymax = expoTime + sd), width = 0.5, lwd = 1) + 
  # add infity for Pr = Ps = 1
  # scale_color_manual(values = sort(scico_palette_data("acton", T))[c(1, 31, 51, 71, 81)]) + 
  labs(x = "Ps (%)", y = "Time to episome expo (generations)", color = "Number of cells") +
  theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))

sup2 <- time_to_expo_summary %>% 
  mutate(n_cells = ifelse(pSeg == 0.5, n_cells - 0.01, 
                          ifelse(pSeg == 0.75, n_cells - 0.005, 
                                 ifelse(pSeg == 0.95, n_cells + 0.005, 
                                        ifelse(pSeg == 1, n_cells + 0.01, n_cells))))) %>% 
  ggplot(aes(n_cells, expoTime, color = as.factor(pSeg*100), group = pSeg)) + 
  geom_line(lwd = 1) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = expoTime - sd, ymax = expoTime + sd), width = 0.5, lwd = 1) + 
  # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 41, 61, 81)]) + 
  # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) + 
  labs(x = "Number of cells", y = "Time to episome expo (generations)", color = "Ps (%)") +
  theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))


sup3 <- averages %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>%
  summarise(avg = mean(avg)) %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps)))) %>%
  ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line(size = 1) +
  facet_wrap(~Pr, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell", color = "Ps (%)") + 
  # add time to expo annotation
  geom_point(data = time_to_expo_summary %>%
               mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                      Pr = factor(Pr, levels = levels(averages$Pr))), 
             aes(expoTime, 0), shape = 4, show.legend = F) +
  geom_rect(aes(xmin = 1850, xmax = 2250, ymin = 0.5, ymax = 3.2), fill = "white", alpha = 0.05, color = NA) +
  geom_text(data = time_to_expo_summary %>%
              mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                     Pr = factor(Pr, levels = levels(averages$Pr)),
                     y = ifelse(Ps == 100, 3, ifelse(Ps == 95, 2.5, ifelse(Ps == 90, 2, ifelse(Ps == 75, 1.5, 1))))), 
            aes(1850, y, label = paste("T[E]:", round(expoTime), "%+-%", round(sd))), 
            parse = T, hjust = 0, vjust = 1, 
            show.legend = F) 

sup <- ((plot_spacer() / sup1 / plot_spacer() / sup2 / plot_spacer() + plot_layout(heights = c(0.5,4,0.5,4,0.5))) | 
          plot_spacer() | 
          sup3) + plot_layout(widths = c(1,0.1,1))

ggsave(here(out_folder, "supplemental_figure.png"), sup , width = 10.5, height = 10)


sup3.2 <- averages %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>%
  summarise(avg = mean(avg)) %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))),
         Pr = fct_rev(Pr)) %>%
  ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line(size = 1) +
  facet_wrap(~Pr, nrow = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell", color = "Ps (%)") + 
  # add time to expo annotation
  geom_point(data = time_to_expo_summary %>%
               mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                      Pr = factor(Pr, levels = levels(averages$Pr))), 
             aes(expoTime, 0), shape = 4, show.legend = F) +
  geom_rect(aes(xmin = 1850, xmax = 2250, ymin = 0.5, ymax = 3.2), fill = "white", alpha = 0.05, color = NA) +
  geom_text(data = time_to_expo_summary %>%
              mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                     Pr = factor(Pr, levels = levels(averages$Pr)),
                     y = ifelse(Ps == 100, 3, ifelse(Ps == 95, 2.5, ifelse(Ps == 90, 2, ifelse(Ps == 75, 1.5, 1))))), 
            aes(1850, y, label = paste("T[E]:", round(expoTime), "%+-%", round(sd))), 
            parse = T, hjust = 0, vjust = 1, 
            show.legend = F)  + 
  theme(legend.position = "bottom")

sup.2 <-( (sup1 + plot_spacer() + sup2 + plot_layout(widths = c(1, 0.1, 1), nrow = 1)) / 
            plot_spacer() / sup3.2) + plot_layout(heights = c(1,0.1,1), ncol = 1)

ggsave(here(out_folder, "supplemental_figure_v2.png"), sup.2 , width = 15, height = 7)


## Calculate shannon index
shannon <- expo_3epi_df_long %>% 
  filter(number_of_cells != 0) %>% 
  mutate(p = number_of_cells/total_cells) %>% 
  group_by(pRep, pSeg, trial, time) %>% 
  summarise(H = -sum(p*log(p)))

shannon %>% 
  mutate(Pr = pRep*100, Ps = pSeg*100) %>% 
  filter(Pr == 100, Ps == 95) %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>% 
  summarise(avg = mean(H)) %>% head
ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line() +
  # facet_wrap(~Ps, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average shannon index", color = "Ps (%)") 



expo_3epi_df_long %>% filter(pSeg == 0.9, pRep == 0.9) %>% 
  filter(episomes_per_cell == 0) %>% 
  mutate(latent_cells = total_cells - number_of_cells) %>% 
  ggplot(aes(time, latent_cells)) + 
  geom_line(aes(group = trial), alpha = 0.1, color = "gray")  + 
  # geom_smooth(formula = (log(y)~x), method = "lm")
  geom_smooth()
# scale_y_log10()
# facet_grid(Pr ~ Ps, labeller = "label_both")

