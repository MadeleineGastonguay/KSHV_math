## This script runs and plots simulations of a cell population with constant size

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))

#####
# Run 100 trials of simulations with varying parameters
#####
ntrials <- 100
n_epi <- 3 

# Indicator if simulations should be rerun (TRUE) or loaded from previous run (FALSE)
rerun <- FALSE

out_folder <- here("results","simulations_updated")

if(!file.exists(out_folder)) dir.create(out_folder, showWarnings = F)

## set up a grid of parameters to simulate with and collect all results:
pReps = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
pSegs = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
param_grid <- expand_grid(pRep = pReps, pSeg = pSegs)

if(!rerun){
  load(here(out_folder, "constant_pop_3epi.RData"))
  
  extinction_3epi_df <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 
  
  extinction_3epi_df_long <- extinction_3epi_df %>% 
    # Thin results for memory issues
    filter(pRep %in% c(0.5, 0.8, 0.9, 0.95, 1), 
           pSeg %in% c(0.5, 0.8, 0.9, 0.95, 1)) %>% 
    pivot_extinction() %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
}else{
  
  extinction_3epi <- param_grid %>% 
    pmap(extinction, nTrials = ntrials, n_epi = n_epi)
  
  save(extinction_3epi,
       file = here(out_folder, "constant_pop_3epi.RData"))
  
  extinction_3epi_df <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 
  
  extinction_3epi_df_long <- extinction_3epi_df  %>% 
    filter(pRep %in% c(0.5, 0.8, 0.9, 0.95, 1), pSeg %in% c(0.5, 0.8, 0.9, 0.95, 1)) %>% 
    pivot_extinction() %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
  
}

# Plot average number of episomes per cell over time
averages <- extinction_3epi_df_long %>% 
  mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) %>% 
  group_by(time = round(time,1), Pr, Ps, trial) %>% 
  summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
  ungroup()

# fill out Pr = Ps = 100%
averages <- rbind(averages, 
                  data.frame(Pr = factor("100%", levels = levels(averages$Pr)), 
                             Ps = factor("100%", levels = levels(averages$Ps)), 
                             avg = 3, expand.grid(time = 41:(max(averages$time) + 30), trial = 1:100))
)

time_to_extinction <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$ExtinctionTime))



# Plot distribution of episomes per cell over time
p1 <- extinction_3epi_df_long %>% 
  ggplot(aes(time, frac, color = episomes_per_cell, group = interaction(trial, episomes_per_cell))) + 
  geom_line(alpha = 0.25) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
  labs(x = "time (generations)", y = "Fraction of\npopulation", color = "episomes per cell")

ggsave(here(out_folder, "extinction_3epi_lines.png"), p1, width = 8, height = 7)


p2.1 <- averages %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color = "gray") +
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
  labs(x = "time (generations)", y = "Average number of episomes per cell") 

ggsave(here(out_folder, "extinction_3epi_average_epi_freex.png"), p2.1, width = 8, height = 7)


p2.3 <- averages %>% 
  filter(Pr == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell")  +
  # theme(axis.title.x = element_blank()) +
  
  averages %>% 
  filter(Pr == "95%") %>% 
  ggplot(aes(time, avg, group = trial)) +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell")  +
  # axis.title.x = element_blank()) +
  
  averages %>% 
  filter(Pr == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell")  +
  # axis.title.x = element_blank()) +
  
  averages %>% 
  filter(Pr == "80%") %>% 
  ggplot(aes(time, avg, group = trial)) +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell")  +
  # axis.title.x = element_blank()) +
  
  averages %>% 
  filter(Pr == "50%") %>% 
  ggplot(aes(time, avg, group = trial)) +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell") +
  
  plot_layout(ncol = 1, axis_titles = "collect") &
  geom_line(alpha = 0.25, color = "gray") &
  facet_grid(Pr ~ Ps, labeller = "label_both") 


ggsave(here(out_folder, "extinction_3epi_average_epi_xaxis_by_Pr.png"), p2.3, width = 8, height = 7)


p2.2 <- averages %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color = "gray") +
  facet_grid(Pr ~ Ps, labeller = "label_both") +
  scale_x_log10() +
  labs(x = "time (generations)", y = "Average number of episomes per cell") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

ggsave(here(out_folder, "extinction_3epi_average_epi_logscale.png"), p2.2, width = 9, height = 7)

# Plot histograms of time to extinction


p3.1 <- time_to_extinction %>% 
  rename(Pr = pRep, Ps = pSeg) %>% 
  filter(Ps*Pr != 1) %>%
  arrange(desc(Pr), Ps) %>% 
  mutate(Pr = paste0(Pr*100, "%"),
         Ps = paste0(Ps*100, "%"),
         Pr = fct_inorder(factor(Pr)),
         Ps = fct_inorder(factor(Ps))) %>% 
  ggplot(aes(ExtinctionTime)) + 
  geom_histogram(bins = 30, aes(y = after_stat(density))) + 
  facet_grid(Pr ~ Ps, labeller = "label_both", scales = "free_y") +
  labs(x = "time to episome extinction (generations)")

ggsave(here(out_folder, "time_to_exinction_3epi.png"), p3.1, width = 7, height = 7)

p3.2 <- time_to_extinction %>% 
  rename(Pr = pRep, Ps = pSeg) %>% 
  filter(Ps*Pr != 1) %>%
  arrange(desc(Pr), Ps) %>% 
  mutate(Pr = paste0(Pr*100, "%"),
         Ps = paste0(Ps*100, "%"),
         Pr = fct_inorder(factor(Pr)),
         Ps = fct_inorder(factor(Ps))) %>% 
  ggplot(aes(ExtinctionTime)) + 
  geom_histogram(bins = 30, aes(y = after_stat(density))) + 
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free", independent = "x") +
  labs(x = "time to episome extinction (generations)")

ggsave(here(out_folder, "time_to_exinction_3epi_freex.png"), p3.2, width = 7, height = 7)

save(extinction_3epi, param_grid, 
     averages, time_to_extinction,
     file = here(out_folder, "constant_pop_3epi.RData"))



## Plot average episome over time for case with moderate defect
simple_plot <-  averages %>% 
  filter(Pr == "100%", Ps == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") +
  labs(title = "100% Replication, 100% Segregation") +
  coord_cartesian(c(0,550)) +
  
  averages %>% 
  filter(Pr == "100%", Ps == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "100% Replication, 90% Segregation") +
  coord_cartesian(c(0,550)) +
  
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
    temp <- time_to_extinction %>% 
      mutate(Pr = paste0(pRep*100, "%"),
             Ps = paste0(pSeg*100, "%")) %>% 
      filter(Pr == PR, Ps == PS) %>% 
      pull(ExtinctionTime) 
    
    mean_time <- round(mean(temp))
    sd_time <- round(sd(temp))
    
    label <- str_interp("T[E]: ${mean_time} %+-% ${sd_time}")
    
    plot <- plot + 
      geom_text(data = data.frame(NA), aes(mean_time, 0, 
                                           label = label),
                inherit.aes = F, parse = T, vjust = -0.3) +
      geom_point(data = data.frame(NA), aes(mean_time, 0),
                 inherit.aes = F) 
  }else{
    plot <- plot + geom_line(alpha = 0.25, color= "gray", linewidth = 2) 
  }
  
  plot +
    labs(title = str_interp("${PR} Replication, ${PS} Segregation")) +
    labs(x = "Time (generations)",y = "Average number of episomes per cell") +
    ylim(c(0,5)) +
    coord_cartesian(c(0,xlim)) 
}

distribution_plot1 <- function(Pr = "100%", Ps = "100%", xlim, facet = F, Pr_facet = "all", abundance = F){
  # browser()
  PR <- Pr
  PS <- Ps
  
  if(facet){
    df <- extinction_3epi_df_long %>% 
      mutate(time = ifelse(time < 0.1, 0.1, round(time ,1))) %>% 
      group_by(Pr, Ps, time, episomes_per_cell) %>% 
      summarise(number_of_cells = mean(number_of_cells))  
    if(!any(Pr_facet == "all")) df <- df %>% filter(Pr %in% Pr_facet)
  }else{
    df <- extinction_3epi_df_long %>% 
      filter(Pr == PR, Ps == PS) %>% 
      mutate(time = ifelse(time < 0.1, 0.1, round(time ,1))) %>% 
      group_by(time, episomes_per_cell) %>% 
      summarise(number_of_cells = mean(number_of_cells)) 
  }
  
  if(facet){
    if(Pr_facet == "100%"){
      df <- df %>% rbind(data.frame(Pr = factor("100%", levels = levels(extinction_3epi_df_long$Pr)), 
                                    Ps = factor("100%", levels = levels(extinction_3epi_df_long$Ps)), 
                                    time = 1:750, 
                                    episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                                    number_of_cells = 1000))
    }
    add_zeros <- function(data, id){
      if(id$Pr == "100%" & id$Ps == "100%") return(data)
      if(max(data$time) > 750) return(data)
      return(rbind(data, 
                   data.frame(time = seq(max(data$time) + 0.1, xlim + 10, by = 0.1), 
                              episomes_per_cell = factor("0", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                              number_of_cells = 1000)))
    } 
    
    add_time_0 <- function(data, id){
      return(rbind(data, 
                   data.frame(time = 0, 
                              episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                              number_of_cells = 1000)))
    } 
    
    df <- df %>% group_by(Pr, Ps) %>% group_modify(add_zeros) %>% group_modify(add_time_0)
  }else if(PR == "100%" & PS == "100%"){
    df <- df %>% rbind(data.frame(time = 1:750, episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), number_of_cells = 1000))
  } 
  
  plot <- df %>% 
    ggplot(aes(time, number_of_cells, fill = episomes_per_cell, color = episomes_per_cell)) + 
    geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + 
    scale_fill_manual(values = safe_colorblind_palette) + 
    scale_color_manual(values = safe_colorblind_palette) + 
    coord_cartesian(c(0,xlim)) +
    guides(color = "none") + 
    theme(panel.grid = element_blank())
  
  if(abundance){
    plot <- plot +
      scale_y_continuous(breaks = c(0, 0.5, 1))  +
      labs(x = "Time (generations)", fill = "Episomes per cell", y = "Relative\nAbundance") 
  }else{
    plot <- plot +
      scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1))  +
      labs(x = "Time (generations)", fill = "Episomes per cell", y = "Percent of\ncell population") 
  }
  
  if(facet) plot <- plot + facet_grid(Pr ~ Ps, labeller = "label_both", scale = "free_x")
  
  plot
}

new_plot <- trajectory_plot(Pr = "100%", Ps = "100%", 500) +
  
  trajectory_plot(Pr = "100%", Ps = "90%", 500) +
  
  distribution_plot1(Pr = "100%", Ps = "100%", 500) + 
  
  distribution_plot1(Pr = "100%", Ps = "90%", 500) +
  
  trajectory_plot(Pr = "90%", Ps = "100%", 100) +
  
  trajectory_plot(Pr = "90%", Ps = "90%", 100) +
  
  distribution_plot1(Pr = "90%", Ps = "100%", 100) +
  
  distribution_plot1(Pr = "90%", Ps = "90%", 100) +
  
  plot_layout(ncol = 2, byrow = T, heights = c(3,1, 3, 1), guides = "collect") &
  # theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", fill = "Episomes per cell") &
  guides(color = "none")

ggsave(here(out_folder, "new_simulations_fig.png"), new_plot, width = 10, height = 10)

poster_plot <- trajectory_plot(Pr = "100%", Ps = "100%", 500) +
  
  trajectory_plot(Pr = "80%", Ps = "90%", 50) +
  
  distribution_plot1(Pr = "100%", Ps = "100%", 500) + 
  
  distribution_plot1(Pr = "80%", Ps = "90%", 50) +
  
  trajectory_plot(Pr = "90%", Ps = "100%", 100) +
  
  trajectory_plot(Pr = "100%", Ps = "90%", 500) +
  
  distribution_plot1(Pr = "90%", Ps = "100%", 100) +
  
  distribution_plot1(Pr = "100%", Ps = "90%", 500) +
  
  plot_layout(ncol = 2, byrow = T, heights = c(3,1, 3, 1), guides = "collect") &
  # theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", fill = "Episomes\nper cell") &
  guides(color = "none") & 
  theme(plot.title = element_text(hjust = 0.5)) & 
  theme(
    # panel.background = element_rect(fill='black'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = "transparent") #transparent legend panel
  )

ggsave(here(out_folder, "poster_simulation_fig.png"), poster_plot , width = 11, height = 9, bg = "transparent")

trajectory_plot(Pr = "100%", Ps = "100%", 500) + 
  theme(
    # panel.background = element_rect(fill='white'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = "transparent") #transparent legend panel
  )



distribution_plot1_v2 <- function(Pr = "100%", Ps = "100%", xlim, facet = F, Pr_facet = "all"){
  PR <- Pr
  PS <- Ps
  
  # browser()
  
  if(facet){
    df <- extinction_3epi_df_long %>% 
      mutate(time = round(time), frac = number_of_cells/total) %>% 
      select(Pr, Ps, time, episomes_per_cell, number_of_cells = frac) 
    if(!any(Pr_facet == "all")) df <- df %>% filter(Pr %in% Pr_facet)
  }else{
    df <- extinction_3epi_df_long %>% 
      filter(Pr == PR, Ps == PS) %>% 
      mutate(time = round(time), frac = number_of_cells/total) %>% 
      select(Pr, Ps, time, episomes_per_cell, number_of_cells = frac) 
  }
  
  if(facet){
    if(Pr_facet == "100%"){
      df <- df %>% rbind(data.frame(Pr = factor("100%", levels = levels(extinction_3epi_df_long$Pr)), 
                                    Ps = factor("100%", levels = levels(extinction_3epi_df_long$Ps)), 
                                    time = 1:750, 
                                    episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                                    number_of_cells = 1000))
    }
    add_zeros <- function(data, id){
      if(id$Pr == "100%" & id$Ps == "100%") return(data)
      return(rbind(data, 
                   data.frame(time = (max(data$time) + 1):750, 
                              episomes_per_cell = factor("0", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                              number_of_cells = 1000)))
    } 
    
    df <- df %>% group_by(Pr, Ps) %>% group_modify(add_zeros)
  }else if(PR == "100%" & PS == "100%"){
    df <- df %>% rbind(data.frame(time = 1:750, episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), number_of_cells = 1000))
  } 
  
  # print(head(df))
  plot <- df %>%
    ggplot(aes(time, log10(number_of_cells + 1/10000), fill = episomes_per_cell, color = episomes_per_cell)) +
    stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
    stat_summary(fun.data="mean_cl_normal", geom="smooth") +
    # geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
    # scale_fill_manual(values = safe_colorblind_palette) +
    scale_color_manual(values = safe_colorblind_palette) +
    coord_cartesian(c(0,xlim)) +
    # scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1))  +
    # labs(x = "Time (generations)", fill = "Episomes per cell", y = "Percent of\ncell population") +
    labs(x = "Time (generations)", color = "Episomes per cell") +
    guides(color = "none") +
    theme(panel.grid = element_blank())
  
  if(facet) plot <- plot + facet_grid(Pr ~ Ps, labeller = "label_both")
  
  plot
}

supplemental_plot <- distribution_plot1(xlim = 700, facet = T, Pr_facet = "100%", abundance = T) + 
  theme(axis.title.x = element_blank()) +
  distribution_plot1(xlim = 250, facet = T, Pr_facet = "95%", abundance = T) + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  distribution_plot1(xlim = 120, facet = T, Pr_facet = "90%", abundance = T) +
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  distribution_plot1(xlim = 60, facet = T, Pr_facet = "80%", abundance = T) + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_vline(data = . %>% filter(Ps == "100%"),
             aes(xintercept = 5), lty = "dashed") +
  distribution_plot1(xlim = 20, facet = T, Pr_facet = "50%", abundance = T)  + 
  theme(strip.background.x = element_blank(), strip.text.x = element_blank()) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.justification = c(0.5, 1)) & 
  labs(fill = "Episomes\nper cell")

ggsave(here(out_folder, "supplemental_episome_distribution.png"), supplemental_plot , width = 12, height = 8)

supplemental_plot_example <- extinction_3epi_df_long %>% 
  filter(Ps == "100%", Pr == "80%") %>% 
  mutate(time = ifelse(time < 0.1, 0.1, round(time ,1))) %>% 
  group_by(Pr, Ps, time, episomes_per_cell) %>% 
  summarise(number_of_cells = mean(number_of_cells)) %>% 
  ungroup %>% 
  filter(time == 5) %>% 
  ggplot(aes("5th generation", number_of_cells, fill = episomes_per_cell)) + 
  geom_bar(stat = 'identity', position = position_fill(reverse = TRUE), show.legend = F) + 
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_text(data = . %>% filter(episomes_per_cell %in% c(0:3)), 
            aes(label = scales::percent(number_of_cells / sum(number_of_cells), accuracy = 0.1),
                y = (cumsum(number_of_cells) - 0.5 * number_of_cells) / sum(number_of_cells)),
            fontface = 2) + 
  scale_y_continuous(breaks = c(0, 0.5, 1))  +
  labs( y = "Relative Abundance") + 
  theme_minimal() +
  theme( axis.title.x = element_blank())


ggsave(here(out_folder, "supplemental_episome_distribution_example.png"), supplemental_plot_example , 
       width = 1.5, height = 3)  

#####
# Plot histograms with distribution of episomes in population at varying time points
#####
distribution_plot2 <- function(generation, PR, PS){
  df <- extinction_3epi_df_long %>% 
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
new_plot2 <- cowplot::plot_grid(
  trajectory_plot(Pr = "100%", Ps = "100%", 500),
  trajectory_plot(Pr = "100%", Ps = "90%", 500),
  plot_epi_distribution(1,1) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  plot_epi_distribution(1,0.9) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  x_label, x_label,
  trajectory_plot(Pr = "90%", Ps = "100%", 100),
  trajectory_plot(Pr = "90%", Ps = "90%", 100),
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

time_to_extinction_summary <- time_to_extinction %>% 
  mutate(Pr = paste0(pRep*100, "%"),
         Ps = paste0(pSeg*100, "%")) %>% 
  filter(!(pRep == 1 & pSeg == 1)) %>% 
  group_by(Pr, Ps, pRep, pSeg) %>% 
  summarise(sd = sd(ExtinctionTime), ExtinctionTime = mean(ExtinctionTime))

# sup1 <- time_to_extinction_summary %>% 
#   ggplot(aes(pSeg*100, ExtinctionTime, color = factor(pRep*100, levels= rev(sort(unique(pRep*100)))), group = pRep)) + 
#   geom_point(size = 2) + geom_line(lwd = 1) + 
#   geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
#   # add infity for Pr = Ps = 1
#   geom_text(data = data.frame(pSeg = 1, pRep = 1, ExtinctionTime = Inf), vjust = 1, aes(label = "+"), size = 7, show.legend = F) +
#   scale_color_manual(values = rev(sort(scico_palette_data("acton", T))[c(1, 21, 31, 41, 51, 71, 81)])) + 
#   labs(x = "Ps (%)", y = "Time to episome extinction (generations)", color = "Pr (%)") +
#   theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))
# 
# sup2 <- time_to_extinction_summary %>% 
#   mutate(pRep = ifelse(pSeg == 0.5, pRep - 0.015, 
#                        ifelse(pSeg == 0.75, pRep - 0.01,
#                               ifelse(pSeg == 0.8, pRep - 0.005, 
#                                      ifelse(pSeg == 0.85, pRep , 
#                                             ifelse(pSeg == 0.9, pRep + 0.005, 
#                                                    ifelse(pSeg == 0.95 , pRep + 0.01, 
#                                                           ifelse(pSeg == 1, pRep + 0.015, pRep)))))))) %>% 
#   ggplot(aes(pRep*100, ExtinctionTime, color = factor(pSeg*100, levels = rev(sort(unique(pSeg*100)))), group = pSeg)) + 
#   geom_line(lwd = 1) +
#   geom_point(size = 2) + 
#   geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
#   # add infity for Pr = Ps = 1
#   geom_text(data = data.frame(pSeg = 1, pRep = 1.01, ExtinctionTime = Inf), vjust = 1, aes(label = "+"), size = 7, show.legend = F) +
#   # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 41, 61, 81)]) + 
#   scale_color_manual(values = rev(sort(scico_palette_data("oslo", T))[c(1, 21, 31, 45, 51, 61, 81)]) )+ 
#   labs(x = "Pr (%)", y = "Time to episome extinction (generations)", color = "Ps (%)") +
#   theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))


sup1 <- time_to_extinction_summary %>% 
  # Format data
  ungroup %>%
  mutate(legend = factor(pRep*100, levels= rev(sort(unique(pRep*100))))) %>% 
  
  # Create plot
  ggplot( aes(x = pSeg * 100, y = ExtinctionTime, color = legend, group = pRep)) +
  geom_point(size = 2) +
  geom_line(lwd = 1) +
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) +
  
  # Add infinity marker
  geom_text(data = data.frame(pSeg = 1, pRep = 1, ExtinctionTime = Inf, legend = "100"), 
            aes(label = "+"), size = 7, show.legend = FALSE) +
  
  # Add dashed line to infinity
  geom_line(data = data_for_plot %>% 
              filter(pRep == 1, pSeg == 0.95) %>%
              add_row(pRep = 1, pSeg = 1, ExtinctionTime = Inf, legend = "100"), 
            linetype = "dotted", lwd = 1) +
  
  # Label lines
  ggrepel::geom_text_repel(data = data_for_plot %>% 
                             group_by(pRep) %>% 
                             filter(pSeg == max(pSeg)),
                           aes(x = 100, y = ifelse(pRep == 1, Inf, ExtinctionTime), label = paste("Pr:", Pr)), 
                           show.legend = FALSE, fontface = 2,
                           min.segment.length = 0, xlim = c(100, 115), direction = "y", nudge_x = 10) +
  
  # Format plot
  labs(x = "Ps (%)", y = bquote(T[E]~(generations)), color = "Pr (%)", shape = "Pr (%)", fill = "Pr (%)") +
  scale_color_manual(values = rev(sort(scico_palette_data("acton", T))[c(1, 21, 31, 41, 51, 71, 81)])) +
  coord_cartesian(clip = 'off', xlim = c(50, 100)) +
  theme(plot.margin = margin(0.1, 2.6, 0.1, 0.1, "cm"),
        legend.position = "none")

# Define color palette
Ps_colors <- rev(sort(scico_palette_data("oslo", T))[c(1, 21, 31, 45, 51, 61, 81)])

sup2 <- time_to_extinction_summary %>% 
  # Format data
  mutate(pRep = case_when(
    pSeg == 0.5  ~ pRep - 0.015,
    pSeg == 0.75 ~ pRep - 0.01,
    pSeg == 0.8  ~ pRep - 0.005,
    pSeg == 0.85 ~ pRep,
    pSeg == 0.9  ~ pRep + 0.005,
    pSeg == 0.95 ~ pRep + 0.01,
    pSeg == 1    ~ pRep + 0.015,
    TRUE         ~ pRep
  )) %>%
  ungroup() %>%
  mutate(legend = factor(pSeg * 100, levels = rev(sort(unique(pSeg * 100))))) %>% 
  
  # Create plot
  ggplot(aes(x = pRep * 100, y = ExtinctionTime, color = legend, fill = legend, group = pSeg)) +
  geom_line(lwd = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) +
  
  # Add infinity marker
  geom_text(data = data.frame(pSeg = 1, pRep = 1.015, ExtinctionTime = Inf, legend = "100"), 
            aes(label = "+"), size = 7, show.legend = FALSE) +
  
  # Add dashed line to infinity
  geom_line(data = sup2_data %>%
              filter(pSeg == 1, pRep == max(pRep)) %>%
              add_row(pRep = 1.015, pSeg = 1, ExtinctionTime = Inf, legend = "100"), 
            linetype = "dotted", lwd = 1) +
  
  # Label lines
  ggrepel::geom_text_repel(data = sup2_data %>%
                             group_by(pSeg) %>%
                             filter(pRep == max(pRep)),
                           aes(x = ifelse(pSeg == 1, 101.5, pRep * 100), 
                               y = ifelse(pSeg == 1, Inf, ExtinctionTime), 
                               label = paste("Ps:", Ps)), 
                           show.legend = FALSE, fontface = 2,
                           min.segment.length = 0, xlim = c(100, 115), direction = "y", nudge_x = 10) +
  
  # Format plot
  labs(x = "Pr (%)", y = bquote(T[E]~(generations)), color = "Ps (%)", shape = "Ps (%)", fill = "Ps (%)") +
  scale_color_manual(values = Ps_colors) +
  coord_cartesian(clip = 'off', xlim = c(48, 100)) +
  theme(plot.margin = margin(0.1, 2.6, 0.1, 0.1, "cm"),
        legend.position = "none")


annot_df <- time_to_extinction_summary %>%
  filter(Pr %in% unique(averages$Pr), 
         Ps %in% unique(averages$Ps)) %>% 
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))),
         Pr = factor(Pr, levels = levels(averages$Pr)))

sup3 <- averages %>%
  # Format Data
  group_by(Pr, Ps, time = round(time, 2)) %>%
  summarise(avg = mean(avg), .groups = 'drop') %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps)))) %>% 
  
  # Create plot
  ggplot( aes(x = time, y = avg, color = factor(Ps, levels = rev(sort(unique(Ps)))), group = Ps)) +
  geom_line(size = 1) +
  facet_wrap(~Pr, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = Ps_colors) +
  
  # Annotations
  geom_point(data = annot_df, aes(x = ExtinctionTime, y = 0), shape = 4, show.legend = FALSE) +
  geom_rect(data = sup3_data %>% filter(Pr != "100%"), aes(xmin = 690, xmax = Inf, ymin = 0.5, ymax = Inf), fill = "white", color = NA) +
  geom_rect(data = sup3_data %>% filter(Pr == "100%"), aes(xmin = 690, xmax = Inf, ymin = 0.5, ymax = 2.55), fill = "white", color = NA) +
  
  # Time to extinction text annotations
  geom_text(data = annot_df %>%
              mutate(y = case_when(
                Ps == 100 ~ 3,
                Ps == 95 ~ 2.5,
                Ps == 90 ~ 2,
                Ps == 80 ~ 1.5,
                TRUE ~ 1
              )),
            aes(x = 700, y = y, label = paste("T[E]:", round(ExtinctionTime), "%+-%", round(sd))),
            parse = TRUE, hjust = 0, vjust = 1, show.legend = FALSE, fontface = "bold") +
  
  # Labels and theme
  labs(x = "Time (generations)", y = "Average number of episomes per cell", color = "Ps (%)") 

sup <- ((plot_spacer() / sup1 / plot_spacer() / sup2 / plot_spacer() + plot_layout(heights = c(0.5,4,0.5,4,0.5))) | 
          plot_spacer() | 
          sup3) + plot_layout(widths = c(1,0.1,1))

ggsave(here(out_folder, "supplemental_figure.png"), sup , width = 11.5, height = 10)


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
  # add time to extinction annotation
  geom_point(data = time_to_extinction_summary %>%
               mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                      Pr = factor(Pr, levels = levels(averages$Pr))), 
             aes(ExtinctionTime, 0), shape = 4, show.legend = F) +
  geom_rect(aes(xmin = 1850, xmax = 2250, ymin = 0.5, ymax = 3.2), fill = "white", alpha = 0.05, color = NA) +
  geom_text(data = time_to_extinction_summary %>%
              mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                     Pr = factor(Pr, levels = levels(averages$Pr)),
                     y = ifelse(Ps == 100, 3, ifelse(Ps == 95, 2.5, ifelse(Ps == 90, 2, ifelse(Ps == 75, 1.5, 1))))), 
            aes(1850, y, label = paste("T[E]:", round(ExtinctionTime), "%+-%", round(sd))), 
            parse = T, hjust = 0, vjust = 1, 
            show.legend = F)  + 
  theme(legend.position = "bottom")

sup.2 <-( (sup1 + plot_spacer() + sup2 + plot_layout(widths = c(1, 0.1, 1), nrow = 1)) / 
            plot_spacer() / sup3.2) + plot_layout(heights = c(1,0.1,1), ncol = 1)

ggsave(here(out_folder, "supplemental_figure_v2_temp.png"), sup.2 , width = 15, height = 7)


## Calculate shannon index
shannon <- extinction_3epi_df_long %>% 
  filter(number_of_cells != 0) %>% 
  mutate(p = number_of_cells/total_cells) %>% 
  group_by(pRep, pSeg, trial, time) %>% 
  summarise(H = -sum(p*log(p)))

shannon %>% 
  group_by(pRep, pSeg, trial) %>% 
  summarise(H1 = max(H), H2 = mean(H)) %>% 
  mutate(Ps = 100*pSeg, Pr = 100*pRep) %>% 
  ggplot(aes(Pr, H2, color = as.factor(Ps))) + 
  geom_point() +
  # geom_linerange(data = . %>% group_by(Pr, Ps) %>% summarise(ymin = min(H), ymax = max(H)), 
  #                aes(ymin = ymin, ymax = ymax), linewidth = 1) +
  # facet_wrap(~Ps, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) 
# labs(x = "Segregation Efficiency (%)", y = "Maximum shannon index", color = "Ps (%)") 


shannon_plot <- shannon %>% 
  mutate(Pr = pRep*100, Ps = pSeg*100) %>% 
  filter(time <= 200) %>% 
  ggplot(aes(round(time,1), H, color = factor(Ps), group  = interaction(Ps, Pr))) +   
  # stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", aes(ymin = after_stat(ymin), ymax = after_stat(ymax)), alpha = 0.2, color = NA) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth") + 
  facet_wrap(~Pr, scales = 'free_x', nrow = 1, labeller = "label_both") + 
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31,  45, 61, 81)])  + 
  labs(x = "Time (generations)", y = "Shannon Index", color = "Ps (%)") + 
  ylim(c(0, max(shannon$H)))

ggsave(here(out_folder, "shannon_index.png"), shannon_plot, width = 10, height = 5)

peak_shannon_plot <- shannon %>% 
  mutate(Pr = pRep*100, Ps = pSeg*100) %>% 
  group_by(Pr, Ps, trial) %>% 
  summarise(H = max(H)) %>% 
  ggplot(aes(factor(Ps), H, color = factor(Ps), group  = interaction(Ps, Pr))) +   
  geom_boxplot() +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31,  45, 61, 81)])  + 
  # labs(x = "Replication Efficiency (%)", y = "Shannon Index", color = "Ps (%)") + 
  ylim(c(0, max(shannon$H))) + 
  facet_wrap(~Pr)

