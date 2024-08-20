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
rerun <- TRUE


out_folder <- here("results","simulations_constant_selection")

if(!file.exists(out_folder)) dir.create(out_folder, showWarnings = F)

if(!rerun){
  load(here(out_folder, "constant_pop_3epi.RData"))
  
  extinction_3epi_df <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 
  
  extinction_3epi_df_long <- pivot_extinction(extinction_3epi_df) %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
}else{
  
  ## set up a grid of parameters to simulate with and collect all results:
  pReps = 0.8
  pSegs = c(0.8, 0.9)
  # First run showed no extinction with n = 50 or higher
  n_cells = seq(5,40, by = 5)
  param_grid <- expand_grid(pRep = pReps, pSeg = pSegs, n_cells)
  
  extinction_3epi <- param_grid %>%
    pmap(extinction, nTrials = ntrials, n_epi = n_epi, selectAgainstZero= T)
  
  extinction_3epi_df <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 
  
  extinction_3epi_df_long <- pivot_extinction(extinction_3epi_df) %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
  
  # Plot average number of episomes per cell over time
  averages <- extinction_3epi_df_long %>% 
    mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) %>% 
    group_by(time = round(time,1), Pr, Ps, trial, n_cells) %>% 
    summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
    ungroup()
  
  save(extinction_3epi, param_grid, averages, 
       file = here(out_folder, "constant_pop_3epi_increasing_ncells_lower_values.RData"))
  
  
  # Plot distribution of episomes per cell over time
  p1 <- extinction_3epi_df_long %>% 
    ggplot(aes(time, frac, color = episomes_per_cell, group = interaction(trial, episomes_per_cell))) + 
    geom_line(alpha = 0.25) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    ggh4x::facet_grid2(n_cells ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
    labs(x = "time (generations)", y = "Fraction of\npopulation", color = "episomes per cell")
  
  ggsave(here(out_folder, "extinction_3epi_lines.png"), p1, width = 8, height = 7)
  
  
  p2.1 <- averages %>% 
    ggplot(aes(time, avg, group = trial)) + 
    geom_line(alpha = 0.25, color = "gray") +
    geom_point(data = . %>% group_by(trial, Pr, Ps, n_cells) %>% filter(time == max(time)), aes(y= 0)) +
    ggh4x::facet_grid2(Ps~n_cells, labeller = "label_both", scales = "free_x", independent = "x") + 
    labs(x = "time (generations)", y = "average number of episomes per cell") 
  
  ggsave(here(out_folder, "extinction_3epi_average_epi_freex.png"), p2.1, width = 12, height = 4)
  
  # p2.2 <- averages %>% 
  #   ggplot(aes(time, avg, group = trial)) + 
  #   geom_line(alpha = 0.25, color = "gray") +
  #   facet_grid(Pr ~ Ps, labeller = "label_both") +
  #   scale_x_log10() +
  #   labs(x = "time (generations)", y = "average number of episomes per cell") +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 0.8))
  # 
  # ggsave(here(out_folder, "extinction_3epi_average_epi_logscale.png"), p2.2, width = 9, height = 7)
  
  # Plot histograms of time to extinction
  time_to_extinction <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$ExtinctionTime))
  
  p3.2 <- time_to_extinction %>% 
    rename(Pr = pRep, Ps = pSeg) %>% 
    # filter(Ps*Pr != 1) %>%
    # arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) %>% 
    ggplot(aes(ExtinctionTime)) + 
    geom_histogram(bins = 30, aes(y = after_stat(density))) + 
    ggh4x::facet_grid2(Ps~n_cells, labeller = "label_both", scales = "free_y", independent = "y") +
    labs(x = "time to episome extinction (generations)")
  
  ggsave(here(out_folder, "time_to_exinction_3epi_freex.png"), p3.2, width = 12, height = 4)
  
  ## Fraction that go extinct:
  fraction_extinct <- time_to_extinction %>% 
    group_by(n_cells, pRep, pSeg = pSeg*100) %>% 
    summarise(frac_extinct = 1-sum(ExtinctionTime >= 700)/n()) %>% 
    ggplot(aes(n_cells, frac_extinct, color = as.factor(pSeg), lty = as.factor(pSeg))) + 
    geom_line(size = 1) + 
    geom_point(show.legend = F, size = 2) +
    labs(x = "Number of cells in population", y = "Fraction of 100 simulated populations\nwith episomal extinction", 
         color = "Segregation\nEfficiency (%)", linetype = "Segregation\nEfficiency (%)",
         title = "Fraction of episomal extinction in simulations with selection",
         subtitle = "80% replication efficieincy, varying segregation efficiency and population size") + 
    theme(legend.position = c(1,1), legend.justification = c(1.1,1.1), legend.key.size = unit(1, "cm")) + 
    scale_color_brewer(palette = "Set1") + 
    scale_linetype_manual(values = c(2,3))
  
  ggsave(here(out_folder, "Fraction_of_extinct_populations.png"), fraction_extinct, width = 6, height = 4)
  
  time_extinct_dist <- time_to_extinction %>% 
    ggplot(aes(as.factor(n_cells), ExtinctionTime, color = as.factor(pSeg*100))) + 
    geom_boxplot(outlier.shape = 1) +
    geom_point(data = . %>% filter(ExtinctionTime < 700, n_cells < 25) %>% 
                 group_by(pRep, pSeg, n_cells) %>% summarise(ExtinctionTime = mean(ExtinctionTime)), 
               shape = 24, size = 1.5, position = position_dodge(0.75),
               aes(fill = as.factor(pSeg)), color = "black", show.legend = F) +
    scale_color_brewer(palette = "Set1") + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Number of cells in population", y = "Number of divisions prior to episomal extinction", 
         color = "Segregation\nEfficiency (%)", 
         title = "Time until episomal extinction in simulations with selection",
         subtitle = "80% replication efficieincy, varying segregation efficiency and population size") + 
    theme(legend.position = c(1,0), legend.justification = c(1.1,-0.1), legend.key.size = unit(1, "cm")) 
    
  ggsave(here(out_folder, "Time_till_extinction_distributions.png"), time_extinct_dist, width = 6, height = 4)
  
}


time_to_extinction %>% 
  mutate(trial = rep(1:100, nrow(param_grid))) %>% 
  filter(n_cells == 25, pSeg == 0.8) %>% 
  filter(ExtinctionTime < 700)


## Plot average episome over time for case with moderate defect
simple_plot <-  averages %>% 
  filter(Pr == "100%", Ps == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") +
  labs(title = "100% Replication, 100% Segregation") +
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
    df <- extinction_3epi_df_long %>% 
      group_by(Pr, Ps, time = round(time), episomes_per_cell) %>% 
      summarise(number_of_cells = mean(number_of_cells))  
    if(!any(Pr_facet == "all")) df <- df %>% filter(Pr %in% Pr_facet)
  }else{
    df <- extinction_3epi_df_long %>% 
      filter(Pr == PR, Ps == PS) %>% 
      group_by(time = round(time), episomes_per_cell) %>% 
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
      return(rbind(data, 
            data.frame(time = (max(data$time) + 1):750, 
                       episomes_per_cell = factor("0", levels = levels(extinction_3epi_df_long$episomes_per_cell)), 
                       number_of_cells = 1000)))
    } 
    
    df <- df %>% group_by(Pr, Ps) %>% group_modify(add_zeros)
  }else if(PR == "100%" & PS == "100%"){
    df <- df %>% rbind(data.frame(time = 1:750, episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), number_of_cells = 1000))
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

time_to_extinction_summary <- time_to_extinction %>% 
  mutate(Pr = paste0(pRep*100, "%"),
         Ps = paste0(pSeg*100, "%")) %>% 
  filter(!(pRep == 1 & pSeg == 1)) %>% 
  group_by(n_cells, Pr, Ps, pRep, pSeg) %>% 
  summarise(min = min(ExtinctionTime), sd = sd(ExtinctionTime), ExtinctionTime = mean(ExtinctionTime))

sup1 <- time_to_extinction_summary %>% 
  ggplot(aes(pSeg*100, ExtinctionTime, color = as.factor(n_cells), group = n_cells)) + 
  geom_point(size = 2) + geom_line(lwd = 1) + 
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
  # add infity for Pr = Ps = 1
  # scale_color_manual(values = sort(scico_palette_data("acton", T))[c(1, 31, 51, 71, 81)]) + 
  labs(x = "Ps (%)", y = "Time to episome extinction (generations)", color = "Number of cells") +
  theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))

sup2 <- time_to_extinction_summary %>% 
  mutate(n_cells = ifelse(pSeg == 0.5, n_cells - 0.01, 
                       ifelse(pSeg == 0.75, n_cells - 0.005, 
                              ifelse(pSeg == 0.95, n_cells + 0.005, 
                                     ifelse(pSeg == 1, n_cells + 0.01, n_cells))))) %>% 
  ggplot(aes(n_cells, ExtinctionTime, color = as.factor(pSeg*100), group = pSeg)) + 
  geom_line(lwd = 1) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
  # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 41, 61, 81)]) + 
  # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) + 
  labs(x = "Number of cells", y = "Time to episome extinction (generations)", color = "Ps (%)") +
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

ggsave(here(out_folder, "supplemental_figure_v2.png"), sup.2 , width = 15, height = 7)


## Calculate shannon index
shannon <- extinction_3epi_df_long %>% 
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



extinction_3epi_df_long %>% filter(pSeg == 0.9, pRep == 0.9) %>% 
  filter(episomes_per_cell == 0) %>% 
  mutate(latent_cells = total_cells - number_of_cells) %>% 
  ggplot(aes(time, latent_cells)) + 
  geom_line(aes(group = trial), alpha = 0.1, color = "gray")  + 
  # geom_smooth(formula = (log(y)~x), method = "lm")
  geom_smooth()
  # scale_y_log10()
  # facet_grid(Pr ~ Ps, labeller = "label_both")

