#####
# Simulating a constant-sized cell population without selection
# Generate Figure 7, S8, S9, S10
### Details of simulations:
# - 100 populations are stochastically simulated (ntrials)
# - populations start with 1,000 cells, each with 3 episomes (n_epi)
# - replication and segregation efficiencies range from 50% to 100%
# - simulations end when episomes are extinct from the population
### Outputs:
# - Number of cells with k episomes at each time point for k = 0,...,9
# - Cells with >9 episomes are grouped at 9 episomes
# - Average number of episomes per cell in the population at each time point
# - Time (in generations) until episomal extinction from the cell population (time to extinction)
#####

### Setup ###################################################################### 

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
theme_set(theme_bw())

source(here("scripts", "functions_simulations.R"))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


### Script inputs ##############################################################
# Run 100 trials
ntrials <- 100
# Start with three episomes per cell
n_epi <- 3 

# Indicator if simulations should be rerun (TRUE) or loaded from previous run (FALSE)
rerun <- FALSE

# Create output folder
out_folder <- here("results","simulations_constant")
if(!file.exists(out_folder)) dir.create(out_folder, showWarnings = F)

## Set up a grid of parameters to simulate with and collect all results:
pReps = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
pSegs = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
param_grid <- expand_grid(pRep = pReps, pSeg = pSegs)

### Run simulations ############################################################

if(!rerun){ 
  # Load results from previous run if they exist
  load(here(out_folder, "constant_pop_3epi.RData"))
}else{ 
  # Run simulations
  extinction_3epi <- param_grid %>% 
    pmap(extinction, nTrials = ntrials, n_epi = n_epi)
  
  save(extinction_3epi,
       file = here(out_folder, "constant_pop_3epi.RData"))
}

### Format outputs and analyze results #########################################

# Format results of simulations
extinction_3epi_df <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 

extinction_3epi_df_long <- extinction_3epi_df %>% 
  # Thin results for memory issues
  filter(pRep %in% c(0.5, 0.8, 0.9, 0.95, 1), 
         pSeg %in% c(0.5, 0.8, 0.9, 0.95, 1)) %>% 
  # Format data
  pivot_extinction() %>% 
  # Calculate fraction of population with k cells
  group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
  mutate(total_cells = sum(number_of_cells)) %>% 
  ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
  arrange(desc(Pr), Ps) %>% 
  mutate(Pr = paste0(Pr*100, "%"),
         Ps = paste0(Ps*100, "%"),
         Pr = fct_inorder(factor(Pr)),
         Ps = fct_inorder(factor(Ps))) 

# Calculate average number of episomes per cell over time
averages_simulated <- extinction_3epi_df_long %>% 
  mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) %>% 
  # Bin times to the nearest tenth
  group_by(time = round(time,1), Pr, Ps, trial) %>% 
  # Calculate average number of episomes per cell
  summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
  ungroup()

# fill out Pr = Ps = 100% to time of longest simulation
start_time <- max(averages_simulated %>% filter(Pr == "100%", Ps == "100%") %>% pull(time)) + 1
averages <- rbind(averages_simulated, 
                  data.frame(Pr = factor("100%", levels = levels(averages_simulated$Pr)), 
                             Ps = factor("100%", levels = levels(averages_simulated$Ps)), 
                             avg = 3, expand.grid(time = start_time:(max(averages_simulated$time) + 30), trial = 1:100))
)

# Extract time to extinction 
time_to_extinction <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$ExtinctionTime))

### Plot results ###############################################################

## Plot distribution of episomes per cell over time
p1 <- extinction_3epi_df_long %>% 
  ggplot(aes(time, frac, color = episomes_per_cell, group = interaction(trial, episomes_per_cell))) + 
  geom_line(alpha = 0.25) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
  labs(x = "time (generations)", y = "Fraction of\npopulation", color = "episomes per cell")

ggsave(here(out_folder, "extinction_3epi_lines.png"), p1, width = 8, height = 7)



## Function to plot trajectories of average number of episomes per cell over time
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

## Function to plot the distribution of episomes per cell over time
# Plotted as stacked barplot
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

# Create plot for MIDAS poster
poster_plot <- trajectory_plot(Pr = "100%", Ps = "100%", 500) +
  
  trajectory_plot(Pr = "80%", Ps = "90%", 50) +
  
  distribution_plot1(Pr = "100%", Ps = "100%", 500) + 
  
  distribution_plot1(Pr = "80%", Ps = "90%", 50) +
  
  trajectory_plot(Pr = "90%", Ps = "100%", 100) +
  
  trajectory_plot(Pr = "100%", Ps = "90%", 500) +
  
  distribution_plot1(Pr = "90%", Ps = "100%", 100) +
  
  distribution_plot1(Pr = "100%", Ps = "90%", 500) +
  
  plot_layout(ncol = 2, byrow = T, heights = c(3,1, 3, 1), guides = "collect") 
  labs(x = "Time (generations)", fill = "Episomes\nper cell") &
  guides(color = "none") & 
  theme(plot.title = element_text(hjust = 0.5)) & 
  theme(
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = "transparent") #transparent legend panel
  )

ggsave(here(out_folder, "poster_simulation_fig.png"), poster_plot , width = 11, height = 9, bg = "transparent")

## Function to plot the distribution of episomes per cell over time
# Alternative from above, plotted as histogram
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

# Plot histograms at specified time points
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
  
  plot
}

# Create figure for paper (Figure 7)
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

ggsave(here(out_folder, "simulations_fig.png"), new_plot2, width = 10, height = 10, bg = "white")



### Supplemental figures of simulation summary #################################

## Figure S8
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

## Example for Figure S8
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

## Plot average number of episomes per cell over time  (Figure S9)
# Use a different x-axis range for each Pr value
p2 <- averages %>% 
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


ggsave(here(out_folder, "extinction_3epi_average_epi_xaxis_by_Pr.png"), p2, width = 8, height = 7)

## Time to extinction plots (Figure S10)
time_to_extinction_summary <- time_to_extinction %>% 
  mutate(Pr = paste0(pRep*100, "%"),
         Ps = paste0(pSeg*100, "%")) %>% 
  filter(!(pRep == 1 & pSeg == 1)) %>% 
  group_by(Pr, Ps, pRep, pSeg) %>% 
  summarise(sd = sd(ExtinctionTime), ExtinctionTime = mean(ExtinctionTime))

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
  geom_line(data = . %>%
              select(pRep, pSeg, ExtinctionTime, legend) %>%
              filter(pRep == 1, pSeg == 0.95) %>%
              add_row(pRep = 1, pSeg = 1, ExtinctionTime = Inf, legend = "100"),
            linetype = "dotted", lwd = 1) +
  
  # Label lines
  ggrepel::geom_text_repel(data = . %>%
                             group_by(pRep) %>%
                             filter(pSeg == max(pSeg)),
                           # aes(x = 100, y = ifelse(pRep == 1, Inf, ExtinctionTime), label = paste("Pr:", Pr)),
                           aes(x = ifelse(pRep == 1, 95, 100), y =  ExtinctionTime, label = paste("Pr:", Pr)),
                           show.legend = FALSE, fontface = 2,
                           min.segment.length = 0, xlim = c(100, 115),
                           direction = "y", nudge_x = 10) +
  
  
  # Format plot
  labs(x = "Ps (%)", y = bquote(T[E]~(generations)), color = "Pr (%)", shape = "Pr (%)", fill = "Pr (%)") +
  scale_color_manual(values = rev(sort(scico_palette_data("acton", T))[c(1, 21, 31, 41, 51, 71, 81)])) +
  scale_x_continuous(breaks = seq(50, 100, by = 10), limits = c(50, 110)) +
  coord_cartesian(clip = 'off') +
  theme_minimal() +
  theme(#plot.margin = margin(0.1, 2.6, 0.1, 0.1, "cm"),
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
  geom_line(data = . %>%
              select(pRep, pSeg, ExtinctionTime, legend) %>% 
              filter(pSeg == 1) %>%
              filter(pRep == max(pRep)) %>% 
              add_row(pRep = 1.015, pSeg = 1, ExtinctionTime = Inf, legend = "100"), 
            linetype = "dotted", lwd = 1) +
  
  # Label lines
  ggrepel::geom_text_repel(data = . %>%
                             group_by(pSeg) %>%
                             filter(pRep == max(pRep)),
                           aes(x = ifelse(pSeg == 1, 101.5, pRep * 100), 
                               y = ifelse(pSeg == 1, 500, ExtinctionTime), 
                               label = paste("Ps:", Ps)), 
                           show.legend = FALSE, fontface = 2,
                           min.segment.length = 0, xlim = c(100, 115), direction = "y", nudge_x = 10) +
  
  # Format plot
  labs(x = "Pr (%)", y = bquote(T[E]~(generations)), color = "Ps (%)", shape = "Ps (%)", fill = "Ps (%)") +
  scale_color_manual(values = Ps_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(50, 100, by = 10)) +
  coord_cartesian(clip = 'off', xlim = c(48, 110)) +
  theme(#plot.margin = margin(0.1, 2.6, 0.1, 0.1, "cm"),
        legend.position = "none")


annot_df <- time_to_extinction_summary %>%
  filter(Pr %in% unique(averages$Pr), 
         Ps %in% unique(averages$Ps)) %>% 
  mutate(Ps = factor(Ps, levels = rev(levels(averages$Ps))),
         Pr = factor(Pr, levels = levels(averages$Pr))) %>% 
  left_join(
    averages %>% group_by(Pr) %>% summarise(xmin = max(time)) %>% select(Pr, xmin)
  )
  


sup3 <- averages %>%
  # Format Data
  group_by(Pr, Ps, time) %>%
  # filter(!(Pr == "100%" & Ps == "100%")) %>% 
  summarise(avg = sum(avg)/100, .groups = 'drop') %>%
  # summarise(avg = mean(avg), .groups = 'drop') %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps)))) %>% 
  
  # Create plot
  ggplot( aes(x = time, y = avg, color = factor(Ps, levels = rev(sort(unique(Ps)))), group = Ps)) +
  geom_line(size = 1) +
  facet_wrap(~Pr, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = Ps_colors) +
  
  # Annotations
  geom_point(data = annot_df, aes(x = ExtinctionTime, y = 0), shape = 4, show.legend = FALSE) +
  geom_rect(data = . %>% filter(Pr != "100%"), aes(xmin = 690, xmax = Inf, ymin = 0.5, ymax = Inf), fill = "white", color = NA) +
  geom_rect(data = . %>% filter(Pr == "100%"), aes(xmin = 690, xmax = Inf, ymin = 0.5, ymax = 2.55), fill = "white", color = NA) +
  
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

ggsave(here(out_folder, "supplemental_figure_v1.png"), sup , width = 11.5, height = 10)


# Version of supplemental figure with kaplan-meier-type curve
sup3.2 <- averages %>% 
  # calculate the number of populations that still have episomes at each time point
  group_by(Pr, Ps, time) %>% count() %>% 
  # calculate the proportion of populations at each time point 
  # note that at early time points, there may not be a record for every population 
  # due to the size of the simulation time step, making it seem like there are 
  # fewer than 100 populations that still have episomes. Therefore, we calculate 
  # these timepoints differently so that they will be 100%
  mutate(p = ifelse(time < 1 & n < 100, n/n*100, n),
         Ps = fct_rev(Ps)) %>% 
  ggplot(aes(time, p, color = Ps)) + 
  geom_line(linewidth = 1) + 
  
  # Annotations
  geom_point(data = annot_df, aes(x = ExtinctionTime, y = 0), size = 2, show.legend = FALSE) +
  geom_rect(data = . %>% group_by(Pr) %>% summarise(xmin = 0.78*max(time)), 
            aes(xmin = xmin, xmax = Inf, ymin = 23, ymax = 95), fill = "white", color = NA, inherit.aes = F) +

  # Time to extinction text annotations
  geom_text(data = annot_df %>%
              mutate(y = case_when(
                Ps == "100%" ~ 95,
                Ps == "95%" ~ 80,
                Ps == "90%" ~ 65,
                Ps == "80%" ~ 50,
                TRUE ~ 35
              )),
            aes(x = 0.8*xmin, y = y, label = paste("T[E]:", round(ExtinctionTime), "%+-%", round(sd))),
            parse = TRUE, hjust = 0, vjust = 1, show.legend = FALSE, fontface = "bold") +
  
  # Labels and theme
  facet_wrap(~Pr, scales = "free_x", ncol = 1, labeller = "label_both") + 
  scale_color_manual(values = Ps_colors) +
  labs(x = "Time (generations)", y = "Simulations with episomes remaining (%)", color = "Segregation\nEfficiency (Ps)") 
  

sup <- ((plot_spacer() / sup1 / plot_spacer() / sup2 / plot_spacer() + plot_layout(heights = c(0.5,4,0.5,4,0.5))) | 
          plot_spacer() |
          sup3.2) + plot_layout(widths = c(1,0.1,1))

ggsave(here(out_folder, "supplemental_figure.png"), sup , width = 11.5, height = 10)



