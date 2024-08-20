# This script creates supplemental figures for all experimental conditions
#####
# Setup
#####
library(here)
library(tidyverse)
library(cluster)
library(factoextra)
library(scico)

source(here("scripts", "functions_inference.R"))
source(here("scripts", "functions_run_pipeline.R"))

theme_set(theme_bw())

#####
## Load in results

# results_folder_suffix <- "updated"
results_folder_suffix <- "updated_pdf_n_prior"
supplemental_folder <- paste0("supplemental_figures_", results_folder_suffix)

# Figure 2
results_folder <- here("results", str_interp("Figure_2_${results_folder_suffix}"))
daughter_cell_file <- "Figure 2 dots new version.xlsx"
mother_cell_file <- "Figure 2 non dividing cells 2.xlsx"
Figure2_data <- load_data(mother_cell_file, daughter_cell_file)
daughter_cell_data <- Figure2_data$daughter_cell_data
mother_cell_data <- Figure2_data$mother_cell_data
Figure2_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

# Figure 3
results_folder <- here("results", str_interp("Figure_3_${results_folder_suffix}"))
daughter_cell_file <- "Fig3 dividing cells.xlsx"
mother_cell_file <- "Fig3 Non dividing cells.xlsx"
Figure3_data <- load_data(mother_cell_file, daughter_cell_file)
daughter_cell_data <- Figure3_data$daughter_cell_data
mother_cell_data <- Figure3_data$mother_cell_data
Figure3_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

# Figure 5
results_folder <- here("results", str_interp("Figure_5_${results_folder_suffix}"))
daughter_cell_file <- "Fig5 dividing cells.xlsx"
mother_cell_file <- "Fig5 Non dividing cells.xlsx"
Figure5_data <- load_data(mother_cell_file, daughter_cell_file)
daughter_cell_data <- Figure5_data$daughter_cell_data
mother_cell_data <- Figure5_data$mother_cell_data
Figure5_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

# Figure 6
if(results_folder_suffix == "updated"){
  results_folder <- here("results", str_interp("Figure_6_different_mean"))
}else{
  results_folder <- here("results", str_interp("Figure_6_${results_folder_suffix}"))  
}
daughter_cell_file <- "Fig6 dividing cells.xlsx"
mother_cell_file <- "Fig6 Non dividing cells.xlsx"
Figure6_data <- load_data(mother_cell_file, daughter_cell_file)
daughter_cell_data <- Figure6_data$daughter_cell_data
mother_cell_data <- Figure6_data$mother_cell_data
Figure6_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, same_mu = F)

#####


MCMC_vs_data <-  function(pipeline_output, data, results_folder, dataset,  same_mu = T, normalize = T){
  
  all_chains <- pipeline_output$all_chains
  daughter_cell_samples <- pipeline_output$daughter_cell_samples
  mother_cell_samples <- pipeline_output$mother_cell_samples
  lambda <- pipeline_output$X0_lambda
  convergence <- pipeline_output$convergence_metrics
  MLE_grid <- pipeline_output$MLE_grid
  
  mother_cell_data <- data$mother_cell_data
  daughter_cell_data <- data$daughter_cell_data
  
  intensity_data <- rbind(daughter_cell_data %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "daughter"), 
                          mother_cell_data  %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "mother"))
  
  ## Color-blind friendly color palette
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
  # How do these inferences relate to the observed intensity data?
  cluster_modes <- all_chains %>% 
    filter(chain == "chain1") %>% 
    select(grep("[0-9]", names(all_chains))) %>% 
    pivot_longer(everything(), names_to = "cluster_id", values_to = "n_episomes") %>% 
    group_by(cluster_id) %>% 
    mutate(sd = sd(n_episomes)) %>% 
    count(cluster_id, n_episomes, sd) %>% 
    group_by(cluster_id) %>%
    mutate(p = n/sum(n)) %>% 
    filter(n == max(n)) %>% 
    ungroup %>% 
    select(-n) %>% 
    rename(mode = n_episomes) %>% 
    mutate(set = ifelse(grepl("d", cluster_id), "daughter", "mother"))
  
  # browser()
  if(!normalize){
    dat <- intensity_data %>% 
      merge(cluster_modes) %>% 
      mutate(set = str_to_title(paste(set, "cell")),
             min_episome_in_cluster = factor(min_episome_in_cluster, levels = 1:4))
  }else{
    if(!same_mu){
      daughter_mean_intensity <- DescTools::Mode(round(all_chains$mu_d,1))
      mother_mean_intensity <- DescTools::Mode(round(all_chains$mu_m,1))
      
      dat <- intensity_data %>% 
        merge(cluster_modes) %>% 
        mutate(total_cluster_intensity = ifelse(set == "daughter", total_cluster_intensity/daughter_mean_intensity,
                                                total_cluster_intensity/mother_mean_intensity),
               set = str_to_title(paste(set, "cell")),
               min_episome_in_cluster = factor(min_episome_in_cluster, levels = 1:4))
    }else{
      mean_intensity <- DescTools::Mode(round(all_chains$mu,1))
      
      dat <- intensity_data %>% 
        merge(cluster_modes) %>%
        mutate(total_cluster_intensity = total_cluster_intensity/mean_intensity,
               set = str_to_title(paste(set, "cell")),
               min_episome_in_cluster = factor(min_episome_in_cluster, levels = 1:4))
    }
    
  }
  
  x_label <- ifelse(normalize, 
                    "Total intensity of LANA dot normalized by mean intensity of a single episome",
                    "Total intensity of LANA dot")
  
  cat("max sd: ", max(dat$sd))
  
  MCMC_vs_intensity_data <- dat %>% 
    ggplot(aes(total_cluster_intensity, as.factor(mode))) + 
    # geom_jitter(aes(shape = set, color = sd)) +
    geom_jitter(aes(shape = set, color = p)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(y = "Mode of the posterior distribution", 
         # color = "Standard deviation\nof the\nposterior distribution",
         color = "Posterior\nProbability\nof the Mode",
         shape = "Cell set",
         x = x_label) + 
    theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1)) + 
    scale_shape_manual(values = c(2,19)) + 
    # scale_color_viridis_c(limits = c(0, 1.15)) + 
    scale_color_viridis_c(limits = c(0, 1), direction = -1) + 
    guides(shape = "none")
  
  if(normalize) MCMC_vs_intensity_data <- MCMC_vs_intensity_data + scale_x_continuous(breaks = 0:20)
  
  # browser()
  intensity_data_plot <- intensity_data %>% 
    mutate(set = str_to_title(paste(set, "cell"))) %>% 
    ggplot(aes(total_cluster_intensity, shape = set, y = set)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() + 
    scale_shape_manual(values = c(2,19)) + 
    labs(x = "Total intensity of LANA dot",
         title = dataset) +
    theme(legend.position = "none", axis.title.y = element_blank()) + 
    scale_y_discrete(labels = c("Daughter Cells", "Non-dividing Cells"))
  
  if(!same_mu){
    intensity_data_plot <- intensity_data_plot + 
      geom_point(data = data.frame(set = c("Daughter Cell", "Mother Cell"), mean = c(daughter_mean_intensity, mother_mean_intensity)),
                 aes(mean, set), color = "red", size = 2, show.legend = F) + 
      geom_vline(xintercept = daughter_mean_intensity, color = "red", alpha = 0.5) + 
      geom_vline(xintercept = mother_mean_intensity, color = "red", alpha = 0.5)
  }else{
    intensity_data_plot <- intensity_data_plot + 
      geom_point(data = . %>% group_by(set) %>% 
                   summarise(mean = mean_intensity),
                 aes(mean, set), color = "red", size = 2, show.legend = F) + 
      geom_vline(xintercept = mean_intensity, color = "red", alpha = 0.5)
  }
  
  # mode_plot <- cluster_modes %>%
  #   mutate(set = str_to_title(paste(set, "cell"))) %>%
  #   ggplot(aes(mode, shape = set, x = set, color= sd)) +
  #   geom_boxplot(outlier.shape = NA) +
  #   geom_jitter() +
  #   scale_shape_manual(values = c(2,19)) +
  #   labs(y = "Mode of the posterior distribution",
  #        color = "Standard deviation of\nthe posteiror distribution") +
  #   theme(legend.position = "none", axis.title.y = element_blank()) +
  #   scale_y_continuous(breaks = 0:10)
  
  # browser()
  
  plot <- intensity_data_plot + #plot_spacer() +
    MCMC_vs_intensity_data  +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) + 
    plot_layout(ncol = 1, heights = c(1,2)) & 
    theme(legend.position = "right")
  
  # mode_plot +
  # plot_layout(ncol = 2, heights = c(1,3), widths = c(3,1), guides = "collect")
  
  file <- switch(tolower(dataset),
                 "fixed 8tr" = "Figure2",
                 "live 8tr" = "Figure3",
                 "fixed kshv" = "Figure5",
                 "live kshv" = "Figure6")
  file <- paste0(file, "_MCMC_inference.png")
  
  ggsave(here(results_folder, file), plot, width = 8, height = 5)
  
}

file_path <- here("results", supplemental_folder)
if(!dir.exists(file_path ))dir.create(file_path)

MCMC_vs_data(Figure2_results, Figure2_data, file_path, "Fixed 8TR", normalize = T)  
MCMC_vs_data(Figure3_results, Figure3_data, file_path, "Live 8TR", normalize = T)  
MCMC_vs_data(Figure5_results, Figure5_data, file_path, "Fixed KSHV", normalize = T)  
MCMC_vs_data(Figure6_results, Figure6_data, file_path, "Live KSHV", same_mu = F, normalize = T)  

certain_example <- Figure5_results$all_chains %>% filter(chain == "chain1") %>% 
  ggplot(aes(d1, after_stat(count/sum(count)))) + 
  geom_bar(aes(fill = max(after_stat(count/sum(count)))), show.legend = F) + 
  labs(x = "Number of episomes in LANA dot", y = "Posterior probability") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  ylim(c(0,1)) + 
  scale_fill_viridis_c(limits = c(0, 1), direction = -1) 

uncertain_example <- Figure5_results$all_chains %>% filter(chain == "chain1") %>% 
  ggplot(aes(d47, after_stat(count/sum(count)))) + 
  geom_bar(aes(fill = max(after_stat(count/sum(count)))), show.legend = F) + 
  labs(x = "Number of episomes in LANA dot", y = "Posterior probability") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = 1:10) +
  scale_fill_viridis_c(limits = c(0, 1), direction = -1)

uncertain_example2 <- Figure5_results$all_chains %>% filter(chain == "chain1") %>% 
  ggplot(aes(d28, after_stat(count/sum(count)))) + 
  geom_bar(aes(fill = max(after_stat(count/sum(count)))), show.legend = F) + 
  labs(x = "Number of episomes in LANA dot", y = "Posterior probability") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = 0:13) +
  scale_fill_viridis_c(limits = c(0, 1), direction = -1)

ggsave(here(file_path, "certain_example.png"), certain_example, width = 3, height = 2)
ggsave(here(file_path, "uncertain_example.png"), uncertain_example, width = 3, height = 2)
ggsave(here(file_path, "intermediate_example.png"), uncertain_example2, width = 3, height = 2)

# Plot MLE grids for all conditions
plot_grid <- function(results, dataset, file_path){
  plot <- plot_grid_search(results$MLE_grid, simulation = F, prob = T, error_bars = F) & 
    plot_annotation(title = element_blank()) & 
    theme(plot.caption = element_blank(), plot.background = element_blank())
  
  file <- switch(tolower(dataset),
                 "fixed 8tr" = "Figure2",
                 "live 8tr" = "Figure3",
                 "fixed kshv" = "Figure5",
                 "live kshv" = "Figure6")
  file <- paste0(file, "_MLE_grid.png")
  
  
  ggsave(here(file_path, file), plot, width = 4.5, height = 3.5, bg = "white") 
  
}

plot_grid(Figure2_results, "Fixed 8TR", file_path)
plot_grid(Figure3_results, "Live 8TR", file_path)
plot_grid(Figure5_results, "Fixed KSHV", file_path)
plot_grid(Figure6_results, "Live KSHV", file_path)


#####
# Silhuoette plots to assess cluster quality
#####

pdf(here(file_path, "silhouette_plot.pdf"))
intensity_data <- rbind(select(Figure2_data$daughter_cell_data, cluster_id, total_cluster_intensity),
                        select(Figure2_data$mother_cell_data, cluster_id, total_cluster_intensity))
Figure2_sil <- calc_sil(intensity_data, Figure2_results$all_chains)
fviz_silhouette(Figure2_sil) + 
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  theme(panel.grid = element_blank())

intensity_data <- rbind(select(Figure3_data$daughter_cell_data, cluster_id, total_cluster_intensity),
                        select(Figure3_data$mother_cell_data, cluster_id, total_cluster_intensity))
Figure3_sil <- calc_sil(intensity_data, Figure3_results$all_chains)
fviz_silhouette(Figure3_sil) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  theme(panel.grid = element_blank())

intensity_data <- rbind(select(Figure5_data$daughter_cell_data, cluster_id, total_cluster_intensity),
                        select(Figure5_data$mother_cell_data, cluster_id, total_cluster_intensity))
Figure5_sil <- calc_sil(intensity_data, Figure5_results$all_chains)
fviz_silhouette(Figure5_sil) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  theme(panel.grid = element_blank())

# intensity_data <- rbind(select(Figure6_data$daughter_cell_data, cluster_id, total_cluster_intensity),
#                         select(Figure6_data$mother_cell_data, cluster_id, total_cluster_intensity))
Figure6_sil <- calc_sil(Figure6_data$daughter_cell_data, Figure6_results$all_chains)
fviz_silhouette(Figure6_sil) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  theme(panel.grid = element_blank())

Figure6_sil <- calc_sil(Figure6_data$mother_cell_data, Figure6_results$all_chains)
fviz_silhouette(Figure6_sil)+
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1))  +
  theme(panel.grid = element_blank())

dev.off()

#####
# Make supplemental figures to show the chain's performance

MCMC_performance <- function(results, results_folder, dataset, same_mu = T){
  
  # browser()
  
  if(!same_mu) results$all_chains <- results$all_chains %>% rename(mu = mu_d, tau = tau_d)
  
  inferred_mu <- round(DescTools::Mode(round(results$all_chains$mu, 1)))
  inferred_sigma <- round(DescTools::Mode(round(sqrt(1/results$all_chains$tau), 1)))
  
  p1 <- results$all_chains %>% 
    ggplot(aes(iteration, mu, color = chain)) + 
    geom_line() +
    geom_point(shape = NA) + 
    geom_hline(yintercept = inferred_mu) + 
    geom_text(data = data.frame(NA), aes(-Inf, inferred_mu), label = paste0("mu: ", inferred_mu),
              hjust = -0.1, vjust = -1, inherit.aes = F) +
    ylim(c(0, max(results$all_chains$mu))) + 
    theme(legend.position = c(0,0), legend.justification = c(-0.1,-0.1)) + 
    labs(x = "Iteration", y = bquote(mu))
  
  p1 <- ggMarginal(p1, groupColour = T, margins = "y")
  
  p2 <- results$all_chains %>% 
    ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
    geom_line() +
    geom_point(shape = NA) + 
    geom_hline(yintercept = inferred_sigma) + 
    geom_text(data = data.frame(NA), aes(-Inf, inferred_sigma), label = paste0("sigma: ", inferred_sigma),
              hjust = -0.1, vjust = -1, inherit.aes = F) +
    ylim(c(0, max(sqrt(1/results$all_chains$tau)))) + 
    theme(legend.position = "none") + 
    labs(x = "Iteration", y = bquote(sigma))
  
  p2 <- ggMarginal(p2, groupColour = T, margins = "y")
  
  df <- results$convergence_metrics
  if(!same_mu) df <- df %>% filter(grepl("d", name))
  p3 <- df %>% 
    mutate(ESS_bulk = ESS_bulk/(100000-5000),
           ESS_tail = ESS_tail/(100000-5000)) %>% 
    pivot_longer(!name, names_to = "metric") %>%
    mutate(name = factor(name, levels = rev(sort(unique(name))))) %>%
    arrange(name) %>%
    ggplot(aes(value, metric)) + 
    geom_text(aes(label = name), position = position_jitter(seed = 1), color = "gray") +
    geom_boxplot(fill = NA, outlier.shape = NA) +
    labs(x = "Metric", y = "") + 
    scale_y_discrete(labels = rev(c("Rhat", "Tail\nESS\nRatio", "Bulk\nESS\nRatio"))) + 
    coord_cartesian(c(0, 2)) 
  
  plot <- cowplot::plot_grid(p1, p2, p3, nrow = 1)
  
  file <- switch(tolower(dataset),
                 "fixed 8tr" = "Figure2",
                 "live 8tr" = "Figure3",
                 "fixed kshv" = "Figure5",
                 "live kshv" = "Figure6")
  if(!same_mu) file <- paste0(file, "_daughter")
  file <- paste0(file, "_MCMC_performance.png")
  
  ggsave(here(results_folder, file), plot, width = 10, height = 3, bg = "white")
  
  if(!same_mu){
    # browser()
    results$all_chains <- results$all_chains %>% select(-mu, -tau) %>% rename(mu = mu_m, tau = tau_m)
    
    inferred_mu <- round(DescTools::Mode(round(results$all_chains$mu,1)))
    inferred_sigma <- round(DescTools::Mode(round(sqrt(1/results$all_chains$tau),1)))
    
    p1 <- results$all_chains %>% 
      ggplot(aes(iteration, mu, color = chain)) + 
      geom_line() +
      geom_point(shape = NA) + 
      geom_hline(yintercept = inferred_mu) + 
      geom_text(data = data.frame(NA), aes(-Inf, inferred_mu), label = paste0("mu: ", inferred_mu),
                hjust = -0.1, vjust = -1, inherit.aes = F) +
      ylim(c(0, max(results$all_chains$mu))) + 
      theme(legend.position = c(0,0), legend.justification = c(-0.1,-0.1)) + 
      labs(x = "Iteration", y = bquote(mu))
    
    p1 <- ggMarginal(p1, groupColour = T, margins = "y")
    
    p2 <- results$all_chains %>% 
      ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
      geom_line() +
      geom_point(shape = NA) + 
      geom_hline(yintercept = inferred_sigma) + 
      geom_text(data = data.frame(NA), aes(-Inf, inferred_sigma), label = paste0("sigma: ", inferred_sigma),
                hjust = -0.1, vjust = -1, inherit.aes = F) +
      ylim(c(0, max(sqrt(1/results$all_chains$tau)))) + 
      theme(legend.position = "none") + 
      labs(x = "Iteration", y = bquote(sigma))
    
    p2 <- ggMarginal(p2, groupColour = T, margins = "y")
    
    df <- results$convergence_metrics
    if(!same_mu) df <- df %>% filter(name != "mu_d") %>% filter(grepl("m", name))
    p3 <- df %>% 
      mutate(ESS_bulk = ESS_bulk/(100000-5000),
             ESS_tail = ESS_tail/(100000-5000)) %>% 
      pivot_longer(!name, names_to = "metric") %>%
      mutate(name = factor(name, levels = rev(sort(unique(name))))) %>%
      arrange(name) %>%
      ggplot(aes(value, metric)) + 
      geom_text(aes(label = name), position = position_jitter(seed = 1), color = "gray") +
      geom_boxplot(fill = NA, outlier.shape = NA) +
      labs(x = "Metric", y = "") + 
      scale_y_discrete(labels = rev(c("Rhat", "Tail\nESS\nRatio", "Bulk\nESS\nRatio"))) + 
      coord_cartesian(c(0, 2))
    
    plot <- cowplot::plot_grid(p1, p2, p3, nrow = 1)
    
    file <- switch(tolower(dataset),
                   "fixed 8tr" = "Figure2",
                   "live 8tr" = "Figure3",
                   "fixed kshv" = "Figure5",
                   "live kshv" = "Figure6")
    if(!same_mu) file <- paste0(file, "_mother")
    file <- paste0(file, "_MCMC_performance.png")
    
    ggsave(here(results_folder, file), plot, width = 10, height = 3, bg = "white")
  }
  
}

MCMC_performance(Figure2_results, file_path, "Fixed 8TR")
MCMC_performance(Figure3_results, file_path, "Live 8TR")
MCMC_performance(Figure5_results, file_path, "Fixed KSHV")
MCMC_performance(Figure6_results, file_path, "Live KSHV", same_mu = F)
