## 
# This script contains two functions:
# 1) run_pipeline() runs the Analysis including MCMC and MLE
# 2) make_plots() creates figures for the results

# load libraries 
require(tidyverse)
require(here)
# For plotting
require(patchwork)
require(RColorBrewer)
require(ggExtra)
require(bayesplot)
# for markov chain convergence statistics
require(rstan)
# for fitting a poisson distribution
require(MASS)
require(furrr)
require(purrr)
require(foreach)


select <- dplyr::select


options(dplyr.summarise.inform = FALSE)

## Color-blind friendly color palette
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#####
# Function to read in data
#####

load_data <- function(mother_cell_file, daughter_cell_file){
    
  # Daughter cell data
  daughter_cell_data <- readxl::read_excel(here("data", "derived", daughter_cell_file))
  # Mother cell data
  mother_cell_data <- readxl::read_excel(here("data", "derived", mother_cell_file))
  
  # reformat the daughter cell data
  names(daughter_cell_data) <- gsub(" ", "_", tolower(names(daughter_cell_data)))
  daughter_cell_data <- rbind(
    select(daughter_cell_data, mother_cell_id, contains("daughter1")) %>% rename_with(~gsub("_daughter1",  "", .x)) %>% mutate(daughter_cell = 1),
    select(daughter_cell_data, mother_cell_id, contains("daughter2")) %>% rename_with(~gsub("_daughter2",  "", .x)) %>% mutate(daughter_cell = 2)
  ) %>% 
    select(mother_cell_id, daughter_cell, cluster, everything()) %>% 
    arrange(mother_cell_id, daughter_cell, cluster) %>% 
    filter(!is.na(total_cluster_intensity), total_cluster_intensity !=0) %>% 
    mutate(cluster_id  = paste0("d", 1:nrow(.)),
           cell_id = paste(mother_cell_id, daughter_cell, sep = "_"))
  
  # reformat the mother cell data
  names(mother_cell_data) <- gsub(" ", "_", tolower(names(mother_cell_data)))
  mother_cell_data <- mother_cell_data %>% 
    # remove clusters with no intensity data
    filter(!is.na(total_cluster_intensity)) %>%  
    mutate(cell_id = paste(image, cell, sep = "_"),
           cluster_id = paste0("m", 1:nrow(.))) 
  
  return(list(daughter_cell_data = daughter_cell_data, mother_cell_data = mother_cell_data))
}

#####
# Function to run statistical inference framework
#####


run_pipeline <- function(daughter_cell_data, mother_cell_data, results_folder, 
                         n_iterations = 100000, burn_in = 5000, MLE_n_samples = 100, 
                         same_mu = T, overwrite = F, n_prior = list("pois", 1), parallel = F,
                         just_Pr = F){
  
  
  # Make results folder
  if(!file.exists(results_folder)) dir.create(results_folder, showWarnings = F)
  
  #####
  # Gibbs sampling to infer the number of episomes per cluster in each daughter cell
  #####
  
  #  set initial conditions:
  # We will assume that the value of mu is bounded between the ranges observed in the data. 
  # Thus, we will select initial conditions within this range. We will only try one value of tau that is small enough 
  # (meaning sigma is large enough) that the chain will converge based on chains run in simulated data. 
  # We will determine the initial conditions for nk based on the observed data and the initial condition for mu.
  tau0 <- 1e-5 
  # Combine intensity data for both mother cells and daughter cells
  intensity_data <- rbind(daughter_cell_data %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "daughter"), 
                          mother_cell_data  %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "mother"))
  
  
  file <- here(results_folder, "MCMC_samples_per_cluster.RData")
  if(file.exists(file) & !overwrite){
    load(file)
  }else{
    if(same_mu){
      ICs <- intensity_data %>% 
        mutate(ratio = total_cluster_intensity/min_episome_in_cluster) %>% filter(!is.na(ratio)) %>% pull(ratio) %>% summary
      mu0 <- ICs[c(2,5)]
      
      cat("\nRunning Gibbs with the following parameters", 
          "\ntau0:", tau0, "\nmu0:", mu0, 
          "\nnumber_of_iterations:", n_iterations,
          "\nn_clusters:", nrow(intensity_data),
          "\n===")
      
      # Case when the mean intensity of a single episome is the same for mother and daughter cells
      # run 2 chains
      intensities <- intensity_data$total_cluster_intensity
      names(intensities) <- intensity_data$cluster_id
      chain1 <- run_gibbs(tau0, mu0[1], intensities, n_iterations, n_prior = n_prior)
      chain2 <- run_gibbs(tau0, mu0[2], intensities, n_iterations, n_prior = n_prior)
    }else{
      # Case when the mean intensity differs between mother and daughter cells
      cat("\nRunning Gibbs separately for mother and daughter cells")
      ICs <- daughter_cell_data %>% 
        mutate(ratio = total_cluster_intensity/min_episome_in_cluster) %>% filter(!is.na(ratio)) %>% pull(ratio) %>% summary
      mu0 <- ICs[c(2,5)]
      cat("\nDaughter cell parameters:", 
          "\ntau0:", tau0, "\nmu0:", mu0, 
          "\nnumber_of_iterations:", n_iterations,
          "\nn_clusters:", nrow(daughter_cell_data),
          "\n===")
      
      # start with daughters
      intensities <- daughter_cell_data$total_cluster_intensity
      names(intensities) <- daughter_cell_data$cluster_id
      chain1_d <- run_gibbs(tau0, mu0[1], intensities, n_iterations, n_prior = n_prior)
      chain2_d <- run_gibbs(tau0, mu0[2], intensities, n_iterations, n_prior = n_prior)
      
      ICs <- mother_cell_data %>% 
        mutate(ratio = total_cluster_intensity/min_episome_in_cluster) %>% filter(!is.na(ratio)) %>% pull(ratio) %>% summary
      mu0 <- ICs[c(2,5)]
      cat("\nMother cell parameters:", 
          "\ntau0:", tau0, "\nmu0:", mu0, 
          "\nnumber_of_iterations:", n_iterations,
          "\nn_clusters:", nrow(mother_cell_data),
          "\n===")
      
      # start with daughters
      intensities <- mother_cell_data$total_cluster_intensity
      names(intensities) <- mother_cell_data$cluster_id
      chain1_m <- run_gibbs(tau0, mu0[1], intensities, n_iterations, n_prior = n_prior)
      chain2_m <- run_gibbs(tau0, mu0[2], intensities, n_iterations, n_prior = n_prior)
      
      chain1 <- merge(chain1_d, chain1_m, by = "iteration", suffixes = c("_d", "_m"))
      chain2 <- merge(chain2_d, chain2_m, by = "iteration", suffixes = c("_d", "_m"))
    }
    
    
    save(chain1, chain2, file = file)
  }
  
  # combine all chains together and remove burn-in period
  all_chains <- rbind(chain1, chain2) %>%
    mutate(chain = as.factor(rep(c("chain1", "chain2"), each = n_iterations)))  %>%
    filter(iteration > burn_in)
  
  # Find inferred value of mu and sigma based off of the mode of their joint posterior
  if(same_mu){
    inferred_mu <- median(all_chains$mu)
    inferred_sigma2 <- median(1/all_chains$tau)
    inferred_sigma <- median(sqrt(1/all_chains$tau))
    
    cat("\nResults of Gibbs", '\nmean:', inferred_mu, '\nvariance:', inferred_sigma2, '\nstandard deviation:', inferred_sigma,
        "\n===")
  }else{
    inferred_mu <- median(all_chains$mu_d)
    inferred_sigma2 <- median(1/all_chains$tau_d)
    inferred_sigma <- median(sqrt(1/all_chains$tau_d))
    
    cat("\nResults of Gibbs in daughter cells", '\nmean:', inferred_mu, '\nvariance:', inferred_sigma2, '\nstandard deviation:', inferred_sigma,
        "\n===")
    
    inferred_mu <- median(all_chains$mu_m)
    inferred_sigma2 <- median(1/all_chains$tau_m)
    inferred_sigma <- median(sqrt(1/all_chains$tau_m))
    
    cat("\nResults of Gibbs in mother cells", '\nmean:', inferred_mu, '\nvariance:', inferred_sigma2, '\nstandard deviation:', inferred_sigma,
        "\n===")
  }
  
  
  #####
  # Get episome per cell at each iteration
  #####
  
  file <- here(results_folder, "MCMC_samples_per_cell.RData")
  if(file.exists(file) & !overwrite){
    load(file)
  }else{
    cell_samples_long <- all_chains %>% 
      select(iteration, chain, grep("[0-9]", names(all_chains), perl = T)) %>% 
      pivot_longer(!c(iteration, chain), names_to = "cluster_id", values_to = "nk") %>% 
      left_join(intensity_data, by = "cluster_id") %>% 
      group_by(chain, iteration, set, cell_id) %>% 
      summarise(number_of_episomes = sum(nk)) %>% 
      ungroup()
    
    daughter_cell_samples <- cell_samples_long %>% filter(set == "daughter") %>% 
      select(-set) %>% 
      pivot_wider(names_from = cell_id, values_from = number_of_episomes)
    
    mother_cell_samples <- cell_samples_long %>% filter(set == "mother") %>% 
      select(-set) %>% 
      pivot_wider(names_from = cell_id, values_from = number_of_episomes)
    
    save(cell_samples_long, daughter_cell_samples, mother_cell_samples, file = file) 
  }
  
  #####
  # Convergence Statistics for MCMC
  #####
  cat("\nChecking convergence of MCMC","\n===")
  file <- here(results_folder, "MCMC_convergence.RData")
  if(file.exists(file)  & !overwrite){
    load(file)
  }else{
    convergence <- convergence_results(all_chains)
    save(convergence, file = file)  
  }
  
  #####
  # Fitting distribution of initial number of episomes
  #####
  
  # we will construct the distribution of initial number of episomes using the full set of samples from 
  # the posterior distribution
  cat("\nFitting PMF for X0")
  lambda <- fitdistr(cell_samples_long %>% filter(set == "mother") %>% pull(number_of_episomes), "Poisson")$estimate
  cat("\nlambda:", lambda,"\n===")
  
  #####
  # Maximum likelihood estimation with grid search
  #####
  
  cat("\nRunning Maximum Likelihood on", MLE_n_samples, "samples from MCMC")
  
  # randomly choose samples from the markov chain to use as inferred value of number of episomes per cell:
  samples <- sample(unique(all_chains$iteration), MLE_n_samples, replace = FALSE)
  
  # reformat samples so that we can apply the grid search function: 
  # need a column for mother cell id, number of cells in daughter cell 1, and number of cells in daughter cecll 2
  sampled_experimental_data <- cell_samples_long %>%
    filter(chain == "chain1", set == "daughter", iteration %in% samples) %>% 
    # pull in mother and daughter cell ids
    left_join(distinct(daughter_cell_data, cell_id, mother_cell_id), by = c("cell_id")) %>%
    group_by(mother_cell_id, iteration) %>% 
    summarise(X1 = max(number_of_episomes), X2 = ifelse(n() > 1, min(number_of_episomes), 0)) %>% 
    ungroup()
  
  file <- here(results_folder, "MLE_with_uncertainty.RData")
  if(file.exists(file) & !overwrite){
    load(file)
  }else{
    if(parallel){
      #create the cluster
      my.cluster <- parallel::makeCluster(2, type = "PSOCK")
      
      #check cluster definition (optional)
      print(my.cluster)
      
      #register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      #check if it is registered (optional)
      cat("cluster registered:")
      print(foreach::getDoParRegistered())
      
    }
    
    # iterate through the sampled data and apply a grid search to each one
    j <- 1
    MLE_with_uncertainty <- foreach(i = samples, .packages = c("tidyverse", "here")) %dopar% {
      cat("\n", j)
      source(here("scripts","likelihood_functions.R"))
      j <- j+1
      temp <- sampled_experimental_data %>% 
        filter(iteration == i) %>%
        select(id =  mother_cell_id,  X1, X2) %>%
        #  run a  grid search with a grid ranging by 0.01 and a PMF inferred from the experimental data
        run_grid_search(viz = F, increment = 0.01,
                        lambda = lambda,
                        known_X0 = F, CI = F)
      
      # scale the log likelihood according to the maximum  value and then convert to a probability
      temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                          probability = likelihood/sum(likelihood))
      temp
    }
    save(MLE_with_uncertainty, file = file)
    
    if(parallel & !just_Pr) parallel::stopCluster(my.cluster)
  }
  
  cat("\n===")
  
  MLE_with_uncertainty_df <- MLE_with_uncertainty %>% bind_rows() %>% 
    select(Pr, Ps, probability) %>% 
    group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
    ungroup %>% 
    mutate(probability = probability/sum(probability),
           log_likelihood = log(probability))
  
  MLE_with_uncertainty_CI <- calculate_CI(MLE_with_uncertainty_df)
  
  MLE_uncertainty_grid <- list(grid_search = MLE_with_uncertainty_df, 
                               estimates = MLE_with_uncertainty_CI$estimates,
                               top_95 = MLE_with_uncertainty_CI$probs)
  
  cat("\nDone") 
  
  if(!just_Pr){
    return(list(all_chains = all_chains, 
                daughter_cell_samples = daughter_cell_samples, 
                mother_cell_samples = mother_cell_samples, 
                convergence_metrics = convergence, 
                X0_lambda = lambda, 
                MLE_grid = MLE_uncertainty_grid))
  }else{
    if(parallel & !exists("my.cluster")){
      #create the cluster
      my.cluster <- parallel::makeCluster(2, type = "PSOCK")
      
      #check cluster definition (optional)
      print(my.cluster)
      
      #register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      #check if it is registered (optional)
      cat("cluster registered:")
      print(foreach::getDoParRegistered())
    }
    
    file <- here(results_folder, "MLE_Pr_with_uncertainty.RData")
    if(file.exists(file) & !overwrite){
      load(file)
    }else{
      print(dim(sampled_experimental_data))
      ## Estimate Pr only from sum of daughter cells
      j <- 1
      MLE_with_uncertainty <- foreach(i = samples, .packages = c("tidyverse", "here")) %dopar% {
        cat("\n", j)
        print("in loop")
        source(here("scripts","likelihood_functions.R"))
        j <- j+1
        temp <- sampled_experimental_data %>% 
          filter(iteration == i) %>%
          select(id =  mother_cell_id,  X1, X2) %>%
          #  run a  grid search with a grid ranging by 0.01 and a PMF inferred from the experimental data
          run_grid_search(viz = F, increment = 0.01,
                          lambda = lambda,
                          known_X0 = F, CI = F, just_Pr = T)
        
        # scale the log likelihood according to the maximum  value and then convert to a probability
        temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                            probability = likelihood/sum(likelihood))
        temp
      }
      save(MLE_with_uncertainty, file = file)
      
      if(parallel) parallel::stopCluster(my.cluster)
      
    }
    
    MLE_with_uncertainty_df <- MLE_with_uncertainty %>% bind_rows() %>% 
      select(Pr, Ps, probability) %>% 
      group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
      ungroup %>% 
      mutate(probability = probability/sum(probability),
             log_likelihood = log(probability))
    
    MLE_with_uncertainty_CI <- calculate_CI(MLE_with_uncertainty_df)
    
    MLE_Pr_uncertainty_grid <- list(grid_search = MLE_with_uncertainty_df, 
                                    estimates = MLE_with_uncertainty_CI$estimates,
                                    top_95 = MLE_with_uncertainty_CI$probs)
    
    return(list(all_chains = all_chains, 
                daughter_cell_samples = daughter_cell_samples, 
                mother_cell_samples = mother_cell_samples, 
                convergence_metrics = convergence, 
                X0_lambda = lambda, 
                MLE_grid = MLE_uncertainty_grid,
                MLE_Pr_grid = MLE_Pr_uncertainty_grid, 
                sampled_data = sampled_experimental_data))
  }
  
  
}

#####
# Function to generate figures
#####

make_plots <- function(pipeline_output, daughter_cell_data, mother_cell_data, results_folder){
  
  all_chains <- pipeline_output$all_chains
  daughter_cell_samples <- pipeline_output$daughter_cell_samples
  mother_cell_samples <- pipeline_output$mother_cell_samples
  lambda <- pipeline_output$X0_lambda
  convergence <- pipeline_output$convergence_metrics
  MLE_grid <- pipeline_output$MLE_grid
  
  intensity_data <- rbind(daughter_cell_data %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "daughter"), 
                          mother_cell_data  %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "mother"))
  
  
  #####
  # Plot intensity data
  #####
  
  intensity_plots <- mother_cell_data %>%
    group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>%
    ggplot(aes(total_cluster_intensity, y = after_stat(count / sum(count)))) +
    geom_histogram(aes(fill = "non-dividing cells", color = "non-dividing cells"), alpha = 0.5) +
    geom_histogram(aes(x = 2*total_cluster_intensity, fill = "perfect replication\nfrom mother cells",
                       color = "perfect replication\nfrom mother cells"), alpha = 0.5) +
    geom_histogram(data = daughter_cell_data %>% group_by(mother_cell_id) %>%
                     mutate(int = sum(total_cluster_intensity)),
                   aes(int, fill = "daughter cells", color = "daughter cells"), alpha = 0.5) +
    labs(x = "cluster intensity", fill = "cell set", y= "frequency") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = c(1,1), legend.justification = c(1.1,1.1))  +
    guides(color = "none") +
    labs(title = "Distribution of whole cell intensities") +
    
    mother_cell_data %>%
    group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>%
    ggplot() +
    geom_boxplot(aes(x = total_cluster_intensity, fill = "non-dividing cells", y = "non-dividing cells"), alpha = 0.5) +
    geom_boxplot(aes(x = 2*total_cluster_intensity, fill = "perfect replication from mother cells", y = "perfect replication\nfrom mother cells"), alpha = 0.5) +
    geom_boxplot(data = daughter_cell_data %>% group_by(mother_cell_id) %>%
                   mutate(int = sum(total_cluster_intensity)),
                 aes(int, fill = "daughter cells", y = "daughter cells"), alpha = 0.5) +
    labs(x = "cluster intensity", fill = "cell set") +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    
    plot_layout(ncol = 1, heights = c(3,1)) &
    labs(x = "Intensity in cell(s)") &
    scale_color_manual(values = safe_colorblind_palette) &
    scale_fill_manual(values = safe_colorblind_palette)
  
  ggsave(here(results_folder, "intensity_data.pdf"), intensity_plots, width = 7, height = 7)
  
  #####
  # Plot results of Gibbs Sampling
  #####
  same_mu  <- !any(grepl("mu_d", names(all_chains)))
  if(same_mu){
    
    inferred_mu <- median(all_chains$mu)
    inferred_sigma2 <- median(1/all_chains$tau)
    inferred_sigma <- median(sqrt(1/all_chains$tau))
    
    mu <- all_chains %>% 
      ggplot(aes(iteration, mu, color = chain)) + 
      geom_line() +
      geom_point(shape = NA) + 
      geom_hline(yintercept = inferred_mu, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_mu, label = str_interp("${round(inferred_mu,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") + 
      ylim(c(0, max(all_chains$mu))) +
      labs(title = "Trace of the mean (\u00b5)", 
           y = bquote("mean of the distribution of\nfluorescence intensity for one episome (\u00b5)"))
    
    mu <- ggMarginal(mu, margins = "y", groupColour = T)
    
    sigma2 <- all_chains %>% 
      ggplot(aes(iteration, 1/tau, color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma2, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma2, label = str_interp("${round(inferred_sigma2,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(1/all_chains$tau))) + 
      labs(title = "Trace of the variance (\u03c3\u00b2)", 
           y = bquote("variance of the distribution of\nfluorescence intensity for one episome (\u03c3\u00b2)"))
    
    sigma2 <- ggMarginal(sigma2, margins = "y", groupColour = T)
    
    sigma <- all_chains %>% 
      ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma, label = str_interp("${round(inferred_sigma,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(sqrt(1/all_chains$tau)))) + 
      labs(title = "Trace of the standard deviation (\u03c3)", 
           y = bquote("standard deviation of the distribution of\nfluorescence intensity for one episome (\u03c3)"))
    
    sigma <- ggMarginal(sigma, margins = "y", groupColour = T)
    
    
    # Based on the median mu and sigma above, the inferred distribution of intensities for a single cluster is
    
    x <- seq(0, max(intensity_data$total_cluster_intensity), length.out = 500)
    intensity_distribution <- ggplot(tibble(x, y = dnorm(x, inferred_mu, inferred_sigma)), aes(x, y)) + 
      geom_line() + 
      labs(x = "Fluorescence Intensity", 
           y = "probability density", 
           title = "Inferred distribution of fluorescence intensity\nof a single cluster",
           subtitle = str_interp("mean = ${round(inferred_mu,2)}, standard deviation = ${round(inferred_sigma, 2)}"))
  }else{
    ## Daughter cell plots
    inferred_mu = median(all_chains$mu_d)
    inferred_sigma2 = median(1/all_chains$tau_d)
    inferred_sigma = median(sqrt(1/all_chains$tau_d))
    
    mu_d <- all_chains %>% 
      ggplot(aes(iteration, mu_d, color = chain)) + 
      geom_line() + 
      geom_point(shape = NA) + 
      geom_hline(yintercept = inferred_mu, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_mu, label = str_interp("${round(inferred_mu,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") + 
      ylim(c(0, max(all_chains$mu_d))) +
      labs(title = "Trace of the mean in daughter cells (\u00b5)", 
           y = bquote("mean of the distribution of\nfluorescence intensity for one episome (\u00b5)")) 
    
    mu_d <- ggMarginal(mu_d, margins = "y", groupColour = T)
    
    sigma2_d <- all_chains %>% 
      ggplot(aes(iteration, 1/tau_d, color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma2, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma2, label = str_interp("${round(inferred_sigma2,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(1/all_chains$tau_d))) + 
      labs(title = "Trace of the variance in daughter cells (\u03c3\u00b2)", 
           y = bquote("variance of the distribution of\nfluorescence intensity for one episome (\u03c3\u00b2)"))
    
    sigma2_d <- ggMarginal(sigma2_d, margins = "y", groupColour = T)
    
    sigma_d <- all_chains %>% 
      ggplot(aes(iteration, sqrt(1/tau_d), color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma, label = str_interp("${round(inferred_sigma,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(sqrt(1/all_chains$tau_d)))) + 
      labs(title = "Trace of the standard deviation in daughter cells (\u03c3)", 
           y = bquote("standard deviation of the distribution of\nfluorescence intensity for one episome (\u03c3)"))
    
    sigma_d <- ggMarginal(sigma_d, margins = "y", groupColour = T)
    
    x <- seq(0, max(daughter_cell_data$total_cluster_intensity), length.out = 500)
    intensity_distribution_d <- ggplot(tibble(x, y = dnorm(x, inferred_mu, inferred_sigma)), aes(x, y)) + 
      geom_line() + 
      labs(x = "Fluorescence Intensity", 
           y = "probability density", 
           title = "Inferred distribution of fluorescence intensity\nof a single cluster in daughter cells",
           subtitle = str_interp("mean = ${round(inferred_mu,2)}, standard deviation = ${round(inferred_sigma, 2)}"))
    
    
    ## Mother cell plots
    inferred_mu = median(all_chains$mu_m)
    inferred_sigma2 = median(1/all_chains$tau_m)
    inferred_sigma = median(sqrt(1/all_chains$tau_m))
    
    mu_m <- all_chains %>% 
      ggplot(aes(iteration, mu_m, color = chain)) + 
      geom_line() + 
      geom_point(shape = NA) + 
      geom_hline(yintercept = inferred_mu, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_mu, label = str_interp("${round(inferred_mu,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") + 
      ylim(c(0, max(all_chains$mu_m))) +
      labs(title = "Trace of the mean in mother cells (\u00b5)", 
           y = bquote("mean of the distribution of\nfluorescence intensity for one episome (\u00b5)")) 
    
    mu_m <- ggMarginal(mu_m, margins = "y", groupColour = T)
    
    sigma2_m <- all_chains %>% 
      ggplot(aes(iteration, 1/tau_m, color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma2, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma2, label = str_interp("${round(inferred_sigma2,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(1/all_chains$tau_m))) + 
      labs(title = "Trace of the variance in mother cells (\u03c3\u00b2)", 
           y = bquote("variance of the distribution of\nfluorescence intensity for one episome (\u03c3\u00b2)"))
    
    sigma2_m <- ggMarginal(sigma2_m, margins = "y", groupColour = T)
    
    sigma_m <- all_chains %>% 
      ggplot(aes(iteration, sqrt(1/tau_m), color = chain)) + 
      geom_line() +
      geom_point(shape = NA) +
      geom_hline(yintercept = inferred_sigma, lty = "dashed") + 
      geom_text(data = data.frame(NA), aes(min(all_chains$iteration), inferred_sigma, label = str_interp("${round(inferred_sigma,2)}")),
                inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
      theme(legend.position = "bottom") +
      ylim(c(0, max(sqrt(1/all_chains$tau_m)))) + 
      labs(title = "Trace of the standard deviation in mother cells (\u03c3)", 
           y = bquote("standard deviation of the distribution of\nfluorescence intensity for one episome (\u03c3)"))
    
    sigma_m <- ggMarginal(sigma_m, margins = "y", groupColour = T)
    
    x <- seq(0, max(mother_cell_data$total_cluster_intensity), length.out = 500)
    intensity_distribution_m <- ggplot(tibble(x, y = dnorm(x, inferred_mu, inferred_sigma)), aes(x, y)) + 
      geom_line() + 
      labs(x = "Fluorescence Intensity", 
           y = "probability density", 
           title = "Inferred distribution of fluorescence intensity\nof a single cluster in mother cells",
           subtitle = str_interp("mean = ${round(inferred_mu,2)}, standard deviation = ${round(inferred_sigma, 2)}"))
  }
  
  
  # How do these inferences relate to the observed intensity data?
  cluster_modes <- all_chains %>% 
    filter(chain == "chain1") %>% 
    select(grep("[0-9]", names(all_chains))) %>% 
    pivot_longer(everything(), names_to = "cluster_id", values_to = "n_episomes") %>% 
    count(cluster_id, n_episomes) %>% 
    group_by(cluster_id) %>% 
    filter(n == max(n)) %>% 
    ungroup %>% 
    select(-n) %>% 
    rename(mode = n_episomes) %>% 
    mutate(set = ifelse(grepl("d", cluster_id), "daughter", "mother"))
  
  
  MCMC_vs_intensity_data <- intensity_data %>% 
    merge(cluster_modes) %>% 
    ggplot(aes(total_cluster_intensity, as.factor(mode))) + 
    geom_jitter(aes(color = as.factor(min_episome_in_cluster))) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(title = "Distribution of cluster intensities",
         y = "Mode of the posterior distribution", 
         color = "Min # episome in cluster",
         x = "Total  Cluster Intensity") + 
    theme(legend.position = c(1,0), legend.justification = c(1,0), legend.background = element_blank()) + 
    facet_wrap(~set) +
    scale_color_manual(values = safe_colorblind_palette)
  
  # look at posterior probabilities
  temp_df <- merge(daughter_cell_samples, mother_cell_samples, by = c("chain", "iteration")) %>% 
    filter(chain == "chain1") %>% 
    pivot_longer(!c(iteration, chain), names_to = "cell_id", values_to = "n_episomes") %>% 
    group_by(cell_id) %>% add_count(n_episomes) %>% 
    mutate(mode = n_episomes[which.max(n)]) %>% 
    ungroup 
  
  check_cell_posteriors <- temp_df %>% 
    ggplot(aes(n_episomes, group = cell_id)) + 
    geom_density() + facet_wrap(~mode, scales = "free_y", labeller = "label_both") + 
    labs(x = "number of episomes per cell",  title = "Posterior probabilities for number of episomes per cell")
  
  # Distribution of episomes per cell in non-dividing cells
  mother_cell_samples_long <- mother_cell_samples %>% 
    filter(chain == "chain1") %>% 
    pivot_longer(!c(iteration, chain), names_to = "cell_id", values_to = "number_of_episomes") 
  
  X0_distribution  <- mother_cell_samples_long %>% 
    count(number_of_episomes) %>% mutate(freq = n/sum(n)) %>% 
    ggplot(aes(number_of_episomes, freq)) + 
    geom_bar(stat = "identity") +
    geom_point(data = tibble(number_of_episomes = 1:max(mother_cell_samples_long$number_of_episomes), 
                             y= dpois(number_of_episomes, lambda = lambda)),
               aes(y  = y)) + 
    labs(title = "Sampled distribtion of initial number of episomes",
         subtitle = str_interp("X0 ~ poisson(lambda = ${round(lambda, 2)}), median = ${median(mother_cell_samples_long$number_of_episomes)}"),
         x = "Number of episomes in Non-dividing cells",
         y = "Frequency") + 
    scale_x_continuous(breaks = 0:100)
  
  # Compare number of episomes per mother cells to total number of episomes in daughter cells
  temp_df <- daughter_cell_samples %>%
    filter(chain == "chain1") %>%
    pivot_longer(!c(chain, iteration), names_to = "cell_id", values_to = "n_episomes") %>% 
    count(cell_id, n_episomes) %>% 
    group_by(cell_id) %>% 
    filter(n == max(n)) %>% 
    select(-n) %>% 
    ungroup %>% 
    merge(distinct(daughter_cell_data, mother_cell_id, cell_id, daughter_cell)) %>% 
    group_by(mother_cell_id) %>% 
    summarise(X1 = max(n_episomes), X2 = ifelse(n() > 1, min(n_episomes), 0), total = sum(n_episomes)) 
  
  
  pair_histogram <- temp_df %>% mutate(pair = paste0("(", X1, ",", X2, ")")) %>%
    count(pair) %>%
    arrange(desc(n)) %>%
    slice(1:10) %>%
    ggplot(aes(pair, n)) +
    geom_bar(stat = "identity") +
    labs(x = "Number of of episomes in each daughter cell (X1, X2)",
         y = "Number of daughter cell pairs")
  
  
  distribution_df <- rbind(
    temp_df %>%
      select(!starts_with("X")) %>%
      pivot_longer(!mother_cell_id, names_to = "type", values_to = "number_of_episomes") %>%
      rename(cell_id = mother_cell_id),
    
    mother_cell_samples_long %>% filter(chain == "chain1") %>%
      count(cell_id, number_of_episomes) %>%
      group_by(cell_id) %>% filter(n == max(n)) %>%
      ungroup() %>% select(-n)  %>%
      mutate(type = "mother cell"),
    
    mother_cell_samples_long %>% filter(chain == "chain1") %>%
      count(cell_id, number_of_episomes) %>%
      group_by(cell_id) %>% filter(n == max(n)) %>%
      ungroup() %>% select(-n)  %>%
      mutate(number_of_episomes = 2*number_of_episomes) %>% 
      mutate(type = "perfect replication")
  ) %>%
    mutate(type = factor(type, levels= c("mother cell", "perfect replication", "total"),
                         labels =  c("non-dividing cells", "perfect replication", "daughter cell pairs"))) %>%
    count(type, number_of_episomes) %>%
    group_by(type) %>%
    mutate(freq = n/sum(n))
  
  
  cell_histograms <- distribution_df %>%
    ggplot(aes(number_of_episomes, freq, fill = type, color = type)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.5, position = "identity") +
    labs(x = "number of episomes", y = "Frequncy (%)", fill = "Cell Set") +
    scale_x_continuous(breaks = c(0:100)) + 
    scale_y_continuous(labels = scales::percent) + 
    scale_color_manual(values = safe_colorblind_palette) +
    scale_fill_manual(values = safe_colorblind_palette) + 
    guides(color = "none")
  
  
  cell_boxplots <- distribution_df %>% 
    ggplot(aes(number_of_episomes, type, fill = type)) +
    geom_boxplot(alpha = 0.5) +
    labs(x = "number of episomes", y = "Cell Set") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0:100)) + 
    scale_fill_manual(values = safe_colorblind_palette) 
  
  inference_per_cell <- cell_histograms + cell_boxplots + 
    plot_layout(ncol = 1, heights = c(3,1))
  
  ggsave(here(results_folder, "inference_per_cell.pdf"), inference_per_cell, width = 7, height = 7)
  
  #####
  # Convergence plots
  #####
  
  # Rhat values near 1 indicate equilibrium
  Rhat <- mcmc_rhat(convergence$Rhat) +
    ggtitle("Rhat")
  
  n_samples <- all_chains %>% filter(chain == "chain1") %>% nrow()
  
  ESS_bulk <- mcmc_neff(convergence$ESS_bulk/n_samples) + 
    ggtitle("Bulk Effective Sample Size Ratio")
  
  ESS_tail <- mcmc_neff(convergence$ESS_tail/n_samples) + 
    ggtitle("Tail Effective Sample Size Ratio")
  
  #####
  # Maximum Likelihood Estimation Results
  #####
  
  # browser()
  MLE_plot <- plot_grid_search(MLE_grid, simulation = F, prob = T) 
  
  ggsave(here(results_folder, "MLE_grid.pdf"), MLE_plot, width = 7, height = 7)
  
  if(same_mu){
    plot_list <- list(mu, sigma, sigma2, intensity_distribution, MCMC_vs_intensity_data,
                      check_cell_posteriors, X0_distribution, pair_histogram,
                      Rhat, ESS_bulk, ESS_tail)  
  }else{
    plot_list <- list(mu_d, sigma_d, sigma2_d, intensity_distribution_d, 
                      mu_m, sigma_m, sigma2_m, intensity_distribution_m,
                      MCMC_vs_intensity_data, check_cell_posteriors, X0_distribution, pair_histogram,
                      Rhat, ESS_bulk, ESS_tail)
  }
  
  
  ggsave(here(results_folder, "analysis_figures.pdf"), 
         gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
         width = 7, height = 7)
  
}

#####
# Function for figure plots
#####

fancy_figures <- function(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder){
  
  
  # Compare number of episomes per mother cells to total number of episomes in daughter cells
  temp_df <- daughter_cell_samples %>%
    filter(chain == "chain1") %>%
    pivot_longer(!c(chain, iteration), names_to = "cell_id", values_to = "n_episomes") %>% 
    separate(cell_id, into = c("mother_cell_id", "daughter"), sep = "_") %>% 
    pivot_wider(names_from = "daughter", values_from = "n_episomes", values_fill = 0) %>% 
    count(mother_cell_id, `1`, `2`) %>% 
    group_by(mother_cell_id) %>% 
    mutate(p = n/sum(n)) %>% 
    filter(n == max(n)) %>% 
    select(-n) %>% 
    # ungroup %>% 
    # group_by(mother_cell_id) %>%
    pivot_longer(c(`1`, `2`), values_to = "n_episomes") %>% 
    summarise(X1 = max(n_episomes), X2 = ifelse(n() > 1, min(n_episomes), 0), total = sum(n_episomes), p = unique(p)) %>%
    mutate(pair = paste0("(", X1, ",", X2, ")")) 
  
  pair_levels <- temp_df %>% count(pair) %>% arrange(desc(n)) %>% pull(pair)
  if(length(pair_levels) > 10){
    pair_labels <- c(pair_levels[1:10], rep("other", length(pair_levels)-10))  
  }else{
    pair_labels <- pair_levels
  }
  
  
  
  intensity_scatter <- daughter_cell_data %>% 
    group_by(mother_cell_id, cell_id, daughter_cell) %>% 
    summarise(total_intensity = sum(total_cluster_intensity)) %>% 
    group_by(mother_cell_id) %>% 
    summarise(daughter_1 = max(total_intensity), daughter_2= ifelse(n() > 1, min(total_intensity), 0)) %>% 
    left_join(temp_df, by = "mother_cell_id") %>% 
    mutate(pair = replace_na(pair, "other"),
           pair = factor(pair, levels = pair_levels, labels = pair_labels)) %>% 
    ggplot(aes(daughter_1, daughter_2, color = pair, size = p)) + 
    geom_abline(color = "gray", lty = "dashed") + 
    geom_point() + 
    labs(x = "Total intensity of daughter cell with more LANA dots",
         y = "Total intensity of daughter cell with fewer LANA dots") + 
    tune::coord_obs_pred() + 
    scale_color_manual(values = c(safe_colorblind_palette[2:12], rep("gray", 10) ))
  # theme(legend.position = "none")
  
  # ggsave(here(results_folder, "intensity_scatter_plot.png"), intensity_scatter, width = 5, height = 5)
  
  intensity_distribution <- mother_cell_data %>%
    group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>%
    ggplot(aes(total_cluster_intensity, y = after_stat(count / sum(count)))) +
    geom_histogram(aes(fill = "non-dividing cells", color = "non-dividing cells"), alpha = 0.5) +
    geom_histogram(data = daughter_cell_data %>% group_by(mother_cell_id) %>%
                     mutate(int = sum(total_cluster_intensity)),
                   aes(int, fill = "daughter cells", color = "daughter cells"), alpha = 0.5) +
    labs(x = "cluster intensity", fill = "cell set", y= "frequency") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = c(1,1), legend.justification = c(1.1,1.1))  +
    guides(color = "none") +
    labs(x = "Intensity in cell(s)") +
    scale_color_manual(values = safe_colorblind_palette) +
    scale_fill_manual(values = safe_colorblind_palette)
  
  # ggsave(here(results_folder, "intensity_distribution_plot.png"), intensity_distribution, width = 5, height = 4)
  
  
  pair_histogram <- temp_df  %>%
    count(pair, X1, X2) %>%
    arrange(desc(n)) %>%
    slice(1:10) %>%
    arrange(n, desc(X1), desc(X2)) %>% 
    mutate(pair = fct_inorder(pair)) %>% 
    ggplot(aes(pair, n, fill = pair)) +
    geom_bar(stat = "identity") +
    labs(x = "Number of of episomes in each daughter cell (X1, X2)",
         y = "Number of daughter\ncell pairs") + 
    scale_fill_manual(values = rev(safe_colorblind_palette[2:11])) +
    coord_flip() + 
    theme(legend.position = "none")
  
  figure <- intensity_distribution /
    (intensity_scatter + pair_histogram + 
       plot_layout(ncol = 2, widths = c(2,1))) + 
    plot_layout(heights = c(1,1.5))
  
  ggsave(here(results_folder, "fancy_figure.png"), figure, width = 10.5, height = 11)
  
  
}

figures <- function(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder){
  
  ## Color-blind friendly color palette
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#661100","#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255",  "#6699CC", "#888888")
  
  # Compare number of episomes per mother cells to total number of episomes in daughter cells
  temp_df <- daughter_cell_samples %>%
    filter(chain == "chain1") %>%
    pivot_longer(!c(chain, iteration), 
                 names_to = c("mother_cell_id",".value"), 
                 names_pattern = "(.*)_(.)") %>% 
    mutate(across(everything(), ~replace_na(.x, 0))) %>% 
    count(mother_cell_id, `1`, `2`) %>% 
    group_by(mother_cell_id) %>% 
    mutate(p = n/sum(n)) %>% 
    filter(n == max(n)) %>% 
    select(-n) %>% 
    pivot_longer(c(`1`, `2`), values_to = "n_episomes") %>% 
    summarise(X1 = max(n_episomes), X2 = ifelse(n() > 1, min(n_episomes), 0), total = sum(n_episomes), p = unique(p)) %>%
    mutate(pair = paste0("(", X1, ",", X2, ")")) 
  
  print(temp_df %>% count(X1 == X2))
  
  pair_levels <- temp_df %>% count(pair) %>% arrange(desc(n)) %>% pull(pair)
  if(length(pair_levels) > 10){
    pair_labels <- c(pair_levels[1:10], rep("other", length(pair_levels)-10))  
  }else{
    pair_labels <- pair_levels
  }
  
  intensity_scatter <- daughter_cell_data %>% 
    group_by(mother_cell_id, cell_id, daughter_cell) %>% 
    summarise(total_intensity = sum(total_cluster_intensity)) %>% 
    group_by(mother_cell_id) %>% 
    summarise(daughter_1 = max(total_intensity), daughter_2= ifelse(n() > 1, min(total_intensity), 0)) %>% 
    ggplot(aes(daughter_1, daughter_2)) + 
    geom_abline(color = "gray", lty = "dashed") + 
    geom_point(size = 3) + 
    labs(x = "Total dot intensity of daughter cell with greater intensity",
         y = "Total dot intensity of daughter cell with less intensity") + 
    tune::coord_obs_pred() 
  
  
  intensity_distribution <- mother_cell_data %>%
    group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>%
    ggplot(aes(total_cluster_intensity, y = after_stat(count / sum(count)))) +
    geom_histogram(aes(fill = "non-dividing cells", color = "non-dividing cells"), alpha = 0.5) +
    geom_histogram(data = daughter_cell_data %>% group_by(mother_cell_id) %>%
                     mutate(int = sum(total_cluster_intensity)),
                   aes(int, fill = "daughter cell pairs", color = "daughter cell pairs"), alpha = 0.5) +
    labs(x = "cluster intensity", fill = "cell set", y= "frequency") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = c(1,1), legend.justification = c(1,0.9),
          legend.background = element_blank(), legend.title = element_blank())  +
    guides(color = "none") +
    labs(x = "Total dot intensity") +
    scale_color_manual(values = safe_colorblind_palette) +
    scale_fill_manual(values = safe_colorblind_palette) 
  
  
  pair_histogram <- temp_df  %>%
    count(pair, X1, X2) %>%
    arrange(desc(n)) %>%
    slice(1:10) %>%
    arrange(desc(n), X1, X2) %>% 
    mutate(pair = fct_inorder(pair)) %>% 
    ggplot(aes(pair, n)) +
    geom_bar(stat = "identity") +
    labs(x = "Estimated number of of episomes in daughter cell pairs",
         y = "Number of occurrences") + 
    # scale_fill_manual(values = rev(safe_colorblind_palette[2:11])) +
    theme(legend.position = "none")
  
  figure <- intensity_distribution +
    intensity_scatter + pair_histogram + 
    plot_layout(ncol = 1, heights = c(1,2,1)) 
  
  # ggsave(here(results_folder, "figure.png"), figure, width = 5, height = 10)
  
  
}

