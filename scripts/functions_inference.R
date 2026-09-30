## This script has functions required to calculate the likelihood, estimate its parameters, and calculate uncertainty

#####
# Likelihood function
#####

# The likelihood itself depends on initial number of episomes, episomes in each daughter cell, and the probability of replication and division
likelihood <- function(X1, X2, X0, Pr, Ps){
  R = X1 + X2 - X0
  # Portion of likelihood that isn't in summation
  if(R > X0 | R < 0){
    L = 0
  }else{
    L = choose(X0, R)*Pr^R*(1-Pr)^(X0-R)*2^(X1 != X2)  
    # Sum over possible values of s (number of replicated episome pairs that segregate) and k (number of replicated episome pairs that don't segregate):
    sums = 0
    for(s in 0:R){
      for(k in 0:((X1-s)/2)){
        sums = sums + 0.5^(X0-s)*choose(R, s)*Ps^s*(1-Ps)^(R-s)*choose(R-s, k)*choose(X0-R, X1-s- 2*k)
      }
    }
    L = L*sums
  }
  return(c("likelihood" = L))
}

# Replication-only likelihood: the probability of observing daughter counts
# X1 and X2 from a mother with X0 episomes, given replication efficiency Pr,
# ignoring segregation.
#
# The number of successful replication events is R = X1 + X2 - X0, which is
# binomial with X0 trials and success probability Pr. Used by the just_Pr
# branch of the pipeline, where segregation is assumed perfect.
#
# Arguments:
#   X1, X2   Episome counts in the two daughter cells.
#   X0       Episome count in the mother cell before replication.
#   Pr       Replication efficiency, in [0, 1].
#   Ps       Ignored. Present so this can be swapped for likelihood() without
#            changing the call site.
#
# Returns:
#   A named numeric of length 1, "likelihood".
likelihood_Pr <- function(X1, X2, X0, Pr, Ps = NA){
  R = X1 + X2 - X0
  # Portion of likelihood that isn't in summation
  L = dbinom(R, X0, Pr)
  return(c("likelihood" = L))
}

#####
# Grid search functions
#####

# Function to find the parameters Pr and Ps that maximize the likelihood of the observed data
calculate_maximum_likelihood <- function(data, Pr_values, Ps_values){
  # data will have a column for id, X0, X1, and X2
  
  # Function to calculate likelihood for a given observed data
  get_likelihood <- function(data, Pr_values, Ps_values){
    
    
    # Set up parameters (Pr, Ps, and X0) to test data with
    parameter_grid <- expand_grid(Pr = Pr_values, Ps = Ps_values)
    
    # Set up a dataframe of likelihoods to calculate
    # Note that we only calculate the likelihood of the observed data once given each Pr and Ps since all observations in this "chunk" are the same
    parameter_grid <- expand_grid(parameter_grid, data %>% distinct(X0, X1, X2))
    
    likelihoods <- parameter_grid %>% 
      # Calculate likelihood at each parameter combination
      bind_cols(pmap_df(., likelihood)) %>% 
      group_by(Pr, Ps) %>%
      summarise(log_likelihood = log(sum(likelihood))*nrow(data)) %>%
      ungroup
    likelihoods
  }
  
  # Nest the data frame to observations with the same values
  data %>% 
    # make sure X1 is larger number of episomes
    mutate(temp_X1 = ifelse(X1 > X2, X1, X2), temp_X2 = ifelse(X1 > X2, X2, X1)) %>% 
    select(-X1, -X2) %>% 
    rename(X1 = temp_X1, X2 = temp_X2) %>% 
    mutate(outcome = paste(X0, X1, X2, sep = "_")) %>% 
    nest(data = c(X0, X1, X2, id)) %>% 
    # For each set of observed data, calculate the likelihood at all combinations of Pr and Ps
    mutate(l = map(data, get_likelihood, Pr_values, Ps_values)) %>%  pull(l) %>% 
    # Pull together the results for each observed value
    bind_rows %>% group_by(Pr, Ps) %>% 
    # Add up log likelihoods to get likelihood of observing all the data given Pr and Ps
    summarise(log_likelihood = sum(log_likelihood))%>% 
    ungroup()
  
}


# Grid search for the replication (Pr) and segregation (Ps) efficiency that
# maximize the likelihood of the observed daughter-cell counts, marginalizing
# over the unobserved mother-cell episome count X0.
#
# X0 is never observed directly, so for each daughter pair the likelihood is
# summed over every X0 consistent with the data, weighted by a prior. The
# support runs from ceiling((X1 + X2) / 2) -- the fewest mother episomes that
# could yield the observed daughters, achieved when every episome replicates --
# up to X1 + X2, which is the case where none did. The prior over that support
# is Poisson(lambda), renormalized to sum to 1.
#
# Observations sharing the same (X1, X2) are collapsed and the likelihood is
# computed once per distinct pair, then weighted by frequency. This is a
# substantial speed-up for the low-count data here, where many pairs repeat.
#
# Arguments:
#   data        Data frame with columns id, X1, X2 -- daughter-cell episome
#               counts. X0 is not required and is ignored if present.
#   Pr_values   Vector of replication efficiencies to evaluate, typically
#               seq(0, 1, by = increment).
#   Ps_values   Vector of segregation efficiencies, same convention.
#   lambda      Rate of the Poisson prior over X0.
#   just_Pr     If TRUE, estimate replication efficiency alone, assuming perfect
#               segregation; Ps is set to NA and the Ps grid is not searched.
#               Default FALSE.
#
# Returns:
#   A data frame with one row per (Pr, Ps) combination and a log_likelihood
#   column. Pass to calculate_CI() for the confidence region.
calculate_maximum_likelihood_unknownX0 <- function(data, Pr_values, Ps_values, lambda, just_Pr = F){
  # data will have a column for id, X1, and X2
  
  # Function to calculate likelihood for a given observed data
  get_likelihood <- function(data, Pr_values, Ps_values, lambda){
    
    # Build PMF for X0:
    max_X0 <- unique(data$X1 + data$X2)
    min_X0 <- ceiling(max_X0/2)
    PMF <- tibble(X0 = min_X0:max_X0, prob = dpois(X0, lambda)) %>% 
      # normalize probability to sum to 1
      mutate(prob = prob/sum(prob))
    
    # Set up parameters (Pr, Ps, and X0) to test data with
    parameter_grid <- expand_grid(Pr = Pr_values, Ps = Ps_values, X0 = PMF$X0)
    
    # Set up a dataframe of likelihoods to calculate
    # Note that we only calculate the likelihood of the observed data once given each Pr and Ps since all observations in this "chunk" are the same
    parameter_grid <- expand_grid(parameter_grid, data %>% distinct(X1, X2))
    
    if(!just_Pr){
      likelihoods <- parameter_grid %>% 
        # Calculate likelihood at each parameter combination
        bind_cols(pmap_df(., likelihood)) %>% 
        merge(PMF) %>% 
        mutate(likelihood = likelihood*prob) %>% 
        group_by(Pr, Ps) %>%
        summarise(log_likelihood = log(sum(likelihood))*nrow(data)) %>%
        ungroup
    }else{
      likelihoods <- parameter_grid %>% 
        mutate(Ps = NA) %>% 
        distinct() %>% 
        # Calculate likelihood at each parameter combination
        bind_cols(pmap_df(., likelihood_Pr)) %>% 
        merge(PMF) %>% 
        mutate(likelihood = likelihood*prob) %>% 
        group_by(Pr, Ps) %>%
        summarise(log_likelihood = log(sum(likelihood))*nrow(data)) %>%
        ungroup
    }
    likelihoods
  }
  
  # Nest the data frame to observations with the same values
  data %>% 
    # make sure X1 is larger number of episomes
    mutate(temp_X1 = ifelse(X1 > X2, X1, X2), temp_X2 = ifelse(X1 > X2, X2, X1)) %>% 
    select(-X1, -X2) %>% 
    rename(X1 = temp_X1, X2 = temp_X2) %>% 
    mutate(outcome = paste(X1, X2, sep = "_")) %>% 
    nest(data = c(X1, X2, id)) %>% 
    # For each set of observed data, calculate the likelihood at all combinations of Pr and Ps
    mutate(l = map(data, get_likelihood, Pr_values, Ps_values, lambda)) %>%  pull(l) %>% 
    # Pull together the results for each observed value
    bind_rows %>% group_by(Pr, Ps) %>% 
    # Add up log likelihoods to get likelihood of observing all the data given Pr and Ps
    summarise(log_likelihood = sum(log_likelihood))%>% 
    ungroup()
  
}


# Computes the 95% confidence region for the replication (Pr) and segregation
# (Ps) efficiency from a likelihood surface produced by
# calculate_maximum_likelihood() or calculate_maximum_likelihood_unknownX0().
#
# Log likelihoods are shifted by their maximum, exponentiated, and normalized to
# sum to 1 across the grid, giving a probability for each (Pr, Ps) box. Boxes
# are then ranked by probability and accumulated from the most probable
# downward until the cumulative mass reaches 0.95; the region is the set of
# boxes included. This is a highest-density region over the normalized
# likelihood, equivalent to a Bayesian credible region under a flat prior on the
# grid, rather than a likelihood-ratio confidence region 
#
# Arguments:
#   likelihoods  Data frame with one row per (Pr, Ps) combination and a
#                log_likelihood column, as returned by the
#                calculate_maximum_likelihood* functions.
#   marginal_CI  If FALSE (default), return the joint region: accumulate over
#                the two-dimensional grid and take the range of Pr and of Ps
#                among the included boxes. 
#                If TRUE, sum the probability over one parameter at a time and
#                accumulate within each margin separately, giving intervals that
#                ignore the correlation between Pr and Ps.
#
# Returns:
#   A list of two elements:
#     probs      The ranked probability table. In the joint case this is
#                filtered to the boxes inside the region; in the marginal case
#                it is the full grid, unfiltered.
#     estimates  A one-row data frame with MLE_Pr, MLE_Ps, min_Pr, min_Ps,
#                max_Pr, max_Ps. The MLE is the single most probable box, so its
#                resolution is limited by the grid increment.
#
calculate_CI <- function(likelihoods, marginal_CI = F){
  probabilities <- likelihoods %>% 
    # Convert to probability
    mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
           probability = likelihood/sum(likelihood)) %>% 
    # Rank by probability 
    arrange(desc(probability)) 
  
  
  if(marginal_CI){
    CI_Pr <- probabilities %>% 
      group_by(Pr) %>%
      summarise(probability = sum(probability)) %>% 
      arrange(desc(probability)) %>% 
      mutate(cum_sum = cumsum(probability)) %>% 
      filter(cum_sum <= 0.95) %>% 
      pull(Pr) %>% range
    
    CI_Ps <- probabilities %>% 
      group_by(Ps) %>%
      summarise(probability = sum(probability)) %>% 
      arrange(desc(probability)) %>% 
      mutate(cum_sum = cumsum(probability)) %>% 
      filter(cum_sum <= 0.95) %>% 
      pull(Ps) %>% range
    
  }else{
    probabilities <- probabilities %>% 
      # Find boxes that sum to 0.95
      mutate(cumulative_prob = cumsum(probability),
             t = abs(cumulative_prob - 0.95))
    
    threshold <- probabilities %>% filter(cumulative_prob >= 0.95) %>% 
      filter(t == min(t)) %>% 
      pull(cumulative_prob)
    probabilities <- probabilities %>% filter(cumulative_prob <= threshold) 
    
    # 95% CI for Pr:
    CI_Pr <- range(probabilities$Pr)
    # 95% CI for Ps:
    CI_Ps <- range(probabilities$Ps)
  }
  
   
  # MLE:
  mle <- probabilities %>% slice(1) %>% select(Pr, Ps) %>% unlist
  
  list(probs = probabilities, estimates = data.frame(MLE_Pr = mle[1], MLE_Ps = mle[2], min_Pr = CI_Pr[1], min_Ps = CI_Ps[1], max_Pr = CI_Pr[2], max_Ps =CI_Ps[2]))
}

# Wrapper that runs the Pr/Ps grid search, computes a confidence region, and
# optionally plots the result. The single entry point used by both the
# benchmarking scripts and run_pipeline().
#
# Arguments:
#   simulated_data  Data frame with columns id, X1, X2 and (when known_X0 is
#                   TRUE) X0. Named for its use with synthetic data from
#                   simulate_multiple_cells(), but real data in the same shape
#                   works identically.
#   viz             If TRUE (default), print a plot_grid_search() figure.
#   increment       Spacing of the Pr and Ps grid over [0, 1]. Default 0.01, so
#                   101 x 101 combinations. Cost scales with the square, so
#                   coarsen this first when a search is too slow.
#   known_X0        If TRUE (default), use the mother-cell count in the data
#                   directly. If FALSE, marginalize over X0 using lambda.
#   lambda          Rate of the Poisson prior over X0. Required when
#                   known_X0 = FALSE, ignored otherwise.
#   CI              If TRUE (default), compute the confidence region via
#                   calculate_CI() and include it in the return value.
#   just_Pr         If TRUE, estimate replication efficiency only. Default FALSE.
#   marginal_CI     If TRUE, return marginal intervals for Pr and Ps
#                   separately rather than the joint region. Default FALSE.
#                   The joint region is what the manuscript reports.
#
# Returns:
#   A list with grid_search (the likelihood surface) and simulated_data; plus
#   estimates and top_95 when CI = TRUE.
run_grid_search <- function(simulated_data, viz = T, increment = 0.01, known_X0 = T, lambda = NA, CI = T, just_Pr = F, marginal_CI = F){
  # simulated_data is a data frame simulated with the columns X0, X1, X2, and id (outcome of simulate_multiple_cells)
  # viz is a logical indicating if the results should be visualized
  # increment is a parameter that determines how fine-grained the grid search is
  
  parameter_values <- seq(0,1,by = increment)
  if(known_X0){
    grid_search <- calculate_maximum_likelihood(simulated_data, parameter_values, parameter_values)  
  }else{
    grid_search <- calculate_maximum_likelihood_unknownX0(simulated_data, parameter_values, parameter_values, lambda, just_Pr = just_Pr)
  }
  
  out <- list(grid_search = grid_search, simulated_data = simulated_data)
  
  if(CI){
    CIs <- calculate_CI(grid_search, marginal_CI)$estimates
    top_95 <- calculate_CI(grid_search, marginal_CI)$probs   
    out <- list(grid_search = grid_search, estimates = CIs, top_95 = top_95, simulated_data = simulated_data)
  }
  
  if(viz){
    print(plot_grid_search(out, known_X0))
  }
  
  # out <- list(grid_search = grid_search, simulated_data = simulated_data)
  return(out)
}

# Function to plot the output of run_grid_search()
plot_grid_search <- function(run_grid_search_out, simulation = T, prob = F, error_bars = T, poster = F){
  grid_search <- run_grid_search_out$grid_search
  CIs <- run_grid_search_out$estimates
  top_95 <- run_grid_search_out$top_95
  simulated_data <- run_grid_search_out$simulated_data
  
  fill_label <- ifelse(prob, "probability", "Log Likelihood")
  
  if(prob){
    if(!poster){
      grid_search_plot <- grid_search %>% 
        mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
               probability = likelihood/sum(likelihood)) %>%  
        ggplot(aes(Pr, Ps)) + 
        geom_raster(aes(fill = probability)) 
    }else{
      grid_search_plot <- grid_search %>% 
        mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
               probability = likelihood/sum(likelihood)) %>%  
        ggplot(aes(Pr*100, Ps*100)) + 
        geom_raster(aes(fill = probability)) 
    }
  }else{
    grid_search_plot <- grid_search %>% 
      ggplot(aes(Pr, Ps)) + 
      geom_raster(aes(fill = log_likelihood)) 
  }
  
  if(!poster){
    Pr_confidence_boundary <- top_95 %>% group_by(Ps) %>% summarise(min_pr = min(Pr), max_pr = max(Pr)) %>% pivot_longer(contains("pr"), values_to = "Pr") %>% select(-name)
    Ps_confidence_boundary <- top_95 %>% group_by(Pr) %>% summarise(min_ps = min(Ps), max_ps = max(Ps)) %>% pivot_longer(contains("ps"), values_to = "Ps") %>% select(-name)
    confidence_boundary <- rbind(Pr_confidence_boundary, Ps_confidence_boundary)
    
    grid_search_plot <- grid_search_plot +
      geom_raster(data = confidence_boundary, fill = "white", alpha = 0.7) +
      # geom_point(data = CIs, aes(MLE_Pr, MLE_Ps, color = "MLE")) +
      # scale_color_manual(values = c("black", "red")) + 
      labs(fill = fill_label, 
           caption = "Error Bars show 95% CI for Pr and Ps", 
           title = "Grid Search for Maximum Likelihood Estimate of Pr and Ps",
           x = "Replication Efficiency", y= "Segregation Efficiency")  + 
      theme_classic() + 
      scale_fill_viridis_c(option = "magma")
  }else{
    
    grid_search_plot <- grid_search_plot +
      # geom_point(data = CIs, aes(MLE_Pr, MLE_Ps, color = "MLE")) +
      # scale_color_manual(values = c("black", "red")) + 
      labs(fill = fill_label, 
           title = "Grid Search for Maximum Likelihood Estimate of Pr and Ps",
           x = "Replication Efficiency (%)", y= "Segregation Efficiency (%)")  + 
      theme_classic() + 
      scale_fill_viridis_c(option = "magma") 
  }
  
  if(error_bars){
    grid_search_plot <- grid_search_plot + 
      geom_errorbarh(data = CIs, aes(y = MLE_Ps, xmin = min_Pr, xmax = max_Pr), height = 0, inherit.aes = F, color = "white", alpha= 0.7) +
      geom_errorbar(data = CIs, aes(x = MLE_Pr, ymin = min_Ps, ymax = max_Ps), width = 0, inherit.aes = F,  color =  "white", alpha = 0.7) 
  }
  
  if(simulation){
    grid_search_plot <- grid_search_plot + geom_point(data = distinct(simulated_data, Pr, Ps),  aes(color = "Parameter Values")) + # Add actual parameter values
      scale_color_manual(values = c("red")) + 
      labs(color = "")
  }else{
    grid_search_plot <- grid_search_plot + guides(color = "none")
  }
  
  Pr_marginal_likelihood <- grid_search %>% 
    mutate(likelihood = exp(log_likelihood - max(log_likelihood))) %>% 
    group_by(Pr) %>% summarise(likelihood = sum(likelihood)) %>% 
    ggplot(aes(Pr, likelihood)) + geom_line() + labs(y = "Marginal\nLikelihood") + 
    theme_classic()
  
  Ps_marginal_likelihood <- grid_search %>% 
    mutate(likelihood = exp(log_likelihood - max(log_likelihood))) %>% 
    group_by(Ps) %>% summarise(likelihood = sum(likelihood)) %>% 
    ggplot(aes(Ps, likelihood)) + geom_line() + labs(y = "Marginal\nLikelihood") + 
    theme_classic()
  
  if(poster){
    grid_search_plot + 
      labs(caption = "") + 
      theme(
        # panel.background = element_rect(fill='black'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent', color = "transparent"), #transparent legend panel
        plot.title = element_blank()
      ) + 
      geom_hline(yintercept = 50, lty = "dashed", color = "red") + 
      geom_vline(xintercept = 50, lty = "dashed", color = "red")
    # plot_layout(guides = "collect", heights = c(0.5,2), widths = c(2,0.5)) + 
    # plot_annotation(title = "Grid Search for Maxmimum Likelihood Estimate of Pr and Ps") 
  }else{
    Pr_marginal_likelihood + 
      theme(axis.title.x = element_blank(), plot.margin = margin(0,0,0,0), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
      plot_spacer() + 
      grid_search_plot + theme(plot.margin = margin(0,0,0,0), plot.title = element_blank()) +
      Ps_marginal_likelihood + coord_flip() + 
      theme(axis.title.y = element_blank(), plot.margin = margin(0,0,0,0),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
      plot_layout(guides = "collect", heights = c(0.5,2), widths = c(2,0.5)) + 
      plot_annotation(title = "Grid Search for Maxmimum Likelihood Estimate of Pr and Ps")
  }
  
  
}


#####
# Gibbs Sampling Functions
#####


# function to get the un-normalized probability of n
log_likelihood_n <- function(n, mu, sigma2, I, n_prior){ 
  prior_fam <- n_prior[[1]]
  prior_param <- n_prior[[2]]
  prior_prob <- ifelse(prior_fam == "pois", 
                       dpois(n, prior_param)/sum(dpois(1:100, prior_param)),
                       dgeom(n, prior_param)/sum(dgeom(1:100, prior_param))
  )
  log(dnorm(I, n*mu, sqrt(n*sigma2))) + log(prior_prob)
}

# Runs Gibbs sampling to infer the number of episomes per LANA cluster (n), the
# mean fluorescence intensity of a single episome (mu), and the precision of
# that intensity (tau), from observed cluster intensities.
#
# Each iteration cycles through three conditional draws: n for every cluster in
# turn, sampled from a categorical distribution over 1 to 100 episomes built
# from log_likelihood_n(); then mu from its normal conditional; then tau from
# its gamma conditional. Cluster likelihoods are computed on the log scale and
# shifted by their maximum before exponentiating, to avoid underflow.
#
# Arguments:
#   tau0          Initial value for tau, the precision (inverse variance) of
#                 intensity per episome. run_pipeline() uses 1e-5, chosen so the
#                 implied variance is wide enough for the chain to converge from
#                 a poor starting point.
#   mu0           Initial value for mu, the intensity of a single episome.
#                 run_pipeline() draws this from the range of observed
#                 intensities.
#   I             Named numeric vector of observed total cluster intensities,
#                 one element per cluster, names being cluster_id.
#   n_iterations  Number of iterations to run. No burn-in is discarded here --
#                 the caller does that.
#   ns            Optional named vector of starting episome counts per cluster.
#                 Defaults to NA, in which case starting values are round(I/mu0)
#                 with zeros raised to 1. Names must match those of I; the
#                 vector is reordered to match.
#   n_prior       Prior on the number of episomes per cluster, as
#                 list(distribution, parameter) -- e.g. list("pois", 1). Passed
#                 through to log_likelihood_n(). See
#                 benchmark_prior_sensitivity.R for the alternatives tested.
#
# Returns:
#   A data frame with n_iterations rows and columns iteration, mu, tau, and one
#   column per cluster holding that cluster's sampled episome count. Bind
#   multiple chains together and pass to convergence_results() for diagnostics.
run_gibbs <- function(tau0, mu0, I, n_iterations, ns = NA, n_prior){
  
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
      nk_like <- nks %>% sapply(log_likelihood_n, mu = mu[j-1], sigma2 = 1/tau[j-1], I = I[k], n_prior) # log likelihood
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


#####
# Functions to assess convergence of markov chains
#####

# To assess the convergence of the chains, I calculated $\hat{R}$, bulk effective sample size (ESS), and tail effective 
# sample size based on this [vignette](https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html). 
#  $\hat{R}$ measures the ratio of the average variance of draws within each chain to the variance of the pooled draws across chains; 
# if all chains are at equilibrium, these will be the same and $\hat{R}$ will be one. The effective sample size is an estimate of the 
# number of independent draws from the posterior distribution of the estimand of interest. The larger the ratio of $n_{eff}$ to $N$, 
# the better. Bulk-ESS is related to efficiency of mean and median estimates and tail-ESS is related to the efficiency of the variance and 
# quantile estimates of the posterior distribution.

get_convergence_stats <- function(data){
  sims <- data %>%
    pivot_wider(names_from = chain, values_from = value) %>%
    select(-iteration) %>% as.matrix
  
  return(tibble(Rhat = Rhat(sims), ESS_bulk = ess_bulk(sims), ESS_tail = ess_tail(sims)))
}


# Computes MCMC convergence diagnostics for every parameter in a set of chains.
#
# Reshapes the chains to long form, groups by parameter, and applies
# get_convergence_stats() to each. Parallelized with future_map, so set a
# future plan beforehand if you want it to actually run in parallel.
#
# Arguments:
#   all_chains   Data frame of posterior draws with columns chain, iteration,
#                and one further column per parameter.
#
# Returns:
#   One row per parameter, with the diagnostics from get_convergence_stats()
#   (Rhat, bulk ESS, tail ESS).
convergence_results <- function(all_chains){
  convergence <- all_chains %>% 
    pivot_longer(!c(chain, iteration)) %>% 
    group_by(name) %>% nest() %>% 
    mutate(data = future_map(data, ~get_convergence_stats(.x))) %>% 
    unnest(cols = c(data)) %>% 
    ungroup
}




