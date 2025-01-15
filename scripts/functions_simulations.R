#####
# Simulation functions
#####

## Color palette for simulations

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#####
# Functions for Simulation
#####

# Simulates the fate of episomes during cell divison
makeChildren<- function(pRep, pSeg, numEpisomes){
  if(numEpisomes == 0){
    return(c(0, 0))
  }
  else{
    outcomes = c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2)
    refMat = matrix(outcomes, nrow= 5, byrow=T) #possible episome fates
    # probVec = c(1/2*(1-pRep), 1/2*(1-pRep), pRep*(pSeg+1/2*(1-pSeg)), pRep*1/4*(1-pSeg), pRep*1/4*(1-pSeg)) # probabilities when pSeg = probability of tethering
    probVec = c(1/2*(1-pRep), 1/2*(1-pRep), pRep*pSeg, 1/2*pRep*(1-pSeg), 1/2*pRep*(1-pSeg)) # probabilities when pSeg = probability of segregation regardless of mechanism
    result = rmultinom(n=1, prob=probVec, size=numEpisomes) #simulate each independently and add the results in the return statement
    return(t(result)%*%refMat)
  }
}

# Simulates one step of the cell population dynamics (either cell birth or death) and advances the time
simStepFlex <- function(pRep, pSeg, cells, birthVec, deathVec, selectAgainstZero = T, max_epi){
  #all FUN arguments are functions
  #cells is a compressed vector of the number of cells with 0, 1, 2, ... episomes
  if(sum(cells) == 0){
    return(cells)
  }
  else{
    type <- sample(length(cells), size=1, replace=TRUE, prob = cells*(birthVec+deathVec))
    b = birthVec[type]
    d = deathVec[type]
    timeAdvance=rexp(n=1, rate=sum(cells*(birthVec+deathVec)))
    test = runif(n=1)
    if(test < d/(b+d)){
      #print("DEATH")
      cells[type] = cells[type] - 1
      if(selectAgainstZero == TRUE){
        cells[1] <- 0
      }
      return(list(timeAdvance, cells))
    }
    else{
      #print("BIRTH")
      children <- makeChildren(pRep, pSeg, type-1)
      if(children[1] > max_epi){
        #print("OVERFLOW")
        children[1] = max_epi
      }
      if(children[2] > max_epi){
        #print("OVERFLOW")
        children[2] = max_epi
      }
      #print(children)
      #print(type)
      cells[type] = cells[type]-1
      cells[children[1]+1] = cells[children[1]+1] + 1
      cells[children[2]+1] = cells[children[2]+1] + 1
      if (selectAgainstZero==TRUE){
        cells[1] <- 0
      }
      return(list(timeAdvance, cells))
    }
  }
}

# Simulates multiple populations (nTrials) of cells with a constant size that may or may not be under selection (selectAgainstZero)
extinction <- function(pRep, pSeg, nTrials, n_epi, selectAgainstZero = F, n_cells = 1000, n_cells_start = NULL,
                       d = 1, b = 3, stop_time = NULL, growth_advantage = NULL){
  
  results = 1:nTrials*0
  times = c()
  totals = c()
  trials = c()
  zero = c()
  one = c()
  two = c()
  three = c()
  four = c()
  five = c()
  six = c()
  seven = c()
  eight = c()
  nine = c()
  
  if(is.null(n_cells_start)) n_cells_start <- n_cells
  
  i = 1
  j = 1
  for(z in 1:nTrials){
    print(z)
    time = 0
    deathVec = rep(d, 10)
    cells = rep(0, 10)
    cells[n_epi + 1] <- n_cells_start
    total = sum((0:9)*cells)
    indicator <- T
    while(indicator){
      if(time > 700 & pRep == 1 & pSeg == 1) break
      birthVec = rep((b-d)*(1-sum(cells)/n_cells) + d, 10)
      if(!is.null(growth_advantage)) birthVec = birthVec*growth_advantage
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selectAgainstZero = selectAgainstZero, max_epi = 9)
      cells = result[[2]]
      time=time + result[[1]]
      #print(cells)
      total = sum((0:9)*cells)
      if(j %% 100 == 0 | j == 1){ # report out every 100 iterations
        times[i] <- time
        totals[i] <- total
        trials[i] <- z
        zero[i] <- cells[1]
        one[i] <- cells[2]
        two[i] <- cells[3]
        three[i] <- cells[4]
        four[i] <- cells[5]
        five[i] <- cells[6]
        six[i] <- cells[7]
        seven[i] <- cells[8]
        eight[i] <- cells[9]
        nine[i] <- cells[10]
        i <- i + 1
      }
      j <- j + 1
      
      # if end time provided, stop then. Otherwise:
      # if simulating with selection, run for 700 generations
      # if simulating without selection, run until there are no more episomes
      indicator <- ifelse(!is.null(stop_time), time <= stop_time, 
                          ifelse(selectAgainstZero, time <= 700 & total > 0, total > 0))
      # indicator <- total > 0
    }
    results[z]= time
    
    # Add last timepoint
    times[i] <- time
    totals[i] <- total
    trials[i] <- z
    zero[i] <- cells[1]
    one[i] <- cells[2]
    two[i] <- cells[3]
    three[i] <- cells[4]
    four[i] <- cells[5]
    five[i] <- cells[6]
    six[i] <- cells[7]
    seven[i] <- cells[8]
    eight[i] <- cells[9]
    nine[i] <- cells[10]
    
    i <- i+1
    
  }
  total_df <- data.frame(trial = trials, time = times, total = totals,
                         zero, one, two, three, four, five, six, seven, eight, nine)
  return(list(ExtinctionTime = tibble(ExtinctionTime = results), Totals = total_df))
}

# Simulates multiple populations (nRuns) of exponentially growing cells that may or may not be under selection (selectAgainstZero)
exponential_growth <- function(pRep, pSeg, nIts, nRuns, n_cells_start = 1, n_epi_start = 3, selection, max_epi, 
                               stop_size = NULL, d = 0, b = 1, growth_advantage = NULL, initial_conditions = NULL, 
                               start_times = 0, stop_time = NULL){
  
  #try using matrix first and then convert to data frame
  # rm(data)
  dataCUT = 1
  start = 1:100
  mid = 11:100*10
  end = seq(2*1000, 1e5, by = 1000)
  recordList = c(start, mid, end)
  if(stop_size > max(end)){
    recordList = c(recordList, seq(1e5, 1.25e5, by = 2000), seq(1.25e5, min(1e6, stop_size), by = 5000))
    if(stop_size > 1e6) recordList = c(recordList, seq(1e6, stop_size, by = 10000))
  }
  # data <- matrix(NA, nrow=500*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5)
  data <- matrix(NA, nrow=100*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5)
  # data <- matrix(NA, nrow=nIts*nRuns, ncol=5)
  print(dim(data))
  names(data) = c("run", "time", "episomes", "frac", "total")
  z = 1
  for(run in 1:nRuns){
    cat("z:", z, "\n")
    cat("run", run, "\n")
    deathVec = rep(d, max_epi + 1)
    if(is.null(initial_conditions)){
      cells = rep(0, max_epi + 1)
      # Start with n_cells_start cells each with n_epi_start episomes
      cells[n_epi_start+1] <- n_cells_start
    }else{
      cells = unname(initial_conditions[run,])
    }
    if(length(start_times) == 1) start_times = rep(start_times, nRuns)
    times = rep(start_times[run], nIts + 1)
    totalEps = append(c(sum((0:max_epi)*cells)), rep(0, nIts))
    for(j in 1:length(cells)){
      tmp = c(run, times[1], j-1, cells[j]/sum(cells), sum(cells))
      try({data[z,] <- tmp})
      z=z+1
    }
    
    tmp = c(run, times[1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells))
    try({data[z,] <- tmp})
    z= z+1
    #data <- rbind(data, data.frame(run=run, time=times[1], episomes=-1, count = sum((0:9)*cells)))
    for(i in 1:nIts){
      
      birthVec = rep(b, max_epi + 1)
      if(!is.null(growth_advantage)) birthVec = birthVec*c(1, rep(1 + growth_advantage, max_epi))
      
      cells2 = cells
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selection, max_epi)
      cells = round(result[[2]])
      times[i+1] = times[i] + result[[1]]
      
      if(any(cells*(birthVec + deathVec) < 0)){
        print(cells)
        print(birthVec)
        print(deathVec)
        
        cat("\nprior:::\n")
        
        print(cells2)
      }
      
      #print(cells)
      totalEps[i+1] = sum((0:max_epi)*cells)
      # if(i%%1000 == 0){
      #   print(i)
      # }
      if(sum(cells) %in% recordList){
        #print(i)
        for(j in 1:length(cells)){
          tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
          try({data[z,] <- tmp})
          z = z+1
        }
        tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
        try({data[z,] <- tmp})
        z= z+1
        #data <- rbind(data, data.frame(run=run, time=times[i+1], episomes=-1, count = sum((0:9)*cells)))
        
      }
      
      # print(sum(cells))
      
      if(sum(cells) == 0){ # add final record and break
        for(j in 1:length(cells)){
          tmp = c(run, times[i+1], j-1, 0, sum(cells)) #added sum cells
          try({data[z,] <- tmp})
          z = z+1
        }
        tmp = c(run, times[i+1], -1, 0, sum(cells)) #totalEps[i+1])
        try({data[z,] <- tmp})
        break
      } 
      
      if(!is.null(stop_time)){
        if(times[i+1] >= stop_time + start_times[run]){
          for(j in 1:length(cells)){
            tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
            try({data[z,] <- tmp})
            z = z+1
          }
          tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
          try({data[z,] <- tmp})
          break
        } 
      }
      
      if(!is.null(stop_size)){
        if(sum(cells) >= stop_size){
          for(j in 1:length(cells)){
            tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
            try({data[z,] <- tmp})
            z = z+1
          }
          tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
          try({data[z,] <- tmp})
          break
        } 
      }
      
      
      # End if we run out of space in data
      if(any(!is.na(data[nrow(data),]))) break
      
      
    }
    
    #plot(times, totalEps)
  }
  # Remove NAs
  data <- data[complete.cases(data), ]
  data <- data.frame(data)
  names(data) <- c("run", "time", "episomes", "frac", "total")
  return(data)
}


# Simulates multiple KSHV-dependent tumors (nRuns) growing from one cell, with a reduction in segregation and/or replication efficiency once
# the tumor reaches a defined size
PEL_simulations <- function(pRep, pSeg, pRep_reduced = NULL, pSeg_reduced = NULL, nRuns, n_cells_start = 1, n_epi_start = 3, selection, max_epi, 
                            stop_size = 1.5e5, treatment_size = 1e5, d = 0, b = 1, growth_advantage = NULL, add_to = NULL, stop_time = NULL,
                            keep_extinct = FALSE, save_baseline_simulations = NULL, save_treatment_simulations = NULL){
  
  if(is.null(add_to)){
  
      uncontrolled_growth <- exponential_growth(pRep, pSeg, 1e8, nRuns, 
                                                n_cells_start = n_cells_start, n_epi_start = n_epi_start, 
                                                selection, max_epi, stop_size = treatment_size, 
                                                d = d, b = b, growth_advantage = growth_advantage)
      
      extinct_runs <- uncontrolled_growth %>% group_by(run) %>% filter(max(total) < treatment_size) %>% pull(run) %>% unique
      
      print(extinct_runs)
      
      if(keep_extinct){
        extinct_run_dict <- nRuns+1:length(extinct_runs) %>% setNames(extinct_runs)
        extinct_df <- uncontrolled_growth %>% filter(run %in% extinct_runs) %>% 
          mutate(pRep = pRep, pSeg = pSeg, run = recode(run, !!!extinct_run_dict))
      }
      
      while(length(extinct_runs) > 0){
        
        nRuns2 <- length(extinct_runs) 
        
        uncontrolled_growth2 <- exponential_growth(pRep, pSeg, 1e8, nRuns2, 
                                                   n_cells_start = n_cells_start, n_epi_start = n_epi_start, 
                                                   selection, max_epi, stop_size = treatment_size, 
                                                   d = d, b = b, growth_advantage = growth_advantage)
        
        run_dict <- extinct_runs %>% setNames(1:length(.))
        
        uncontrolled_growth2 <- uncontrolled_growth2 %>% mutate(run = recode(run, !!!run_dict))
        
        uncontrolled_growth <- rbind(uncontrolled_growth %>% filter(!run %in% extinct_runs),
                                     uncontrolled_growth2) %>% arrange(run)
        
        extinct_runs <- uncontrolled_growth %>% group_by(run) %>% filter(max(total) < treatment_size) %>% pull(run) %>% unique
        
        if(keep_extinct){
          extinct_run_dict <- max(extinct_run_dict)+1:length(extinct_runs) %>% setNames(extinct_runs)
          extinct_df <- rbind(extinct_df, 
                              uncontrolled_growth %>% filter(run %in% extinct_runs) %>% 
                                mutate(pRep = pRep, pSeg = pSeg, run = recode(run, !!!extinct_run_dict))
          )
        }
        
    }
    
    initial_conditions <- uncontrolled_growth %>% group_by(run) %>% filter(time == max(time), episomes != -1) %>% 
      mutate(number = frac*total) %>% distinct %>% ungroup %>% 
      select(run, number, episomes) %>% pivot_wider(names_from = episomes, values_from = number) %>% 
      column_to_rownames("run") %>% as.matrix()
    
    start_times <- uncontrolled_growth %>% group_by(run) %>% filter(time == max(time), episomes == 0) %>% 
      distinct() %>% pull(time)
    
    out <- uncontrolled_growth %>% mutate(pRep = pRep, pSeg = pSeg)
    
    if(!is.null(save_baseline_simulations)) write_csv(rbind(out, extinct_df), save_baseline_simulations)
    
    cat("Baseline Simulations Done")
    
  }else{
    
    initial_conditions <- add_to %>% select(-pRep, - pSeg) %>% 
      group_by(run) %>% filter(time == min(time[total == treatment_size]), episomes != -1) %>% 
      mutate(number = frac*total) %>% distinct %>% ungroup %>% 
      select(run, number, episomes) %>% pivot_wider(names_from = episomes, values_from = number) %>% 
      column_to_rownames("run") %>% as.matrix()
    
    start_times <- add_to %>% select(-pRep, - pSeg) %>% group_by(run) %>% 
      filter(time == min(time[total == treatment_size]), episomes == 0) %>% 
      distinct() %>% pull(time)
    
    out <- add_to
  }
  
  
  if(is.null(pRep_reduced)) pRep_reduced <- pRep
  if(is.null(pSeg_reduced)) pSeg_reduced <- pSeg
  
  for(pRep in pRep_reduced){
    for(pSeg in pSeg_reduced){
      cat("\n", "pRep: ", pRep, " pSeg: ", pSeg, "\n")
      intervention <- exponential_growth(pRep, pSeg, 1e8, nRuns, 
                                         selection = selection, max_epi = max_epi, stop_size = stop_size, 
                                         d = d, b = b, growth_advantage = growth_advantage,
                                         initial_conditions = initial_conditions, start_times = start_times, stop_time = stop_time) 
      
      out <- rbind(out, intervention %>% mutate(pRep = pRep, pSeg = pSeg))
      
      if(!is.null(save_treatment_simulations)) write_csv(out, save_treatment_simulations)
      
    }
  }
  
  if(keep_extinct){
    out <- rbind(out, extinct_df)
  }
  
  return(out %>% arrange(run, time))
}

#####
# Plotting functions
#####

# Function to format outputs of above functions
pivot_extinction <- function(extinction_out){
  if(!is.data.frame(extinction_out)) extinction_out <- extinction_out$Totals
  extinction_out %>% 
    pivot_longer(c("zero","one", "two", "three", "four", "five", "six", "seven", "eight", "nine"), 
                 names_to = "episomes_per_cell", values_to = "number_of_cells") %>% 
    mutate(episomes_per_cell = factor(fct_recode(episomes_per_cell, 
                                                 !!!c(`0` = "zero", `1`= "one", `2` = "two", `3`= "three", `4` = "four", 
                                                      `5` = "five", `6`= "six", `7` = "seven", `8` = "eight", `9` = "nine")
    ), levels = 0:9)
    )
}

# Function to plot the number of episomes over time
number_over_time <- function(extinction_long){
  extinction_long %>% ggplot(aes(time, number_of_cells, color = episomes_per_cell)) + 
    geom_point(alpha = 0.25) +  
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    facet_wrap(~episomes_per_cell)
}

fraction_over_time <- function(extinction_long, multiple = F){
  if(!multiple){
    df <- extinction_long %>% 
      group_by(time, trial) %>% 
      mutate(total_cells = sum(number_of_cells)) %>% 
      ungroup %>% mutate(frac = number_of_cells/total_cells) 
  }else{
    df <- extinction_long %>% 
      group_by(time, trial, pRep, pSeg) %>% 
      mutate(total_cells = sum(number_of_cells)) %>% 
      ungroup %>% mutate(frac = number_of_cells/total_cells) 
  }
  df %>% 
    ggplot(aes(time, frac, color = episomes_per_cell)) + 
    geom_point(alpha = 0.25) +  
    guides(color = guide_legend(override.aes = list(alpha = 1))) #+ 
}

