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

# Advances the population by one Gillespie event: samples which episome class
# the event happens to, whether it is a birth or a death, and how much time
# passes.
#
# The event class is drawn with probability proportional to cells * (birth +
# death) rate, and the waiting time is exponential with rate equal to the total
# propensity across all classes. On a birth, makeChildren() determines how the
# parent's episomes are replicated and partitioned; daughters exceeding max_epi
# are truncated to max_epi.
#
# Arguments:
#   pRep, pSeg      Replication and segregation efficiency, each in [0, 1].
#   cells           Integer vector of length max_epi + 1; element i holds the
#                   number of cells carrying i - 1 episomes.
#   birthVec        Per-class birth rates, same length as cells. Callers vary
#                   this to impose density dependence: extinction() passes a
#                   logistic rate, exponential_growth() a constant one.
#   deathVec        Per-class death rates, same length as cells.
#   selectAgainstZero  If TRUE (default), zero the count of episome-free cells
#                   after every event, so they cannot persist.
#   max_epi         Largest episome count represented. Required.
#   record_selection_deaths  If TRUE, return a third element counting cells
#                   removed by selection this step. Default FALSE.
#
# Returns:
#   list(timeAdvance, cells), plus a third element when
#   record_selection_deaths = TRUE. Returns cells unchanged (not a list) if the
#   population is already empty -- callers should check for this.
simStepFlex <- function(pRep, pSeg, cells, birthVec, deathVec, selectAgainstZero = T, max_epi, record_selection_deaths = F){
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
      if(record_selection_deaths){
        out <- list(timeAdvance, cells, 0)
      }else{
        out <- list(timeAdvance, cells)
      }
      return(out)
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
        if(record_selection_deaths){
          selection_deaths <- cells[1]
        }
        cells[1] <- 0
      }
      
      if(record_selection_deaths){
        out <- list(timeAdvance, cells, selection_deaths)
      }else{
        out <- list(timeAdvance, cells)
      }
      return(out)
    }
  }
}

# Simulates nTrials independent cell populations held at approximately constant
# size (via balanced birth and death rates), tracking the distribution of episomes
# per cell over time. This function can simulate selection against cells without any
# episomes.
#
# Each trial runs until the population loses all episomes, or until stop_time,
# or (when pRep and pSeg are both 1) until 700 generations. The episome
# distribution is recorded every 100 events.
#
# Arguments:
#   pRep               Replication efficiency: probability an episome is
#                      duplicated before division. In [0, 1].
#   pSeg               Segregation efficiency: probability a duplicated episome
#                      is partitioned to the intended daughter. In [0, 1].
#   nTrials            Number of independent populations to simulate.
#   n_epi              Episomes per cell at t = 0; every starting cell gets this
#                      many. Must be <= 9 (the hard cap below).
#   selectAgainstZero  If TRUE, cells carrying zero episomes die immediately.
#                      If FALSE (default), they persist and dilute the
#                      population.
#   n_cells            Carrying capacity the population fluctuates around.
#                      Default 1000. Enters the birth rate as a logistic term.
#   n_cells_start      Starting population size. Defaults to n_cells, i.e. the
#                      population begins at carrying capacity.
#   d                  Per-cell death rate. Default 1.
#   b                  Per-cell birth rate at zero density. Default 3. The
#                      realized birth rate is (b - d)(1 - N/n_cells) + d.
#   stop_time          Optional time limit. NULL (default) runs to episome
#                      extinction.
#   growth_advantage   Optional multiplier on the birth rate, applied to all
#                      cells. NULL (default) means no advantage.
#
# Note: episome number per cell is capped at 9 in this function (the state
# vector is length 10). Use exponential_growth() if you need a larger max_epi.
#
# Returns:
#   A data frame in wide form with columns time, trial, total, and zero..nine
#   giving the number of cells carrying each episome count. Pass to
#   pivot_extinction() to reshape to long form.
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

# Simulates nRuns independent, exponentially growing cell populations, tracking
# the distribution of episomes per cell as each population expands.
#
# Growth is a Gillespie birth-death process with no density dependence, so
# populations grow without bound until stop_size or stop_time is reached. The
# episome distribution is recorded at a pre-set list of population sizes rather
# than at fixed time intervals: densely at small sizes (every cell up to 100,
# then every 10 up to 1000, every 1000 up to 1e5) and progressively more sparsely
# above that. This keeps output size manageable across many orders of magnitude
# of population size.
#
# Arguments:
#   pRep, pSeg        Replication and segregation efficiency, each in [0, 1].
#   nIts              Maximum simulation steps per run. Sizes the internal time
#                     vector; make it large enough that stop_size or stop_time
#                     is reached first.
#   nRuns             Number of independent populations to simulate.
#   n_cells_start     Cells at t = 0. Default 1. Ignored when initial_conditions
#                     is supplied.
#   n_epi_start       Episomes in each starting cell. Default 3. Ignored when
#                     initial_conditions is supplied.
#   selection         If TRUE, cells with zero episomes die immediately.
#                     Required, no default.
#   max_epi           Cap on episomes per cell; the state vector has length
#                     max_epi + 1. Required. Daughters exceeding the cap are
#                     truncated to max_epi, so set it above the largest count
#                     the run will plausibly reach or the distribution's upper
#                     tail will pile up at the cap.
#   stop_size         Population size at which a run terminates. Also determines
#                     how far the recording schedule extends, so it must be set
#                     even when stop_time is the binding constraint.
#   d                 Per-cell death rate. Default 0 (immortal cells; deaths
#                     arise only from selection).
#   b                 Per-cell birth rate. Default 1, which makes one time unit
#                     one mean cell generation.
#   growth_advantage  Optional multiplier on the birth rate. NULL for none.
#   initial_conditions  Optional matrix, one row per run, each row a vector of
#                     length max_epi + 1 giving the number of cells carrying
#                     0, 1, ... max_epi episomes at t = 0. Overrides
#                     n_cells_start and n_epi_start. This is how
#                     simulate_BRK219_experiments.R seeds runs from the fitted
#                     negative binomial (see sample_initial_epi()).
#   start_times       Starting time for each run. A single value is recycled
#                     across runs; a vector of length nRuns sets them
#                     individually. Default 0.
#   stop_time         Optional time limit, in the same units as 1/b. NULL runs
#                     until stop_size.
#   record_selection_deaths  If TRUE, add a selection_deaths column counting
#                     cells removed by selection at each recorded step.
#                     Default FALSE.
#
# Returns:
#   A long data frame with columns run, time, episomes, frac, total (plus
#   selection_deaths when requested). Rows with episomes == -1 hold the mean
#   episomes per cell rather than the frequency of a specific count; filter on
#   this to separate the average from the distribution.
exponential_growth <- function(pRep, pSeg, nIts, nRuns, n_cells_start = 1, n_epi_start = 3, selection, max_epi, 
                               stop_size = NULL, d = 0, b = 1, growth_advantage = NULL, initial_conditions = NULL, 
                               start_times = 0, stop_time = NULL, record_selection_deaths = FALSE){
  
  #try using matrix first and then convert to data frame
  # rm(data)
  dataCUT = 1
  start = 1:100
  mid = 11:100*10
  end = seq(2*1000, 1e5, by = 1000)
  recordList = c(start, mid, end)
  if(stop_size > max(end)){
    recordList = c(recordList, seq(1e5, 1.25e5, by = 2000), seq(1.25e5, min(1e6, stop_size), by = 5000))
    if(stop_size > 1e6){
      step = 100000
      start = 1e6
      stop = 1e7
      for(i in 1:(ceiling(log10(stop_size)) - 6)){
        recordList = c(recordList, seq(start, stop, by = step))  
        step = step*10
        start = start*10
        stop = min(stop*10, stop_size)
      }
    } 
  }
  # data <- matrix(NA, nrow=500*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5)
  data <- matrix(NA, nrow=100*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5+1*record_selection_deaths)
  # data <- matrix(NA, nrow=nIts*nRuns, ncol=5)
  print(dim(data))
  if(record_selection_deaths){
    names(data) = c("run", "time", "episomes", "frac", "total", "selection_deaths")
  }else{
    names(data) = c("run", "time", "episomes", "frac", "total")  
  }
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
      if(record_selection_deaths) tmp = c(tmp, 0)
      try({data[z,] <- tmp})
      z=z+1
    }
    
    tmp = c(run, times[1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells))
    if(record_selection_deaths) tmp = c(tmp, 0)
    try({data[z,] <- tmp})
    z= z+1
    
    #data <- rbind(data, data.frame(run=run, time=times[1], episomes=-1, count = sum((0:9)*cells)))
    for(i in 1:nIts){
      birthVec = rep(b, max_epi + 1)
      if(!is.null(growth_advantage)) birthVec = birthVec*c(1, rep(1 + growth_advantage, max_epi))
      
      cells2 = cells
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selection, max_epi, record_selection_deaths)
      
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
          if(record_selection_deaths) tmp = c(tmp, result[[3]])
          try({data[z,] <- tmp})
          z = z+1
        }
        tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
        if(record_selection_deaths) tmp = c(tmp, result[[3]])
        try({data[z,] <- tmp})
        z= z+1
        #data <- rbind(data, data.frame(run=run, time=times[i+1], episomes=-1, count = sum((0:9)*cells)))
        
      }
      
      # print(sum(cells))
      
      if(sum(cells) == 0){ # add final record and break
        for(j in 1:length(cells)){
          tmp = c(run, times[i+1], j-1, 0, sum(cells)) #added sum cells
          if(record_selection_deaths) tmp = c(tmp, result[[3]])
          try({data[z,] <- tmp})
          z = z+1
        }
        tmp = c(run, times[i+1], -1, 0, sum(cells)) #totalEps[i+1])
        if(record_selection_deaths) tmp = c(tmp, result[[3]])
        try({data[z,] <- tmp})
        break
      } 
      
      if(!is.null(stop_time)){
        if(times[i+1] >= stop_time + start_times[run]){
          for(j in 1:length(cells)){
            tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
            if(record_selection_deaths) tmp = c(tmp, result[[3]])
            try({data[z,] <- tmp})
            z = z+1
          }
          tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
          if(record_selection_deaths) tmp = c(tmp, result[[3]])
          try({data[z,] <- tmp})
          break
        } 
      }
      
      if(!is.null(stop_size)){
        if(sum(cells) >= stop_size){
          for(j in 1:length(cells)){
            tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
            if(record_selection_deaths) tmp = c(tmp, result[[3]])
            try({data[z,] <- tmp})
            z = z+1
          }
          tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
          if(record_selection_deaths) tmp = c(tmp, result[[3]])
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
  if(record_selection_deaths) {
    names(data) <- c("run", "time", "episomes", "frac", "total", "selection_deaths")
  }else{
    names(data) <- c("run", "time", "episomes", "frac", "total")
  }
  return(data)
}


# Simulates nRuns KSHV-dependent tumors growing from n_cells_start cells, then
# applies a "treatment" that reduces replication and/or segregation efficiency
# once the tumor reaches treatment_size.
#
# Runs that go extinct before reaching treatment_size are re-simulated until
# nRuns tumors survive to treatment, so the treated cohort is always full size.
# Set keep_extinct = TRUE to retain the discarded trajectories instead.
#
# Arguments:
#   pRep, pSeg           Baseline replication and segregation efficiency, used
#                        for growth up to treatment_size. In [0, 1].
#   pRep_reduced         Post-treatment replication efficiency. May be a vector
#                        to simulate several treatment strengths in one call;
#                        pass c(reduced, pRep) to include an untreated control.
#                        NULL leaves replication unchanged.
#   pSeg_reduced         Post-treatment segregation efficiency, same convention.
#                        NULL leaves segregation unchanged.
#   nRuns                Number of tumors that must survive to treatment.
#   n_cells_start        Cells at t = 0. Default 1.
#   n_epi_start          Episomes in each starting cell. Default 3.
#   selection            If TRUE, cells with zero episomes die immediately (the
#                        KSHV-dependent tumor assumption). Required, no default.
#   max_epi              Cap on episomes per cell; the state vector has length
#                        max_epi + 1. Required. Book-keeping only, so set it
#                        above the largest count the simulation will reach.
#   stop_size            Population size at which a run terminates. Default
#                        1.5e5.
#   treatment_size       Population size at which the reduced efficiencies take
#                        effect. Default 1e5.
#   d                    Per-cell death rate. Default 0 (immortal cells; deaths
#                        come only from selection).
#   b                    Per-cell birth rate. Default 1, so time is in units of
#                        cell generations.
#   growth_advantage     Optional multiplier on the birth rate. NULL for none.
#   add_to               Optional output from a previous call. When supplied,
#                        the baseline growth phase is skipped and treatment is
#                        applied to those existing trajectories, which saves
#                        re-simulating the pre-treatment phase across scenarios.
#   stop_time            Optional time limit, in the same units as 1/b.
#   keep_extinct         If TRUE, retain trajectories that died before reaching
#                        treatment_size, renumbered after the surviving runs.
#                        Default FALSE discards them.
#   save_baseline_simulations   Optional path to save pre-treatment
#                        trajectories to, for reuse via add_to.
#   save_treatment_simulations  Optional path to save post-treatment
#                        trajectories to.
#   record_selection_deaths     If TRUE, record how many cells were removed by
#                        selection at each step. Default FALSE.
#
# Returns:
#   A long data frame with columns run, time, episomes, frac, total, plus pRep
#   and pSeg identifying the treatment scenario.
PEL_simulations <- function(pRep, pSeg, pRep_reduced = NULL, pSeg_reduced = NULL, nRuns, n_cells_start = 1, n_epi_start = 3, selection, max_epi, 
                            stop_size = 1.5e5, treatment_size = 1e5, d = 0, b = 1, growth_advantage = NULL, add_to = NULL, stop_time = NULL,
                            keep_extinct = FALSE, save_baseline_simulations = NULL, save_treatment_simulations = NULL, record_selection_deaths = FALSE){
  
  if(is.null(add_to)){
  
      uncontrolled_growth <- exponential_growth(pRep, pSeg, 1e8, nRuns, 
                                                n_cells_start = n_cells_start, n_epi_start = n_epi_start, 
                                                selection, max_epi, stop_size = treatment_size, 
                                                d = d, b = b, growth_advantage = growth_advantage,
                                                record_selection_deaths = record_selection_deaths)
      
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
                                                   d = d, b = b, growth_advantage = growth_advantage,
                                                   record_selection_deaths = record_selection_deaths)
        
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
                                         initial_conditions = initial_conditions, start_times = start_times, stop_time = stop_time,
                                         record_selection_deaths = record_selection_deaths) 
      
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

# Plots the fraction of the population carrying each episome count over time.
#
# Companion to number_over_time(), which plots absolute cell counts instead.
# Fractions are computed within each timepoint, so they sum to 1 at every time.
#
# Arguments:
#   extinction_long   Long-form simulation output from pivot_extinction(), with
#                     columns time, trial, episomes_per_cell, number_of_cells.
#   multiple          Set TRUE when the data frame pools several parameter
#                     combinations, so that fractions are computed within each
#                     (time, trial, pRep, pSeg) group rather than within
#                     (time, trial) alone. Default FALSE.
#
# Returns:
#   A ggplot object.
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

