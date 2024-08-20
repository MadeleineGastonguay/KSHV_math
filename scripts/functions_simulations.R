#####
# Simulation functions
#####

## Color palette for simulations

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

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

extinction <- function(pRep, pSeg, nTrials, n_epi, selectAgainstZero = F, n_cells = 1000, n_cells_start = NULL){
  
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
    deathVec = rep(1, 10)
    cells = rep(0, 10)
    cells[n_epi + 1] <- n_cells_start
    total = sum((0:9)*cells)
    indicator <- T
    while(indicator){
      if(time > 700 & pRep == 1 & pSeg == 1) break
      birthVec = rep(2*(1-sum(cells)/n_cells) + 1, 10)
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
      
      # if simulating with selection, run for 700 generations
      # if simulating without selection, run until there are no more episomes
      indicator <- ifelse(selectAgainstZero, time <= 700 & total > 0, total > 0)
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

exponential_growth <- function(pRep, pSeg, nIts, nRuns, n_cells_start = 1, n_epi_start = 3, selection, max_epi){
  #try using matrix first and then convert to data frame
  # rm(data)
  dataCUT = 1
  start = 1:100
  mid = 11:100*10
  end = 2:100*1000
  recordList = c(start, mid, end)
  data <- matrix(NA, nrow=2*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5)
  names(data) = c("run", "time", "episomes", "frac", "total")
  z = 1
  for(run in 1:nRuns){
    print(run)
    deathVec = rep(0, max_epi + 1)
    cells = rep(0, max_epi + 1)
    # Start with n_cells_start cells each with n_epi_start episomes
    cells[n_epi_start+1] <- n_cells_start
    times = rep(0, nIts + 1)
    totalEps = append(c(sum((0:max_epi)*cells)), rep(0, nIts))
    for(j in 1:length(cells)){
      tmp = c(run, times[1], j-1, cells[j]/sum(cells), sum(cells))
      data[z,] <- tmp
      z=z+1
    }
    
    tmp = c(run, times[1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells))
    data[z,] <- tmp
    z= z+1
    #data <- rbind(data, data.frame(run=run, time=times[1], episomes=-1, count = sum((0:9)*cells)))
    for(i in 1:nIts){
      
      birthVec = rep(1, max_epi + 1)
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selection, max_epi)
      cells = result[[2]]
      times[i+1] = times[i] + result[[1]]
      #print(cells)
      totalEps[i+1] = sum((0:max_epi)*cells)
      if(i%%1000 == 0){
        print(i)
      }
      if(sum(cells) %in% recordList){
        #print(i)
        for(j in 1:length(cells)){
          tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
          data[z,] <- tmp
          z = z+1
        }
        tmp = c(run, times[i+1], -1, sum((0:max_epi)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
        data[z,] <- tmp
        z= z+1
        #data <- rbind(data, data.frame(run=run, time=times[i+1], episomes=-1, count = sum((0:9)*cells)))
      }
    }
    
    #plot(times, totalEps)
  }
  # Remove NAs
  data <- data[complete.cases(data), ]
  data <- data.frame(data)
  names(data) <- c("run", "time", "episomes", "frac", "total")
  return(data)
}


#####
# Plotting functions
#####

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

