##### RSS FUNCTION #####

my.RSS <- function(population, ranker, k, ni.size){
  # population is a data frame of observations to draw samples from
  # ranker is the variable which will be used to rank the SRS observations
  # k is the set size
  # ni.size is a vector of length k representing the number of observations to sample from each order statistic
  # note that the total sample size is not specified; it is given by the sum of ni's
  
  N = nrow(population) # size of the population
  remaining.ind = 1:N # indices to sample from
  
  ## ERRORS
  if(length(ni.size) != k){ stop("ni.size must be of length k") } # returns error is ni.size is not defined properly
  if(sum(ni.size) > N){ stop("sample size exceeds population size") } # returns error if n is too large
  
  RSS.ind = list() # hold onto sample indices in list
  for(i in 1:k){
    
    ni = ni.size[i] # select sample size for i-th order statistic
    
    if(ni != 0){ # If ni is 0, then skip sampling for that group
      
      sample.ind = c() # holds sample indices for i-th order statistic
      for(j in 1:ni){
        ind <- sample(remaining.ind, k) # take sample of set size k
        sample.ind[j] = ind[order(ranker[ind])][i] # orders SRS indices by ranker value and takes i-th entry
        remaining.ind = remaining.ind[-match(sample.ind[j], remaining.ind)] # remove the chosen observation from population
      }
      RSS.ind[[i]] = sample.ind # store indices in list
    }
  }
  # each list entry in output is the sampled indices for each order statistic
  return(RSS.ind)
}


##### MOVING RSS FUNCTION
##### Date: October 25, 2023

moving.RSS <- function(population, ranker, nrank){
  
  # store each set of indices in a list item
  movingRSS <- list()
  for(r in 1:length(nrank)){
    
    # creating MRSS vector for ni.size
    ni.size <- rep(0, 2*r-1)
    ni.size[r] <- nrank[r]
    
    # Take RSS for each m value
    movingRSS[[r]] <- unlist(my.RSS(population = population, ranker = ranker, 
                                    k = 2*r-1, ni.size = ni.size))
  }
  return(movingRSS)
}





##### FUNCTION FOR ESTIMATING RANKING ERROR PROBABILITIES #####

ranking.prob <- function(response, ranker, k){
  # response: variable of interest from population - used for perfect ranking
  # ranker: variable which will be used to rank the SRS observations - used for imperfect ranking
  # k: set size

  # NOTE: This function provides only one ordered sample and
  # and must be replicated to estimate these probabilities
  
  N = length(response) # size of the population
  
  ind <- sample(1:N, k, replace = TRUE) # take sample of set size k
  
  imperfect.set = ind[order(ranker[ind])] # orders SRS indices by ranker value
  perfect.set = ind[order(response[ind])] # orders SRS indices by true value
  
  p.counts = matrix (NA, k, k) # store the counts of classifications
  
  # For each entry of the p.counts matrix, we check if the 
  # observations in the imperfect set match the perfect set 
  # e.g., check if the first observation in the imperfect set matches the first in the perfect set
  # then check if the first observation in the imperfect set matches the second in the perfect set
  # and so on until all observation pairs are checked
  for(i in 1:k){  
    for(j in 1:k){
      p.counts[i, j] <- sum(imperfect.set[i] == perfect.set[j]) 
    }
  }
  return(p.counts)
}