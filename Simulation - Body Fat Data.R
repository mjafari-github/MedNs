##### REAL DATA: BODY FAT DATA #####
##### February 28, 2024

#################### LOADING PACKAGES / FUNCTIONS ####################
library(MASS)
library(ald)
library(quantreg)
source(here::here("Functions", "M-estimation Functions.R"))
source(here::here("Functions", "RSS Functions.R"))
source(here::here("Functions", "Other Functions.R"))

# Settings
B = 500 # number of replications for simulation
tau = 0.5 # quantile
n = 40 # sample size

#################### START OF POPULATION ####################
{
  # load data
  BodyFat <- read.csv(here::here("Data", "Body Fat Data", "BodyFat - Extended.csv"))
  
  BodyFat <- BodyFat[,-2] # remove 'original' variable
    
  BodyFat[sapply(BodyFat, is.character)] <- 
    lapply(BodyFat[sapply(BodyFat, is.character)], as.factor)
  
  BodyFat[sapply(BodyFat, is.factor)] <- 
    lapply(BodyFat[sapply(BodyFat, is.factor)], as.numeric)
  
  BodyFat[,2] = BodyFat[,2]-1 # 0/1 coding for factors
  
  BodyFat <- cbind(BodyFat, "BMI" = BodyFat$Weight/BodyFat$Height^2)
  
  N = nrow(BodyFat) # population size
  
  y = BodyFat$BodyFat # response
  features = as.matrix(BodyFat[,6:15]) # anthropometric features not used in ranker
  
  R1 = 1.2*BodyFat$BMI + 0.23*BodyFat$Age - 10.8*BodyFat$Sex - 5.4
  R2 = (ifelse(BodyFat$Sex == 1, 1.407*BodyFat$BMI - 21.389, 1.9337*BodyFat$BMI - 26.422)/BodyFat$Weight)*100
  R3 = 1.506*BodyFat$BMI + 0.133*BodyFat$Age - 11.481*BodyFat$Sex - 11.52
  R4 = 1.61*BodyFat$BMI + 0.13*BodyFat$Age - 12.1*BodyFat$Sex - 13.9
  
  rankers = cbind(R1)
  
  X = cbind(1, features)
  population = cbind(y, features)
  d = ncol(X) # number of regression coefficients
  R = ncol(rankers)
}
#################### END OF POPULATION ####################


#################### START OF TAU VALUES ####################
{
  # taking pilot sample
  n.pilot = 40
  pilot.ind <- sample(N, size = n.pilot)
  
  ##### TAU VALUES FOR SRS APPROACH #####
  ### IMPERFECT RANKING
  
  # store results for ranking probabilities in list
  p.matrix=list()
  
  # storing values for tau for k and ranker
  tau.k.imperfect = matrix(NA, nrow = R, ncol = 3)
  #rownames(tau.k.imperfect) = c("Weight")
  colnames(tau.k.imperfect) = c("k=3", "k=5", "k=7")
  
  for(i in 1:R){
    
    ranker = rankers[, i] # select ranker
    
    for(k in c(3, 5, 7)){
      
      # replicate ordering of SRS 
      reps = 1000
      probEst <- replicate(reps, ranking.prob(y[pilot.ind], ranker[pilot.ind], k=k), simplify=FALSE)
      
      # average probabilities
      p.matrix[[i]] <- Reduce("+", probEst)/reps
      
      # vectorizing calculations for H_r
      r = 1:k
      H.terms = choose(k, r)*tau^r*(1-tau)^(k-r)
      H_r = rev(cumsum(rev(H.terms)))
      
      alpha_r = p.matrix[[i]][(k+1)/2,] # estimates of ranking probabilities for median
      
      tau.k.imperfect[i, (k-1)/2] = sum(H_r*alpha_r)
    }
  }
  
}
#################### END OF TAU VALUES ####################


#################### START OF SIMULATION ####################
{
  ### Store results for coefficients
  
  # MR and QR to compare to
  coef.mean <- matrix(NA, nrow = B, ncol = d)
  coef.QR <- matrix(NA, nrow = B, ncol = d)
  
  # Original MRSS Loss results
  coef.K3 <- matrix(NA, nrow = B, ncol = R*d)
  coef.K5 <- matrix(NA, nrow = B, ncol = R*d)
  coef.K7 <- matrix(NA, nrow = B, ncol = R*d)
  coef.K3SRS <- matrix(NA, nrow = B, ncol = d)
  coef.K5SRS <- matrix(NA, nrow = B, ncol = d)
  coef.K7SRS <- matrix(NA, nrow = B, ncol = d)
  
  # Approximate MRSS Loss results
  coef.K3approx <- matrix(NA, nrow = B, ncol = R*d)
  coef.K5approx <- matrix(NA, nrow = B, ncol = R*d)
  coef.K7approx <- matrix(NA, nrow = B, ncol = R*d)
  coef.K3SRSapprox <- matrix(NA, nrow = B, ncol = d)
  coef.K5SRSapprox <- matrix(NA, nrow = B, ncol = d)
  coef.K7SRSapprox <- matrix(NA, nrow = B, ncol = d)
  
  # Changing tau for SRS from induced population
  coef.K3SRSinduced <- matrix(NA, nrow = B, ncol = R*d)
  coef.K5SRSinduced <- matrix(NA, nrow = B, ncol = R*d)
  coef.K7SRSinduced <- matrix(NA, nrow = B, ncol = R*d)
}
##### BEGIN REPLICATIONS
{
  
  start = Sys.time()
  
  for(b in 1:B){
    
    set.seed(b)
    
    #### Setting aside SRS data
    
    ### SRS DATA (Training)
    ind <- sample(1:N, size = n)
    
    # training data (SRS)
    population.trainSRS <- as.data.frame(population[ind,])
    y.trainSRS <- y[ind]
    features.trainSRS <- features[ind,]
    X.trainSRS <- X[ind,]
    
    
    ##### FIT MEAN REGRESSION TO SRS DATA
    MeanReg <- lm(y.trainSRS ~ features.trainSRS)
    coef.mean[b,] <- MeanReg$coefficients
    
    ##### FIT MEDIAN REGRESSION TO SRS DATA
    QRcoef <- mrss.M.estimation(response = y.trainSRS, 
                                features = features.trainSRS, 
                                tau = tau, m = 1, max.diff = 1e-3) # NOTE: k=1 is equivalent to standard check function
    coef.QR[b,] <- QRcoef
    
    ##### MRSS REGRESSION USING SRS DATA AS INPUT
    
    ## k=3
    if(tau <= 2/3 && tau >= 1/3){
      # Original MRSS Loss
      model.coef.3SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=2)
      coef.K3SRS[b, ] <- model.coef.3SRS
      # Approximate Loss
      model.coef.3SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=2)
      coef.K3SRSapprox[b, ] <- model.coef.3SRSapprox
    }
    ## k=5
    if(tau <= 3/5 && tau >= 2/5){
      # Original MRSS Loss
      model.coef.5SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=3)
      coef.K5SRS[b, ] <- model.coef.5SRS
      # Approximate Loss
      model.coef.5SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=3)
      coef.K5SRSapprox[b, ] <- model.coef.5SRSapprox
    }
    ## k=7
    if(tau <= 4/7 && tau >= 3/7){
      # Original MRSS Loss
      model.coef.7SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=4)
      coef.K7SRS[b, ] <- model.coef.7SRS
      # Approximate Loss
      model.coef.7SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=4)
      coef.K7SRSapprox[b, ] <- model.coef.7SRSapprox
    }
    
    
    ##### FIT MEDIAN REGRESSION TO MRSS DATA (with each ranker)
    coef3 <- coef5 <- coef7 <- c()
    coef3approx <- coef5approx <- coef7approx <- c()
    coef3SRSinduced <- coef5SRSinduced <- coef7SRSinduced <- c()
    for(r in 1:R){
      ranker = rankers[,r]
      
      # take training data for each MRSS model
      ind.trainK3 <- unlist(my.RSS(population, ranker, k = 3, 
                                   ni.size = c(0, n, 0))) # m=2
      ind.trainK5 <- unlist(my.RSS(population, ranker, k = 5, 
                                   ni.size = c(0, 0, n, 0, 0))) # m=3
      ind.trainK7 <- unlist(my.RSS(population, ranker, k = 7, 
                                   ni.size = c(0, 0, 0, n, 0, 0, 0))) # m=4
      
      # training data (k=3)
      population.trainK3 <- as.data.frame(population[ind.trainK3,])
      y.trainK3 <- y[ind.trainK3]
      features.trainK3 <- features[ind.trainK3,]
      X.trainK3 <- X[ind.trainK3,]
      
      # training data (k=5)
      population.trainK5 <- as.data.frame(population[ind.trainK5,])
      y.trainK5 <- y[ind.trainK5]
      features.trainK5 <- features[ind.trainK5,]
      X.trainK5 <- X[ind.trainK5,]
      
      # training data (k=7)
      population.trainK7 <- as.data.frame(population[ind.trainK7,])
      y.trainK7 <- y[ind.trainK7]
      features.trainK7 <- features[ind.trainK7,]
      X.trainK7 <- X[ind.trainK7,]
      
      ##### MRSS REGRESSION
      
      ## k=3
      if(tau <= 2/3 && tau >= 1/3){
        # Original
        model.coef.3 <- mrss.M.estimation(response = y.trainK3, 
                                          features = features.trainK3, 
                                          tau = tau, m=2)
        coef3 <- c(coef3, model.coef.3)
        # Approximation
        model.coef.3approx <- mrssAPPROX.M.estimation(response = y.trainK3, 
                                                      features = features.trainK3, 
                                                      tau = tau, m=2)
        coef3approx <- c(coef3approx, model.coef.3approx)
      }
      ## k=5
      if(tau <= 3/5 && tau >= 2/5){
        # Original
        model.coef.5 <- mrss.M.estimation(response = y.trainK5, 
                                          features = features.trainK5, 
                                          tau = tau, m=3)
        coef5 <- c(coef5, model.coef.5)
        # Approximation
        model.coef.5approx <- mrssAPPROX.M.estimation(response = y.trainK5, 
                                                      features = features.trainK5, 
                                                      tau = tau, m=3)
        coef5approx <- c(coef5approx, model.coef.5approx)
      }
      ## k=7
      if(tau <= 4/7 && tau >= 3/7){
        # Original
        model.coef.7 <- mrss.M.estimation(response = y.trainK7, 
                                          features = features.trainK7, 
                                          tau = tau, m=4)
        coef7 <- c(coef7, model.coef.7)
        # Approximation
        model.coef.7approx <- mrssAPPROX.M.estimation(response = y.trainK7, 
                                                      features = features.trainK7, 
                                                      tau = tau, m=4)
        coef7approx <- c(coef7approx, model.coef.7approx)
      }
      
      
      ##### SRS from induced population Approach
      
      # Using quantreg package for quantile regression
      model.3SRSinduced <- rq(y.trainK3 ~ features.trainK3, tau = tau.k.imperfect[r, 1])
      model.5SRSinduced <- rq(y.trainK5 ~ features.trainK5, tau = tau.k.imperfect[r, 2])
      model.7SRSinduced <- rq(y.trainK7 ~ features.trainK7, tau = tau.k.imperfect[r, 3])

      # record coefficients
      coef3SRSinduced <- c(coef3SRSinduced, model.3SRSinduced$coefficients)
      coef5SRSinduced <- c(coef5SRSinduced, model.5SRSinduced$coefficients)
      coef7SRSinduced <- c(coef7SRSinduced, model.7SRSinduced$coefficients)
    }
    
    ### record coefficients
    ## k=3
    if(tau <= 2/3 && tau >= 1/3){
      coef.K3[b, ] <- coef3
      coef.K3approx[b, ] <- coef3approx
    }
    ## k=5
    if(tau <= 3/5 && tau >= 2/5){
      coef.K5[b, ] <- coef5
      coef.K5approx[b, ] <- coef5approx
    }
    ## k=7
    if(tau <= 4/7 && tau >= 3/7){
      coef.K7[b, ] <- coef7
      coef.K7approx[b, ] <- coef7approx
    }
    
    coef.K3SRSinduced[b, ] <- coef3SRSinduced
    coef.K5SRSinduced[b, ] <- coef5SRSinduced
    coef.K7SRSinduced[b, ] <- coef7SRSinduced
    
    if(b %% 100 == 0){ print(paste("iteration", sep = " = ", b)) } # counts iterations
    
  }
  
  end = Sys.time()
  
  end-start # time for simulation
  
}

#################### END OF SIMULATION ####################


#################### START OF STORING COEFFICIENTS ####################

coef.all <- cbind(coef.mean, coef.QR,
                  coef.K3, coef.K5, coef.K7,
                  coef.K3approx, coef.K5approx, coef.K7approx,
                  coef.K3SRSinduced, coef.K5SRSinduced, coef.K7SRSinduced,
                  coef.K3SRS, coef.K5SRS, coef.K7SRS,
                  coef.K3SRSapprox, coef.K5SRSapprox, coef.K7SRSapprox)

# Note: no bias adjustment needed when tau = 0.5

#write.csv(coef.all, here::here("Data", "Body Fat Data", paste0("CoefResults-BodyFat", n, ".csv")))

#################### END OF STORING COEFFICIENTS ####################

