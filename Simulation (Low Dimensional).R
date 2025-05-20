##### SIMULATION: FINAL (with heteroskedasticity & IMSE) #####
##### MRSS loss continuous covariate population
##### Original and approximate loss included
##### Low dimensional Setting
##### FINAL CONTINUOUS SIMULATION SETTING
##### May 9, 2024

#################### LOADING PACKAGES / FUNCTIONS ####################

# load packages
library(MASS)
library(ald)
library(quantreg)
source(here::here("Functions", "M-estimation Functions.R"))
source(here::here("Functions", "RSS Functions.R"))
source(here::here("Functions", "Other Functions.R"))


# Settings
if(F){
  B = 100 # number of replications for simulation
  tau = 0.5 # quantile
  case = 5 # choose 1 through 6 for different error setting in paper
}


#################### START OF GENERATING POPULATION ####################

### Population from "Robust Linear Regression: A Review and Comparison"

##### CONSTANT VARIANCE

##### CH3 SIMULATION
### three predictors
if(T){
  set.seed(70)
  
  beta = c(0, 1, 1, 1) # parameter values
  gamma = c(1, 0, 0, 0)
  N = 10000 # population size
  
  # generate covariates
  X1 <- rnorm(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  X3 <- rnorm(N, 0, 1)
  X <- cbind(1, X1, X2, X3)
  
  # generate y with error and rankers
  rho.val = c(0.9, 0.8, 0.7, 0.6) # selected correlation of rankers and response
  
  if(case == 1){ # Normal
    y <- as.vector( X %*% beta + rnorm(N) )
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 2){ # t(3)
    # Add t(3) error to response
    y <- as.vector( X %*% beta + rt(N, 3))
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 3){ # t(1)
    # Add t(1) error to response
    y <- as.vector( X %*% beta + rt(N, 1) )
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 4){ # Mixture
    # Add contaminated mixture error to response
    y <- as.vector( X %*% beta + 0.95*rnorm(N) + 0.05*rnorm(N, 0, 10^2) )
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 5){ # 10% y outliers
    # Add Normal error to response
    y <- as.vector( X %*% beta + rnorm(N) )
    # Add outliers
    y[1:(0.1*N)] <- 30
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 6){ # 1% x and y outliers
    # Add Normal error to response
    y <- as.vector( X %*% beta + rnorm(N) )
    # Add outliers
    X[1:(0.01*N),2:4] <- 10
    y[1:(0.01*N)] <- 50
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  }
  # Objects used in simulation
  features = as.matrix(X[,-1])
  population = cbind(y, features)
  rankers = cbind(R1, R2, R3, R4)
  d = length(beta)
  n = 200 # sample size
}

##### CH2 SIMULATION:
### single predictor with heteroskedasticity
if(F){
  set.seed(40)
  
  beta = c(2, 3) # parameter values
  gamma = c(1, 1) # parameters for variance
  N = 1000 # population size
  N.out = 0.1*N # number of outliers
  
  # generate covariates
  X <- runif(N, min = 0, max = 5)
  X <- cbind(1, X)
  
  # generate y with error and rankers
  rho.val = c(0.9, 0.8, 0.7, 0.6) # selected correlation of rankers and response
  
  if(case == 1){ # Normal
    y <- as.vector( X %*% beta + X %*% gamma * rnorm(N) )
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 5){ # 10% outliers in y direction
    # Add Normal error to response
    y <- as.vector( X %*% beta +  X %*% gamma * rnorm(N) )
    # Add outliers
    out.data <- mvrnorm(N.out, mu=c(4, 50), Sigma=diag(0.25*c(1, 20)))
    X[1:N.out,2] <- out.data[,1]
    y[1:N.out] <- out.data[,2]
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  } else if(case == 6){ # 10% outliers in both x and y direction
    # Add Normal error to response
    y <- as.vector( X %*% beta +  X %*% gamma * rnorm(N) )
    # Add outliers
    out.data <- mvrnorm(N.out, mu=c(6, 50), Sigma=diag(0.25*c(1, 20)))
    X[1:N.out,2] <- out.data[,1]
    y[1:N.out] <- out.data[,2]
    # Simulate rankers
    R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
    R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
    R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
    R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  }
  # Objects used in simulation
  features = as.matrix(X[,-1])
  population = cbind(y, features)
  rankers = cbind(R1, R2, R3, R4)
  d = length(beta)
  n = 50 # sample size
}

#################### END OF GENERATING POPULATION #################### 


#################### START OF TAU VALUES ####################
{
  ##### TAU VALUES FOR SRS APPROACH #####
  ### IMPERFECT RANKING
  
  # store results for ranking probabilities in list
  p.matrix=list()
  
  # storing values for tau for k and ranker
  tau.k.imperfect = matrix(NA, nrow = ncol(rankers), ncol = 3)
  rownames(tau.k.imperfect) = c("R1", "R2", "R3", "R4")
  colnames(tau.k.imperfect) = c("k=3", "k=5", "k=7")
  
  for(i in 1:ncol(rankers)){
    
    ranker = rankers[, i] # select ranker
    
    for(k in c(3, 5, 7)){
      
      # replicate ordering of SRS 
      reps = 1000
      probEst <- replicate(reps, ranking.prob(y, ranker, k=k), simplify=FALSE)
      
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
  coef.K3 <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5 <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7 <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K3SRS <- matrix(NA, nrow = B, ncol = d)
  coef.K5SRS <- matrix(NA, nrow = B, ncol = d)
  coef.K7SRS <- matrix(NA, nrow = B, ncol = d)
  
  # Approximate MRSS Loss results
  coef.K3approx <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5approx <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7approx <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K3SRSapprox <- matrix(NA, nrow = B, ncol = d)
  coef.K5SRSapprox <- matrix(NA, nrow = B, ncol = d)
  coef.K7SRSapprox <- matrix(NA, nrow = B, ncol = d)
  
  # Changing tau for SRS from induced population
  coef.K3SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  
  # Using MedNS data in quantile regression (check loss)
  coef.K3QR <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5QR <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7QR <- matrix(NA, nrow = B, ncol = 4*d)
  
  ### Store results for IMSEs
  
  # MR and QR to compare to
  IMSE.QR <- c()
  
  # Original MRSS Loss results
  IMSE.K3 <- IMSE.K5 <- IMSE.K7 <- matrix(NA, nrow = B, ncol = 4)
  IMSE.K3SRS <- IMSE.K5SRS <- IMSE.K7SRS <- c()
  
  # Approximate MRSS Loss results
  IMSE.K3approx <- IMSE.K5approx <- IMSE.K7approx <- matrix(NA, nrow = B, ncol = 4)
  IMSE.K3SRSapprox <- IMSE.K5SRSapprox <- IMSE.K7SRSapprox <- c()
  
  # Changing tau for SRS from induced population
  IMSE.K3SRSinduced <- IMSE.K5SRSinduced <- IMSE.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4)
  
  # Using MedNS data in quantile regression (check loss)
  IMSE.K3QR <- IMSE.K5QR <- IMSE.K7QR <- matrix(NA, nrow = B, ncol = 4)
}
##### BEGIN REPLICATIONS
{
  
  start = Sys.time()
  
  for(b in 1:B){
    
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
    IMSE.QR[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% QRcoef)^2)
    
    #sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% QRcoef)^2)
    #sum(diag(X.trainSRS %*% (beta_calc(tau) - QRcoef) %*% t(beta_calc(tau) - QRcoef) %*% t(X.trainSRS)))
    
    ##### MRSS REGRESSION USING SRS DATA AS INPUT
    
    ## k=3
    if(tau <= 2/3 && tau >= 1/3){
      # Original MRSS Loss
      model.coef.3SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=2)
      coef.K3SRS[b, ] <- model.coef.3SRS
      model.coef.3SRS[1] <- model.coef.3SRS[1] + (1-2*tau)
      IMSE.K3SRS[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.3SRS)^2)
      # Approximate Loss
      model.coef.3SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=2)
      coef.K3SRSapprox[b, ] <- model.coef.3SRSapprox
      model.coef.3SRSapprox[1] <- model.coef.3SRSapprox[1] + (1-2*tau)
      IMSE.K3SRSapprox[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.3SRSapprox)^2)
    }
    ## k=5
    if(tau <= 3/5 && tau >= 2/5){
      # Original MRSS Loss
      model.coef.5SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=3)
      coef.K5SRS[b, ] <- model.coef.5SRS
      model.coef.5SRS[1] <- model.coef.5SRS[1] + (1-2*tau)
      IMSE.K5SRS[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.5SRS)^2)
      # Approximate Loss
      model.coef.5SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=3)
      coef.K5SRSapprox[b, ] <- model.coef.5SRSapprox
      model.coef.5SRSapprox[1] <- model.coef.5SRSapprox[1] + (1-2*tau)
      IMSE.K5SRSapprox[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.5SRSapprox)^2)
    }
    ## k=7
    if(tau <= 4/7 && tau >= 3/7){
      # Original MRSS Loss
      model.coef.7SRS <- mrss.M.estimation(response = y.trainSRS, 
                                           features = features.trainSRS, 
                                           tau = tau, m=4)
      coef.K7SRS[b, ] <- model.coef.7SRS
      model.coef.7SRS[1] <- model.coef.7SRS[1] + (1-2*tau)
      IMSE.K7SRS[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.7SRS)^2)
      # Approximate Loss
      model.coef.7SRSapprox <- mrssAPPROX.M.estimation(response = y.trainSRS, 
                                                       features = features.trainSRS, 
                                                       tau = tau, m=4)
      coef.K7SRSapprox[b, ] <- model.coef.7SRSapprox
      model.coef.7SRSapprox[1] <- model.coef.7SRSapprox[1] + (1-2*tau)
      IMSE.K7SRSapprox[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% model.coef.7SRSapprox)^2)
    }
    
    
    ##### FIT MEDIAN REGRESSION TO MRSS DATA (with each ranker)
    coef3 <- coef5 <- coef7 <- c()
    coef3approx <- coef5approx <- coef7approx <- c()
    coef3SRSinduced <- coef5SRSinduced <- coef7SRSinduced <- c()
    coef3QR <- coef5QR <- coef7QR <- c()
    for(r in 1:ncol(rankers)){
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
        model.coef.3[1] <- model.coef.3[1] + (1-2*tau)
        IMSE.K3[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% model.coef.3)^2)
        # Approximation
        model.coef.3approx <- mrssAPPROX.M.estimation(response = y.trainK3, 
                                                      features = features.trainK3, 
                                                      tau = tau, m=2)
        coef3approx <- c(coef3approx, model.coef.3approx)
        model.coef.3approx[1] <- model.coef.3approx[1] + (1-2*tau)
        IMSE.K3approx[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% model.coef.3approx)^2)
      }
      ## k=5
      if(tau <= 3/5 && tau >= 2/5){
        # Original
        model.coef.5 <- mrss.M.estimation(response = y.trainK5, 
                                          features = features.trainK5, 
                                          tau = tau, m=3)
        coef5 <- c(coef5, model.coef.5)
        model.coef.5[1] <- model.coef.5[1] + (1-2*tau)
        IMSE.K5[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% model.coef.5)^2)
        # Approximation
        model.coef.5approx <- mrssAPPROX.M.estimation(response = y.trainK5, 
                                                      features = features.trainK5, 
                                                      tau = tau, m=3)
        coef5approx <- c(coef5approx, model.coef.5approx)
        model.coef.5approx[1] <- model.coef.5approx[1] + (1-2*tau)
        IMSE.K5approx[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% model.coef.5approx)^2)
      }
      ## k=7
      if(tau <= 4/7 && tau >= 3/7){
        # Original
        model.coef.7 <- mrss.M.estimation(response = y.trainK7, 
                                          features = features.trainK7, 
                                          tau = tau, m=4)
        coef7 <- c(coef7, model.coef.7)
        model.coef.7[1] <- model.coef.7[1] + (1-2*tau)
        IMSE.K7[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% model.coef.7)^2)
        # Approximation
        model.coef.7approx <- mrssAPPROX.M.estimation(response = y.trainK7, 
                                                      features = features.trainK7, 
                                                      tau = tau, m=4)
        coef7approx <- c(coef7approx, model.coef.7approx)
        model.coef.7approx[1] <- model.coef.7approx[1] + (1-2*tau)
        IMSE.K7approx[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% model.coef.7approx)^2)
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
      
      IMSE.K3SRSinduced[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% model.3SRSinduced$coefficients)^2)
      IMSE.K5SRSinduced[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% model.5SRSinduced$coefficients)^2)
      IMSE.K7SRSinduced[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% model.7SRSinduced$coefficients)^2)
      
      
      
      ##### Using MedNS data in quantile regression (check loss)
      
      # Using quantreg package for quantile regression
      model.3QR <- rq(y.trainK3 ~ features.trainK3, tau = tau)
      model.5QR <- rq(y.trainK5 ~ features.trainK5, tau = tau)
      model.7QR <- rq(y.trainK7 ~ features.trainK7, tau = tau)
      
      # record coefficients
      coef3QR <- c(coef3QR, model.3QR$coefficients)
      coef5QR <- c(coef5QR, model.5QR$coefficients)
      coef7QR <- c(coef7QR, model.7QR$coefficients)
      
      IMSE.K3QR[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% model.3QR$coefficients)^2)
      IMSE.K5QR[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% model.5QR$coefficients)^2)
      IMSE.K7QR[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% model.7QR$coefficients)^2)
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
    
    coef.K3QR[b, ] <- coef3QR
    coef.K5QR[b, ] <- coef5QR
    coef.K7QR[b, ] <- coef7QR
    
    if(b %% 50 == 0){ print(paste("iteration", sep = " = ", b)) } # counts every 5 iterations
    
  }
  
  end = Sys.time()
  
  end-start # time for simulation
  
}

#################### END OF SIMULATION ####################

#################### START OF STORING COEFFICIENTS ####################

coef.all <- cbind(coef.mean, coef.QR,
                  coef.K3, coef.K5, coef.K7,
                  coef.K3SRS, coef.K5SRS, coef.K7SRS,
                  coef.K3approx, coef.K5approx, coef.K7approx,
                  coef.K3SRSapprox, coef.K5SRSapprox, coef.K7SRSapprox,
                  coef.K3SRSinduced, coef.K5SRSinduced, coef.K7SRSinduced,
                  coef.K3QR, coef.K5QR, coef.K7QR)

#################### END OF STORING COEFFICIENTS ####################

#################### START OF STORING IMSEs ####################

# IMSEs
if(tau >= 4/7 || tau <= 3/7){
  IMSE.K7SRS <- IMSE.K7SRSapprox <- rep(NA, length = B)
}

IMSE = cbind(IMSE.QR, IMSE.K3, IMSE.K5, IMSE.K7, 
             IMSE.K3SRS, IMSE.K5SRS, IMSE.K7SRS,
             IMSE.K3approx, IMSE.K5approx, IMSE.K7approx, 
             IMSE.K3SRSapprox, IMSE.K5SRSapprox, IMSE.K7SRSapprox,
             IMSE.K3SRSinduced, IMSE.K5SRSinduced, IMSE.K7SRSinduced,
             IMSE.K3QR, IMSE.K5QR, IMSE.K7QR)
colnames(IMSE) <- c("QR", 
                    "k=3 R1", "k=3 R2", "k=3 R3", "k=3 R4",
                    "k=5 R1", "k=5 R2", "k=5 R3", "k=5 R4",
                    "k=7 R1", "k=7 R2", "k=7 R3", "k=7 R4", 
                    "k=3 SRS", "k=5 SRS", "k=7 SRS",
                    "k=3 R1 Approx", "k=3 R2 Approx", "k=3 R3 Approx", "k=3 R4 Approx",
                    "k=5 R1 Approx", "k=5 R2 Approx", "k=5 R3 Approx", "k=5 R4 Approx",
                    "k=7 R1 Approx", "k=7 R2 Approx", "k=7 R3 Approx", "k=7 R4 Approx", 
                    "k=3 SRS Approx", "k=5 SRS Approx", "k=7 SRS Approx",
                    "k=3 R1 IP", "k=3 R2 IP", "k=3 R3 IP", "k=3 R4 IP",
                    "k=5 R1 IP", "k=5 R2 IP", "k=5 R3 IP", "k=5 R4 IP",
                    "k=7 R1 IP", "k=7 R2 IP", "k=7 R3 IP", "k=7 R4 IP",
                    "k=3 R1 QR", "k=3 R2 QR", "k=3 R3 QR", "k=3 R4 QR",
                    "k=5 R1 QR", "k=5 R2 QR", "k=5 R3 QR", "k=5 R4 QR",
                    "k=7 R1 QR", "k=7 R2 QR", "k=7 R3 QR", "k=7 R4 QR")

#################### END OF STORING IMSEs ####################
