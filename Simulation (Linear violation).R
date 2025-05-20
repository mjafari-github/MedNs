##### SIMULATION #####
##### Linear Methods on Non-linear Population
##### May 16, 2025


#################### LOADING PACKAGES / FUNCTIONS ####################

# load packages
library(MASS) 
library(ald)
library(quantreg)
source(here::here("Functions", "M-estimation Functions.R"))
source(here::here("Functions", "RSS Functions.R"))
source(here::here("Functions", "Other Functions.R"))


# Settings that may be altered by master file
# Make FALSE if running from "Simulation Results" file
if(F){
  B = 500 # number of replications for simulation
  tau = 0.5 # quantile
  case = 3 # choose 1 through 3
}



# other functions for non-linear setting

# calculates X matrix for a given degree and single covariate x
X.poly <- function(x, degree){
  x.mat = matrix(NA, nrow = length(x), ncol = degree+1)
  x.mat[,1] = 1
  for(i in 1:degree){
    x.mat[,(i+1)] = x^i
  }
  return(x.mat)
}

# calculates ISE under each of the cases in the population
ISE.poly <- function(coef, tau){
  
  sqErr.poly <- function(x){
    # true relationship
    true = 2 + 3*x^2 + qnorm(tau, 0, 1)
    
    # predicted relationship
    pred = X.poly(x, length(coef)-1) %*% coef
    
    # function to integrate for ISE
    return((pred-true)^2)
  }
  
  # integration
  ISE = integrate(sqErr.poly, lower = -2, upper = 2)
  
  return(list("ISE" = ISE$value, "abs.error" = ISE$abs.error))
}

#################### START OF GENERATING POPULATION ####################

if(T){
  set.seed(8)
  
  N = 10000 # population size
  
  x <- rnorm(N, 0, 1) # covariate
  y <- 2 + 3*x^2 + rnorm(N) # response
  
  if(case == 2){ # 10% outliers in y direction
 # Add outliers
    out.data <- mvrnorm(0.1*N, mu=c(0, 40), Sigma=diag(0.25*c(1, 20)))
    x[1:(0.1*N)] <- out.data[,1]
    y[1:(0.1*N)] <- out.data[,2]
  } else if(case == 3){ # 1% outliers in y direction
    # Add outliers
    out.data <- mvrnorm(0.01*N, mu=c(6, 40), Sigma=diag(0.25*c(1, 20)))
    x[1:(0.01*N)] <- out.data[,1]
    y[1:(0.01*N)] <- out.data[,2]
  }
  
  # Simulate rankers
  rho.val = c(0.9, 0.8, 0.7, 0.6) # selected correlation of rankers and response
  R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
  R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
  R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
  R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  
  # Objects used in simulation
  X = cbind(1, x)
  features = as.matrix(X[,-1])
  population = cbind(y, features)
  rankers = cbind(R1, R2, R3, R4)
  n = 100 # sample size
  D = ncol(X)
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
  coef.mean <- matrix(NA, nrow = B, ncol = D)
  coef.QR <- matrix(NA, nrow = B, ncol = D)
  
  # MRSS Loss results with LASSO
  coef.K3 <- matrix(NA, nrow = B, ncol = 4*D)
  coef.K5 <- matrix(NA, nrow = B, ncol = 4*D)
  coef.K7 <- matrix(NA, nrow = B, ncol = 4*D)
  
  # Changing tau for SRS from induced population
  coef.K3SRSinduced <- matrix(NA, nrow = B, ncol = 4*D)
  coef.K5SRSinduced <- matrix(NA, nrow = B, ncol = 4*D)
  coef.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4*D)
  
  ### Store results for ISE
  
  # MR and QR to compare to
  ISE.mean <- ISE.QR <- c()
  
  # Original MRSS Loss results
  ISE.K3 <- ISE.K5 <- ISE.K7 <- matrix(NA, nrow = B, ncol = 4)
  
  # Changing tau for SRS from induced population
  ISE.K3SRSinduced <- ISE.K5SRSinduced <- ISE.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4)
}

##### BEGIN REPLICATIONS
{
  start = Sys.time()
  
  for(b in 1:B){
    
    #### Setting aside SRS data
    
    ### SRS DATA (Training)
    ind <- sample(1:N, size = n)
    
    # training data (SRS)
    y.trainSRS <- y[ind]
    features.trainSRS <- features[ind,]
    X.trainSRS <- X[ind,]
    
    ##### MEAN REG AND QR
    
    # OLS
    
    MeanReg <- lm(y.trainSRS ~ features.trainSRS)
    coef.mean[b,] <- MeanReg$coefficients
    ISE.mean[b] <- ISE.poly(MeanReg$coefficients, tau = tau)$ISE
    
    
    # QR
    
    QRcoef <- mrss.M.estimation(response = y.trainSRS, 
                                features = features.trainSRS, 
                                tau = tau, m = 1, max.diff = 1e-3) # NOTE: k=1 is equivalent to standard check function
    coef.QR[b,] <- QRcoef
    
    ISE.QR[b] <- ISE.poly(QRcoef, tau = tau)$ISE
    
    
    ##### FIT MEDIAN REGRESSION TO MRSS DATA (with each ranker)
    coef3 <- coef5 <- coef7 <- c()
    coef3SRSinduced <- coef5SRSinduced <- coef7SRSinduced <- c()
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
      y.trainK3 <- y[ind.trainK3]
      features.trainK3 <- features[ind.trainK3,]
      X.trainK3 <- X[ind.trainK3,]
      
      # training data (k=5)
      y.trainK5 <- y[ind.trainK5]
      features.trainK5 <- features[ind.trainK5,]
      X.trainK5 <- X[ind.trainK5,]
      
      # training data (k=7)
      y.trainK7 <- y[ind.trainK7]
      features.trainK7 <- features[ind.trainK7,]
      X.trainK7 <- X[ind.trainK7,]
      
      ##### MRSS REGRESSION
      
      ## k=3
      if(tau <= 2/3 && tau >= 1/3){
        model.coef.3 <- mrss.M.estimation(response = y.trainK3, 
                                          features = features.trainK3, 
                                          tau = tau, m=2)
        coef3 <- c(coef3, model.coef.3)
        model.coef.3[1] <- model.coef.3[1] + (1-2*tau)
        
        ISE.K3[b, r] <- ISE.poly(model.coef.3, tau = tau)$ISE
      }
      ## k=5
      if(tau <= 3/5 && tau >= 2/5){
        model.coef.5 <- mrss.M.estimation(response = y.trainK5, 
                                          features = features.trainK5, 
                                          tau = tau, m=3)
        coef5 <- c(coef5, model.coef.5)
        model.coef.5[1] <- model.coef.5[1] + (1-2*tau)

        ISE.K5[b, r] <- ISE.poly(model.coef.5, tau = tau)$ISE
      }
      ## k=7
      if(tau <= 4/7 && tau >= 3/7){
        model.coef.7 <- mrss.M.estimation(response = y.trainK7, 
                                          features = features.trainK7, 
                                          tau = tau, m=4)
        coef7 <- c(coef7, model.coef.7)
        model.coef.7[1] <- model.coef.7[1] + (1-2*tau)
        
        ISE.K7[b, r] <- ISE.poly(model.coef.7, tau = tau)$ISE
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
      
      ISE.K3SRSinduced[b, r] <- ISE.poly(model.3SRSinduced$coefficients, tau = tau)$ISE
      ISE.K5SRSinduced[b, r] <- ISE.poly(model.5SRSinduced$coefficients, tau = tau)$ISE
      ISE.K7SRSinduced[b, r] <- ISE.poly(model.7SRSinduced$coefficients, tau = tau)$ISE
    
    }
    
    ### record coefficients
    ## k=3
    if(tau <= 2/3 && tau >= 1/3){
      coef.K3[b, ] <- coef3
    }
    ## k=5
    if(tau <= 3/5 && tau >= 2/5){
      coef.K5[b, ] <- coef5
    }
    ## k=7
    if(tau <= 4/7 && tau >= 3/7){
      coef.K7[b, ] <- coef7
    }
    
    coef.K3SRSinduced[b, ] <- coef3SRSinduced
    coef.K5SRSinduced[b, ] <- coef5SRSinduced
    coef.K7SRSinduced[b, ] <- coef7SRSinduced
    
    if(b %% 100 == 0){ print(paste("iteration", sep = " = ", b)) } # counts every 5 iterations
    
  }
  
  end = Sys.time()
  
  end-start # time for simulation
  
}

#################### END OF SIMULATION ####################

#################### START OF STORING COEFFICIENTS ####################

coef.all <- cbind(coef.mean, coef.QR,
                  coef.K3, coef.K5, coef.K7,
                  coef.K3SRSinduced, coef.K5SRSinduced, coef.K7SRSinduced)

#write.csv(coef.all, file = here::here("Results", paste("5000Violation-Coef-Case", case, "-tau", tau, ".csv", sep = "")))

#################### END OF STORING COEFFICIENTS ####################

#################### START OF STORING ISEs ####################



ISE = cbind(ISE.mean, ISE.QR, ISE.K3, ISE.K5, ISE.K7,
            ISE.K3SRSinduced, ISE.K5SRSinduced, ISE.K7SRSinduced)
colnames(ISE) <- c("OLS", "QR", 
                   "k=3 R1", "k=3 R2", "k=3 R3", "k=3 R4",
                   "k=5 R1", "k=5 R2", "k=5 R3", "k=5 R4",
                   "k=7 R1", "k=7 R2", "k=7 R3", "k=7 R4",
                   "k=3 R1 SRS", "k=3 R2 SRS", "k=3 R3 SRS", "k=3 R4 SRS",
                   "k=5 R1 SRS", "k=5 R2 SRS", "k=5 R3 SRS", "k=5 R4 SRS",
                   "k=7 R1 SRS", "k=7 R2 SRS", "k=7 R3 SRS", "k=7 R4 SRS")

colMeans(ISE)

#write.csv(ISE, file = here::here("Results", paste("5000Violation-err-Case", case, "-tau", tau, ".csv", sep = "")))

#################### END OF STORING ISEs ####################



