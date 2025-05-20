##### SIMULATION #####
##### High Dimensional Setting, constant variance, with outliers
##### May 22, 2024

#################### LOADING PACKAGES / FUNCTIONS ####################

# load packages
library(MASS)
library(ald)
library(quantreg)
library(glmnet)
library(CVXR)
source(here::here("Functions", "M-estimation Functions.R"))
source(here::here("Functions", "RSS Functions.R"))
source(here::here("Functions", "Other Functions.R"))


# Settings
if(T){
  B = 1000 # number of replications for simulation
  tau = 0.5 # quantile
  case = 3 # choose 1 through 3
}


#################### START OF GENERATING POPULATION ####################

### Population from Inesh Thesis

if(T){
  set.seed(400)
  
  N = 10000 # population size
  
  if(case == 1){
    
    # generate covariates
    Sigma = matrix(NA, nrow = 30, ncol = 30)
    for(i in 1:nrow(Sigma)){
      for(j in 1:ncol(Sigma)){
        Sigma[i, j] = 0.5^abs(i-j)
      }
    }
    X <- mvrnorm(N, mu = rep(0, times = 30), Sigma = Sigma)
    X <- cbind(1, X)
    
    # error
    sigma = 1
    error = rnorm(N, 0, sigma)
    
    # beta values
    beta = c(0, rgamma(10, shape = 1, rate = 10), rep(1, 5), 
             rgamma(10, shape = 1, rate = 10), rep(4, 5))
    
  } else if(case == 2){
    
    # generate covariates
    Sigma = matrix(NA, nrow = 40, ncol = 40)
    for(i in 1:nrow(Sigma)){
      for(j in 1:ncol(Sigma)){
        Sigma[i, j] = 0.5^abs(i-j)
      }
    }
    X <- mvrnorm(N, mu = rep(0, times = 40), Sigma = Sigma)
    X <- cbind(1, X)
    
    # error
    sigma = 1
    error = rnorm(N, 0, sigma)
    
    # beta values
    beta = c(0, rep(3, 5), 
             rgamma(10, shape = 1, rate = 10), rep(1, 10), 
             rgamma(10, shape = 1, rate = 10), rep(2, 5))
    
  } else if(case == 3){
    
    # generate covariates
    Sigma = matrix(NA, nrow = 40, ncol = 40)
    for(i in 1:nrow(Sigma)){
      for(j in 1:ncol(Sigma)){
        Sigma[i, j] = 0.5^abs(i-j)
      }
    }
    X <- mvrnorm(N, mu = rep(0, times = 40), Sigma = Sigma)
    X <- cbind(1, X)
    
    # error
    sigma = 1
    error = rnorm(N, 0, sigma)
    
    # beta values
    beta = c(0, rgamma(10, shape = 1, rate = 10), rep(5, 3), 
             rgamma(10, shape = 1, rate = 10), rep(1, 5),
             rgamma(10, shape = 1, rate = 10), rep(10, 2))
    
  }
  
  # response
  y <- as.vector( X %*% beta + error )
  # Add outliers
  y[1:(0.1*N)] <- rnorm((0.1*N), mean = 80)
  
  # Simulate rankers
  rho.val = c(0.9, 0.8, 0.7, 0.6) # selected correlation of rankers and response
  R1 = rho.val[1]*y + sqrt(1-rho.val[1]^2)*rnorm(N)
  R2 = rho.val[2]*y + sqrt(1-rho.val[2]^2)*rnorm(N)
  R3 = rho.val[3]*y + sqrt(1-rho.val[3]^2)*rnorm(N)
  R4 = rho.val[4]*y + sqrt(1-rho.val[4]^2)*rnorm(N)
  
  # Objects used in simulation
  features = as.matrix(X[,-1])
  population = cbind(y, features)
  rankers = cbind(R1, R2, R3, R4)
  d = length(beta)
  n = 100 # sample size
  
  lambdas = exp(seq(0, 5, length.out = 6))
  gamma = c(1, rep(0, d-1))
  
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
  coef.meanLASSO <- matrix(NA, nrow = B, ncol = d)
  coef.QRLASSO <- matrix(NA, nrow = B, ncol = d)
  
  # MRSS Loss results with LASSO
  coef.K3LASSO <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5LASSO <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7LASSO <- matrix(NA, nrow = B, ncol = 4*d)
  
  # Changing tau for SRS from induced population
  coef.K3SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K5SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  coef.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4*d)
  
  ### Store results for IMSEs
  
  # MR and QR to compare to
  IMSE.QRLASSO <- c()
  
  # Original MRSS Loss results
  IMSE.K3LASSO <- IMSE.K5LASSO <- IMSE.K7LASSO <- matrix(NA, nrow = B, ncol = 4)
  
  # Changing tau for SRS from induced population
  IMSE.K3SRSinduced <- IMSE.K5SRSinduced <- IMSE.K7SRSinduced <- matrix(NA, nrow = B, ncol = 4)
  
  ### coefficients across lambda (path of coefficients)
  
  path.meanLASSO <- path.QRLASSO <- 
    path.K3LASSO <- path.K5LASSO <- path.K7LASSO <- 
    path.K3SRSinducedLASSO <- path.K5SRSinducedLASSO <- path.K7SRSinducedLASSO <- list()
}
##### BEGIN REPLICATIONS
{
  start = Sys.time()
  
  for(b in 1:B){
    
    #### Setting aside SRS data
    
    ### SRS DATA (Training)
    #ind.srs <- sample(1:N, size = n)
    ind <- sample(1:N, size = 2*n) # when using validation set as well
    ind.srs <- ind[1:n]
    ind.val <- ind[(n+1):(2*n)]
    
    # training data (SRS)
    y.trainSRS <- y[ind.srs]
    features.trainSRS <- features[ind.srs,]
    X.trainSRS <- X[ind.srs,]
    
    # training data (SRS)
    y.val <- y[ind.val]
    features.val <- features[ind.val,]
    X.val <- X[ind.val,]
    
    ##### MEAN REG AND QR WITH LASSO
    
    cv.mean <- cv.glmnet(features.trainSRS, y.trainSRS, alpha = 1, nfolds = 5, lambda = lambdas)
    lasso.mean <- glmnet(features.trainSRS, y.trainSRS, alpha = 1, lambda = cv.mean$lambda.min)
    coef.meanLASSO[b,] <- as.vector(coef(lasso.mean))
    path.meanLASSO[[b]] <- coef(cv.mean, s = lambdas)
    
    cv.QR <- QR.LASSOcv(response = y.trainSRS, features = features.trainSRS, tau=tau, lambdas = lambdas)
    lasso.QR <- rq.fit.lasso(X.trainSRS, y.trainSRS, tau = tau, lambda = lambdas[which.min(cv.QR$errors)])
    coef.QRLASSO[b,] <- coef(lasso.QR)
    IMSE.QRLASSO[b] <- sum((X.trainSRS %*% beta_calc(tau) - X.trainSRS %*% coef(lasso.QR))^2)
    path.QRLASSO[[b]] <- QR.LASSOregpath(response = y.trainSRS, features = features.trainSRS, tau=0.5, lambdas = lambdas) # regularization path
    
    ##### FIT MEDIAN REGRESSION TO MRSS DATA (with each ranker)
    coef3LASSO <- coef5LASSO <- coef7LASSO <- c()
    coef3SRSinduced <- coef5SRSinduced <- coef7SRSinduced <- c()
    path.K3LASSO.loop <- path.K5LASSO.loop <- path.K7LASSO.loop <- NULL
    path.K3SRSinducedLASSO.loop <- path.K5SRSinducedLASSO.loop <- path.K7SRSinducedLASSO.loop <- NULL
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
        # using CVXR
        model.coef.3cv <- MedNS.LASSOvalidation(train.response = y.trainK3, train.features = features.trainK3,
                                                test.response = y.val, test.features = features.val,
                                                tau=tau, m=2, lambdas = lambdas)
        model.coef.3LASSO <- MedNS.LASSO(response = y.trainK3, features = features.trainK3,
                                         tau=tau, m=2, lambdas[which.min(model.coef.3cv$errors)])
        model.coef.3LASSO <- as.vector(model.coef.3LASSO)

        coef3LASSO <- c(coef3LASSO, model.coef.3LASSO)
        model.coef.3LASSO[1] <- model.coef.3LASSO[1] + (1-2*tau)
        IMSE.K3LASSO[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% model.coef.3LASSO)^2)
        #path.K3LASSO.loop <- cbind(path.K3LASSO.loop, MedNS.LASSOregpath(y.trainK3, features.trainK3, tau = tau, m=2, lambdas = lambdas))
      }
      ## k=5
      if(tau <= 3/5 && tau >= 2/5){
        # using CVXR
        model.coef.5cv <- MedNS.LASSOvalidation(train.response = y.trainK5, train.features = features.trainK5,
                                                test.response = y.val, test.features = features.val,
                                                tau=tau, m=3, lambdas = lambdas)
        model.coef.5LASSO <- MedNS.LASSO(response = y.trainK5, features = features.trainK5,
                                         tau=tau, m=3, lambda = lambdas[which.min(model.coef.5cv$errors)])
        model.coef.5LASSO <- as.vector(model.coef.5LASSO)

        coef5LASSO <- c(coef5LASSO, model.coef.5LASSO)
        model.coef.5LASSO[1] <- model.coef.5LASSO[1] + (1-2*tau)
        IMSE.K5LASSO[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% model.coef.5LASSO)^2)
        #path.K5LASSO.loop <- cbind(path.K5LASSO.loop, MedNS.LASSOregpath(y.trainK5, features.trainK5, tau = tau, m=3, lambdas = lambdas))
      }
      ## k=7
      if(tau <= 4/7 && tau >= 3/7){
        # using CVXR
        model.coef.7cv <- MedNS.LASSOvalidation(train.response = y.trainK7, train.features = features.trainK7,
                                                test.response = y.val, test.features = features.val,
                                                tau=tau, m=4, lambdas = lambdas)
        model.coef.7LASSO <- MedNS.LASSO(response = y.trainK7, features = features.trainK7,
                                         tau=tau, m=4, lambda = lambdas[which.min(model.coef.7cv$errors)])
        model.coef.7LASSO <- as.vector(model.coef.7LASSO)
      
        coef7LASSO <- c(coef7LASSO, model.coef.7LASSO)
        model.coef.7LASSO[1] <- model.coef.7LASSO[1] + (1-2*tau)
        IMSE.K7LASSO[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% model.coef.7LASSO)^2)
        #path.K7LASSO.loop <- cbind(path.K7LASSO.loop, MedNS.LASSOregpath(y.trainK7, features.trainK7, tau = tau, m=4, lambdas = lambdas))
      }
      
      ##### SRS from induced population Approach
      
      cv.3SRSinduced <- QR.LASSOcv(response = y.trainK3, features = features.trainK3, tau=tau.k.imperfect[r, 1], lambdas = lambdas)
      lasso.3SRSinduced <- rq.fit.lasso(X.trainK3, y.trainK3, tau = tau.k.imperfect[r, 1], lambda = lambdas[which.min(cv.3SRSinduced$errors)])
      
      cv.5SRSinduced <- QR.LASSOcv(response = y.trainK5, features = features.trainK5, tau=tau.k.imperfect[r, 2], lambdas = lambdas)
      lasso.5SRSinduced <- rq.fit.lasso(X.trainK5, y.trainK5, tau = tau.k.imperfect[r, 2], lambda = lambdas[which.min(cv.5SRSinduced$errors)])
      
      cv.7SRSinduced <- QR.LASSOcv(response = y.trainK7, features = features.trainK7, tau=tau.k.imperfect[r, 3], lambdas = lambdas)
      lasso.7SRSinduced <- rq.fit.lasso(X.trainK7, y.trainK7, tau = tau.k.imperfect[r, 3], lambda = lambdas[which.min(cv.7SRSinduced$errors)])
      
      # record coefficients
      coef3SRSinduced <- c(coef3SRSinduced, coef(lasso.3SRSinduced))
      coef5SRSinduced <- c(coef5SRSinduced, coef(lasso.5SRSinduced))
      coef7SRSinduced <- c(coef7SRSinduced, coef(lasso.7SRSinduced))
      
      IMSE.K3SRSinduced[b, r] <- sum((X.trainK3 %*% beta_calc(tau) - X.trainK3 %*% coef(lasso.3SRSinduced))^2)
      IMSE.K5SRSinduced[b, r] <- sum((X.trainK5 %*% beta_calc(tau) - X.trainK5 %*% coef(lasso.5SRSinduced))^2)
      IMSE.K7SRSinduced[b, r] <- sum((X.trainK7 %*% beta_calc(tau) - X.trainK7 %*% coef(lasso.7SRSinduced))^2)
      
      # regularization path
      #path.K3SRSinducedLASSO.loop <- cbind(path.K3SRSinducedLASSO.loop, QR.LASSOregpath(y.trainK3, features.trainK3, tau = tau.k.imperfect[r, 1], lambdas = lambdas))
      #path.K5SRSinducedLASSO.loop <- cbind(path.K5SRSinducedLASSO.loop, QR.LASSOregpath(y.trainK5, features.trainK5, tau = tau.k.imperfect[r, 2], lambdas = lambdas))
      #path.K7SRSinducedLASSO.loop <- cbind(path.K7SRSinducedLASSO.loop, QR.LASSOregpath(y.trainK7, features.trainK7, tau = tau.k.imperfect[r, 3], lambdas = lambdas))
    }
    
    ### record coefficients
    ## k=3
    if(tau <= 2/3 && tau >= 1/3){
      coef.K3LASSO[b, ] <- coef3LASSO
      #path.K3LASSO[[b]] <- path.K3LASSO.loop
    }
    ## k=5
    if(tau <= 3/5 && tau >= 2/5){
      coef.K5LASSO[b, ] <- coef5LASSO
      #path.K5LASSO[[b]] <- path.K5LASSO.loop
    }
    ## k=7
    if(tau <= 4/7 && tau >= 3/7){
      coef.K7LASSO[b, ] <- coef7LASSO
      #path.K7LASSO[[b]] <- path.K7LASSO.loop
    }
    
    coef.K3SRSinduced[b, ] <- coef3SRSinduced
    coef.K5SRSinduced[b, ] <- coef5SRSinduced
    coef.K7SRSinduced[b, ] <- coef7SRSinduced
    
    #path.K3SRSinducedLASSO[[b]] <- path.K3SRSinducedLASSO.loop
    #path.K5SRSinducedLASSO[[b]] <- path.K5SRSinducedLASSO.loop
    #path.K7SRSinducedLASSO[[b]] <- path.K7SRSinducedLASSO.loop
    
    if(b %% 10 == 0){ print(paste("iteration", sep = " = ", b)) } # counts every 5 iterations
    
  }
  
  end = Sys.time()
  
  end-start # time for simulation
  
}

#################### END OF SIMULATION ####################

#################### START OF STORING COEFFICIENTS ####################

coef.all <- cbind(coef.meanLASSO, coef.QRLASSO,
                  coef.K3LASSO, coef.K5LASSO, coef.K7LASSO,
                  coef.K3SRSinduced, coef.K5SRSinduced, coef.K7SRSinduced)

#write.csv(coef.all, file = here::here("Results", paste("TEST1000HighDim-Coef-Case", case, "-tau", tau, ".csv", sep = "")))

#################### END OF STORING COEFFICIENTS ####################

#################### START OF STORING IMSEs ####################



IMSE = cbind(IMSE.QRLASSO, IMSE.K3LASSO, IMSE.K5LASSO, IMSE.K7LASSO,
             IMSE.K3SRSinduced, IMSE.K5SRSinduced, IMSE.K7SRSinduced)
colnames(IMSE) <- c("QR", 
                    "k=3 R1", "k=3 R2", "k=3 R3", "k=3 R4",
                    "k=5 R1", "k=5 R2", "k=5 R3", "k=5 R4",
                    "k=7 R1", "k=7 R2", "k=7 R3", "k=7 R4",
                    "k=3 R1 SRS", "k=3 R2 SRS", "k=3 R3 SRS", "k=3 R4 SRS",
                    "k=5 R1 SRS", "k=5 R2 SRS", "k=5 R3 SRS", "k=5 R4 SRS",
                    "k=7 R1 SRS", "k=7 R2 SRS", "k=7 R3 SRS", "k=7 R4 SRS")

#write.csv(IMSE, file = here::here("Results", paste("TEST1000HighDim-IMSE-Case", case, "-tau", tau, ".csv", sep = "")))

#################### END OF STORING IMSEs ####################


#################### START OF STORING REGULARIZATION PATHS ####################

# Recorded paths in separate file due to computation time

avg.path.meanLASSO = Reduce("+", path.meanLASSO)/B
avg.path.QRLASSO = Reduce("+", path.QRLASSO)/B
avg.path.K3LASSO = Reduce("+", path.K3LASSO)/B
avg.path.K5LASSO = Reduce("+", path.K5LASSO)/B
avg.path.K7LASSO = Reduce("+", path.K7LASSO)/B
avg.path.K3SRSinducedLASSO = Reduce("+", path.K3SRSinducedLASSO)/B
avg.path.K5SRSinducedLASSO = Reduce("+", path.K5SRSinducedLASSO)/B
avg.path.K7SRSinducedLASSO = Reduce("+", path.K7SRSinducedLASSO)/B

avg.path = cbind(as.matrix(avg.path.meanLASSO), avg.path.QRLASSO,
                 avg.path.K3LASSO, avg.path.K5LASSO, avg.path.K7LASSO,
                 avg.path.K3SRSinducedLASSO, avg.path.K5SRSinducedLASSO, avg.path.K7SRSinducedLASSO)

#write.csv(avg.path, file = here::here("Results", paste("1000HighDim-PATH-Case", case, "-tau", tau, ".csv", sep = "")))

#################### END OF STORING REGULARIZATION PATHS ####################
