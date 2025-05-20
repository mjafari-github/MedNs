##### M-estimation in standard quantile regression setting (SRS) #####

### LOSS FUNCTIONS

# standard loss function
rho = function(error, tau){
  w = c()
  for(i in 1:length(error)){
    w[i] = ifelse(error[i]>0, tau*error[i], (tau-1)*error[i])
  }
  return(w)
}

# Loss function for mrss regression
rho_mrss = function(error, m, tau){
  w = c()
  for(i in 1:length(error)){
    w[i] = ifelse(error[i]>0, 
                  tau*error[i] - ((m-1)/(2*m-1))*log((exp(tau*error[i])-(1-tau))/tau), 
                  (tau - 1)*error[i] - ((m-1)/(2*m-1))*log((exp((tau-1)*error[i])-tau)/(1-tau)))
  }
  return(w)
}

# Loss function for mrss regression approximation
rho_mrssAPPROX = function(error, m, tau){
  w = c()
  for(i in 1:length(error)){
    w[i] = ifelse(error[i]>0, tau*error[i], (tau - 1)*error[i]) - 
      ((m-1)/((2*m-1)))*abs(error[i]) + ((m-1)/(2*(2*m-1)))*(error[i]^2)
  }
  return(w)
}


####################################################################




### M-estimation functions
####################################################################


##### Functions for m-estimation in MRSS regression setting #####
##### Updated: May 2, 2024

## Consider MRSS with odd k
## m-estimation derivation is the same but now using a different loss function

# function is replicated from standard m-estimation function
# options are added for m=2k-1 to indicate set size
# loss function is changed to rho_mrss

mrss.M.estimation <- function(response, features, tau=0.5, m=1, max.diff = 1e-3, max.iter = 1000){
  
  # response: vector of response variables
  # features: data frame or matrix with features in columns
  # tau: quantile to estimate
  # m: function (k=2m-1 where k is the set size to guarantee odd set size)
  # max.diff: difference between old and new beta to close loop
  # max.iter: maximum iterations before the loop quits
  
  # organize data
  data = data.frame(cbind(response, features))
  
  # fit original model
  ls.model <- lm(response ~ . , data=data)
  beta.new <- ls.model$coefficients # initial beta
  
  # initiate values for while loop
  counter = 0
  beta.diff = max.diff+1
  
  while(beta.diff > max.diff && counter < max.iter){
    
    counter = counter + 1 # count iterations
    
    # store previous beta values
    beta.old = beta.new
    
    # calculate errors and weights
    Error = ls.model$residuals
    Error[round(Error, 10) == 0] = 0.00001 # avoid zero weights or negative weights from rounding issues
    
    w = rho_mrss(Error, m, tau)/Error^2 # calculate weights
    
    w[is.na(w)] = .001 # inf/inf -> error is large so weight should be small
    w = pmin(w, 10000) # make large or infinite values 10000
    w[!is.finite(w)] = 0.001 # make -Inf values small (weight near 0 but r can't calculate loss accurately)
    
    # fit new model
    ls.model <- lm(response ~ . , data=data, weights = w)
    beta.new <- ls.model$coefficients # new beta values
    
    # calculate difference between estimates to determine convergence
    beta.diff <- sum(abs(beta.new-beta.old)^2, na.rm = TRUE)/sum(beta.new^2, na.rm = TRUE)

    if(counter == max.iter){
      warning("Maximum iterations reached")
    }
  }
  return(beta.new)
}



####################################################################

##### Function for m-estimation in MRSS regression setting (APPROXIMATION) #####
##### Date: October 26, 2023

mrssAPPROX.M.estimation <- function(response, features, tau=0.5, m=1, max.diff = 1e-3, max.iter = 1000){
  
  # response: vector of response variables
  # features: data frame or matrix with features in columns
  # tau: quantile to estimate
  # m: function (k=2m-1 where k is the set size to guarantee odd set size)
  # max.diff: difference between old and new beta to close loop
  # max.iter: maximum iterations before the loop quits
  
  # organize data
  data = data.frame(cbind(response, features))
  
  # fit original model
  ls.model <- lm(response ~ . , data=data)
  beta.new <- ls.model$coefficients # initial beta
  
  # initiate values for while loop
  counter = 0
  beta.diff = max.diff+1
  
  while(beta.diff > max.diff && counter < max.iter){
    
    counter = counter + 1 # count iterations
    
    # store previous beta values
    beta.old = beta.new
    
    # calculate errors and weights
    Error = ls.model$residuals
    Error[round(Error, 10) == 0] = 0.00001 # avoid zero weights or negative weights from rounding issues
    
    w = rho_mrssAPPROX(Error, m, tau)/Error^2 # calculate weights
    
    w[is.na(w)] = .001 # inf/inf -> error is large so weight should be small
    w = pmin(w, 10000) # make large or infinite values 10000
    w[!is.finite(w)] = 0.001 # make -Inf values small (weight near 0 but r can't calculate loss accurately)
    
    # fit new model
    ls.model <- lm(response ~ . , data=data, weights = w)
    beta.new <- ls.model$coefficients # new beta values
    
    # calculate difference between estimates to determine convergence
    beta.diff <- sum(abs(beta.new-beta.old)^2, na.rm = TRUE)/sum(beta.new^2, na.rm = TRUE)
    
    if(counter == max.iter){
      warning("Maximum iterations reached")
    }
  }
  return(beta.new)
}










####################################################################
# LASSO FUNCTIONS
####################################################################

##### Functions for m-estimation in penalized MRSS regression #####
##### Updated: May 22, 2024

####################################################################

# MEDNS LASSO USING CVXR

MedNS.LASSO <- function(response, features, tau=0.5, m=1, lambda = 0, max.diff = 1e-3, max.iter = 100){
  
  # response: vector of response variables
  # features: data frame or matrix with features in columns
  # tau: quantile to estimate
  # m: function (k=2m-1 where k is the set size to guarantee odd set size)
  # lambda: LASSO penalty parameter (higher value = higher penalty)
  # max.diff: difference between old and new beta to close loop
  # max.iter: maximum iterations before the loop quits
  
  # organize data
  data = data.frame(cbind(response, features))
  
  # design matrix to be used for calculating residuals
  X <- as.matrix(cbind(1, features))
  
  # initial estimate from ols
  ls.model <- lm(response ~ . , data=data)
  beta.new <- ls.model$coefficients # initial beta
  
  # initiate values for while loop
  counter = 0
  beta.diff = max.diff+1
  
  while(beta.diff > max.diff && counter < max.iter){
    
    counter = counter + 1 # count iterations
    
    # store previous beta values
    beta.old = beta.new
    
    # obtain estimate of residuals
    error = response - X %*% beta.new
    error[round(error, 10) == 0] = 0.00001 # avoid zero weights or negative weights from rounding issues
    
    # calculate weights
    w = rho_mrss(error, m, tau)/error^2
    
    w[is.na(w)] = .001 # inf/inf -> error is large so weight should be small
    w = pmin(w, 10000) # make large or infinite values 10000
    w[!is.finite(w)] = 0.001 # make -Inf values small (weight near 0 but r can't calculate loss accurately)
    
    # optimize objective function for new parameters
    # CVXR
    beta.var <- Variable(d)
    obj <- sum(w*(response - X %*% beta.var)^2) + lambda * p_norm(beta.var[-1], 1)
    prob <- Problem(Minimize(obj))
    result <- solve(prob)
    beta.new <- result$getValue(beta.var) # extract estimates
    
    # calculate difference between estimates to determine convergence
    beta.diff <- sum(abs(beta.new-beta.old)^2, na.rm = TRUE)/sum(beta.new^2, na.rm = TRUE)
    
    if(counter == max.iter){warning("Maximum iterations reached")}
  }
  return(beta.new)
}

####################################################################


# VALIDATION FUNCTION FOR LASSO
####################################################################

# FUNCTION USING VALIDATION SET

MedNS.LASSOvalidation <- function(train.response, train.features, test.response, test.features, tau=0.5, m=1, lambdas = seq(1, 7, by=1)){
  
  IMSE.lambda = c() # errors for each lambda
  for(l in 1:length(lambdas)){
    # design matrix to be used for calculating residuals
    X <- cbind(1, test.features)
    
    # fit model using glmnet/cvxr
    model.estimates = MedNS.LASSO(response = train.response, features = train.features, tau = tau, m = m, lambda = lambdas[l], max.diff = 1e-2, max.iter = 20)   
    IMSE.lambda[l] = sum((X %*% (beta + gamma*qnorm(tau)) - X %*% model.estimates)^2) # test IMSE
  }
  
  # returns errors for different lambda
  return(list("lambdas" = lambdas, "errors" = IMSE.lambda))
}



# CROSS VALIDATION FUNCTION USING CVXR / GLMNET
# NOT USED IN THESIS TO REDUCE TIME

MedNS.LASSOcv <- function(response, features, tau=0.5, m=1, num.folds = 5, lambdas = seq(1, 7, by=1)){
  
  # make folds for cross validation
  elem = sample(1:NROW(features), NROW(features))
  folds = matrix(elem, ncol = num.folds)
  
  # start of cross validation
  IMSE.lambda = c()
  for(l in 1:length(lambdas)){
    
    IMSE.vec = c() # holds the error for each fold
    for(i in 1:num.folds)
    {
      # separating test and train for each fold
      test.elem = folds[, i]
      train.response = response[-test.elem]
      train.features = features[-test.elem, ]
      test.response = response[test.elem]
      test.features = features[test.elem, ]
      
      # design matrix to be used for calculating residuals
      X <- cbind(1, test.features)
      
      # fit model using glmnet/cvxr
      model.estimates = MedNS.LASSO(response = train.response, features = train.features, tau = tau, m = m, lambda = lambdas[l], max.diff = 1e-2, max.iter = 20)   
      IMSE.vec[i] = sum((X %*% (beta + gamma*qnorm(tau)) - X %*% model.estimates)^2) # test IMSE
    }
    IMSE.lambda[l] = mean(IMSE.vec)
  }
  # returns errors for different lambda
  # use m-estimation on whole sample to get final results
  return(list("lambdas" = lambdas, "errors" = IMSE.lambda))
}



# CROSS VALIDATION FUNCTION FOR QR

QR.LASSOcv <- function(response, features, tau=0.5, num.folds = 5, lambdas = seq(1, 7, by=1)){
  
  # make folds for cross validation
  elem = sample(1:NROW(features), NROW(features))
  folds = matrix(elem, ncol = num.folds)
  
  # start of cross validation
  IMSE.lambda = c()
  for(l in 1:length(lambdas)){
    
    IMSE.vec = c() # holds the error for each fold
    for(i in 1:num.folds)
    {
      # separating test and train for each fold
      test.elem = folds[, i]
      train.response = response[-test.elem]
      train.features = features[-test.elem, ]
      test.response = response[test.elem]
      test.features = features[test.elem, ]
      
      # design matrix to be used for calculating residuals
      X <- cbind(1, test.features)
      
      # fit model
      model = rq.fit.lasso(x = cbind(1, train.features), y = train.response, tau = tau, lambda =lambdas[l])
      IMSE.vec[i] = sum((X %*% (beta + gamma*qnorm(tau)) - X %*% coef(model))^2) # test IMSE
    }
    IMSE.lambda[l] = mean(IMSE.vec)
  }
  # returns errors for different lambda
  # use m-estimation on whole sample to get final results
  return(list("lambdas" = lambdas, "errors" = IMSE.lambda))
}

####################################################################



####################################################################
# FUNCTIONS TO RECORD REGULARIZATION PATH
####################################################################

# FOR QR
QR.LASSOregpath <- function(response, features, tau=0.5, lambdas = seq(1, 7, by=1)){
  coef.mat = matrix(NA, nrow = d, ncol = length(lambdas))
  for(l in 1:length(lambdas)){
    model = rq.fit.lasso(x = cbind(1, features), y = response, tau = tau, lambda = lambdas[l])
    coef.mat[,l] = coef(model)
  }
  return(coef.mat)
}

# FOR MEDNS REGRESSION
MedNS.LASSOregpath <- function(response, features, tau=0.5, m=1, lambdas = seq(1, 7, by=1)){
  coef.mat = matrix(NA, nrow = d, ncol = length(lambdas))
  for(l in 1:length(lambdas)){
    # GLMNET / CVXR
    model = MedNS.LASSO(response, features, tau, m, lambda = lambdas[l], max.diff = 1e-2, max.iter = 50)
    coef.mat[,l] = as.vector(model)
  }
  return(coef.mat)
}

####################################################################














# NON-LINEAR FUNCTIONS
####################################################################

##### Non-linear M-estimation function #####
##### Updated: May 22, 2024

MedNS.POLY <- function(response, features, tau=0.5, m=1, degree = 1, max.diff = 1e-3, max.iter = 1000){
  
  # response: vector of response variables
  # features: data frame or matrix with features in columns
  # tau: quantile to estimate
  # m: function (k=2m-1 where k is the set size to guarantee odd set size)
  # max.diff: difference between old and new beta to close loop
  # max.iter: maximum iterations before the loop quits
  # degree: degree of polynomial
  
  # organize data
  data = data.frame(cbind(response, features))
  X = cbind(1, features)
  
  # fit original model
  ls.model <- lm(response ~ poly(features, degree = degree) , data=data)
  beta.new <- ls.model$coefficients # initial beta
  
  # initiate values for while loop
  counter = 0
  beta.diff = max.diff+1
  
  while(beta.diff > max.diff && counter < max.iter){
    
    counter = counter + 1 # count iterations
    
    # store previous beta values
    beta.old = beta.new
    
    # calculate errors and weights
    Error = ls.model$residuals
    Error[round(Error, 10) == 0] = 0.00001 # avoid zero weights or negative weights from rounding issues
    
    w = rho_mrss(Error, m, tau)/Error^2 # calculate weights
    
    w[is.na(w)] = .001 # inf/inf -> error is large so weight should be small
    w = pmin(w, 10000) # make large or infinite values 10000
    w[!is.finite(w)] = 0.001 # make -Inf values small (weight near 0 but r can't calculate loss accurately)
    
    # fit new model
    ls.model <- lm(response ~ poly(features, degree = degree), data=data, weights = w)
    beta.new <- ls.model$coefficients # new beta values
    
    # calculate difference between estimates to determine convergence
    beta.diff <- sum(abs(beta.new-beta.old)^2, na.rm = TRUE)/sum(beta.new^2, na.rm = TRUE)
    
    if(counter == max.iter){
      warning("Maximum iterations reached")
    }
  }
  return(beta.new)
}

##### Non-linear M-estimation function #####
##### Updated: May 22, 2024

QR.POLY <- function(response, features, tau=0.5, degree = 1, max.diff = 1e-3, max.iter = 1000){
  
  # response: vector of response variables
  # features: data frame or matrix with features in columns
  # tau: quantile to estimate
  # max.diff: difference between old and new beta to close loop
  # max.iter: maximum iterations before the loop quits
  # degree: degree of polynomial
  
  # organize data
  data = data.frame(cbind(response, features))
  X = cbind(1, features)
  
  # fit original model
  ls.model <- lm(response ~ poly(features, degree = degree) , data=data)
  beta.new <- ls.model$coefficients # initial beta
  
  # initiate values for while loop
  counter = 0
  beta.diff = max.diff+1
  
  while(beta.diff > max.diff && counter < max.iter){
    
    counter = counter + 1 # count iterations
    
    # store previous beta values
    beta.old = beta.new
    
    # calculate errors and weights
    Error = ls.model$residuals
    Error[round(Error, 10) == 0] = 0.00001 # avoid zero weights or negative weights from rounding issues
    
    w = rho(Error, tau)/Error^2 # calculate weights
    
    w[is.na(w)] = .001 # inf/inf -> error is large so weight should be small
    w = pmin(w, 10000) # make large or infinite values 10000
    w[!is.finite(w)] = 0.001 # make -Inf values small (weight near 0 but r can't calculate loss accurately)
    
    # fit new model
    ls.model <- lm(response ~ poly(features, degree = degree), data=data, weights = w)
    beta.new <- ls.model$coefficients # new beta values
    
    # calculate difference between estimates to determine convergence
    beta.diff <- sum(abs(beta.new-beta.old)^2, na.rm = TRUE)/sum(beta.new^2, na.rm = TRUE)
    
    if(counter == max.iter){
      warning("Maximum iterations reached")
    }
  }
  return(beta.new)
}


# VALIDATION FUNCTION FOR NON-LINEAR
####################################################################

# FUNCTION USING VALIDATION SET

MedNS.POLYvalidation <- function(train.response, train.features, test.response, test.features, tau=0.5, m=1, degrees = seq(1, 5, by=1)){
  
  err.deg = c() # errors for each lambda
  for(i in 1:length(degrees)){
    # fit model using glmnet/cvxr
    model.estimates = MedNS.POLY(response = train.response, features = train.features, tau, m, degree = degrees[i], max.diff = 1e-3, max.iter = 50)
    model.fit = X.poly(test.features, degrees[i]) %*% as.vector(model.estimates) # predicted values
    err.deg[i] = sum((test.response - model.fit)^2) # test error
  }
  
  # returns errors for different lambda
  return(list("degrees" = degrees, "errors" = err.deg))
}

MeanReg.POLYvalidation <- function(train.response, train.features, test.response, test.features, degrees = seq(1, 5, by=1)){
  
  err.deg = c() # errors for each lambda
  for(i in 1:length(degrees)){
    model = lm(train.response ~ poly(train.features, degree = degrees[i]))
    model.fit = predict(model, newdata = data.frame(test.features)) #X.poly(test.features, degrees[i]) %*% as.vector(model.estimates) # predicted values
    err.deg[i] = sum((test.response - model.fit)^2) # test error
  }
  
  # returns errors for different lambda
  return(list("degrees" = degrees, "errors" = err.deg))
}

QR.POLYvalidation <- function(train.response, train.features, test.response, test.features, tau=0.5, degrees = seq(1, 5, by=1)){
  
  err.deg = c() # errors for each lambda
  for(i in 1:length(degrees)){
    
    # fit model using rq
    #model = rq(train.response ~ poly(train.features, degree = degrees[i]), tau = tau)
    #model.fit = predict(model, newdata = data.frame(test.features)) #X.poly(test.features, degrees[i]) %*% as.vector(model.estimates) # predicted values
    
    # fit model using m-estimation
    model.estimates = QR.POLY(response = train.response, features = train.features, tau, degree = degrees[i], max.diff = 1e-3, max.iter = 50)
    model.fit = X.poly(test.features, degrees[i]) %*% as.vector(model.estimates) # predicted values
    
    err.deg[i] = sum((test.response - model.fit)^2) # test error
  }
  
  # returns errors for different lambda
  return(list("degrees" = degrees, "errors" = err.deg))
}

# CROSS VALIDATION FUNCTION (NOT USED)

MedNS.POLYcv <- function(response, features, tau=0.5, m=1, num.folds = 5, degrees = seq(1, 5, by=1)){
  
  # make folds for cross validation
  elem = sample(1:NROW(features), NROW(features))
  folds = matrix(elem, ncol = num.folds)
  
  # start of cross validation
  error.p = c()
  for(l in 1:length(degrees)){
    
    error = c() # holds the error for each fold
    for(i in 1:num.folds)
    {
      # separating test and train for each fold
      test.elem = folds[, i]
      train.response = response[-test.elem]
      train.features = features[-test.elem, ]
      test.response = response[test.elem]
      test.features = features[test.elem, ]
      
      # design matrix to be used for calculating residuals
      X <- cbind(1, test.features)
      
      # fit model
      model.estimates = MedNS.POLY(train.response, train.features, tau, m, degree = degrees[l], max.diff = 1e-3, max.iter = 50)
      model.fit = X %*% model.estimates # predicted values
      error[i] = sum((test.response - model.fit)^2) # test error
    }
    error.p[l] = mean(error)
  }
  # returns errors for different degrees
  # use m-estimation on whole sample to get final results
  return(list("degrees" = degrees, "errors" = error.p))
}





