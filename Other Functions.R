##### OTHER FUNCTIONS
##### Any other functions used across multiple scripts

# min-max normalization
normalize <- function(x) {
  return((x- min(x)) /(max(x)-min(x)))
}

# calculates true beta to compare estimates to
# beta, gamma, and case should be predefined
beta_calc <- function(tau){
    beta + gamma*qnorm(tau)
}

##### TAU VALUES FOR SRS APPROACH #####
### PERFECT RANKING
tau.k.perfect <- function(k, tau){
  m<- (k+1)/2
  z<- qbeta(tau, m, m)
  return(z)
}

##### FOR RESULTS

Variance.calc <- function(Coefficients, tau){
  
  # Coefficients for each method
  coef.mean <- Coefficients[[match(tau, taus)]][,1:d]
  coef.QR <- Coefficients[[match(tau, taus)]][,(d+1):(2*d)]
  coef.K3 <- Coefficients[[match(tau, taus)]][,(2*d+1):(6*d)]
  coef.K5 <- Coefficients[[match(tau, taus)]][,(6*d+1):(10*d)]
  coef.K7 <- Coefficients[[match(tau, taus)]][,(10*d+1):(14*d)]
  coef.K3SRS <- Coefficients[[match(tau, taus)]][,(14*d+1):(15*d)]
  coef.K5SRS <- Coefficients[[match(tau, taus)]][,(15*d+1):(16*d)]
  coef.K7SRS <- Coefficients[[match(tau, taus)]][,(16*d+1):(17*d)]
  coef.K3approx <- Coefficients[[match(tau, taus)]][,(17*d+1):(21*d)]
  coef.K5approx <- Coefficients[[match(tau, taus)]][,(21*d+1):(25*d)]
  coef.K7approx <- Coefficients[[match(tau, taus)]][,(25*d+1):(29*d)]
  coef.K3SRSapprox <- Coefficients[[match(tau, taus)]][,(29*d+1):(30*d)]
  coef.K5SRSapprox <- Coefficients[[match(tau, taus)]][,(30*d+1):(31*d)]
  coef.K7SRSapprox <- Coefficients[[match(tau, taus)]][,(31*d+1):(32*d)]
  coef.K3SRSinduced <- Coefficients[[match(tau, taus)]][,(32*d+1):(36*d)]
  coef.K5SRSinduced <- Coefficients[[match(tau, taus)]][,(36*d+1):(40*d)]
  coef.K7SRSinduced <- Coefficients[[match(tau, taus)]][,(40*d+1):(44*d)]
  
  coef.K3QR <- Coefficients[[match(tau, taus)]][,(44*d+1):(48*d)]
  coef.K5QR <- Coefficients[[match(tau, taus)]][,(48*d+1):(52*d)]
  coef.K7QR <- Coefficients[[match(tau, taus)]][,(52*d+1):(56*d)]
  
  # variance 
  var.coef.mean <- apply(coef.mean, 2, var)
  var.coef.QR <- apply(coef.QR, 2, var)
  var.coef.K3 <- matrix(apply(coef.K3, 2, var), ncol = 4)
  var.coef.K5 <- matrix(apply(coef.K5, 2, var), ncol = 4)
  var.coef.K7 <- matrix(apply(coef.K7, 2, var), ncol = 4)
  var.coef.K3SRS <- apply(coef.K3SRS, 2, var)
  var.coef.K5SRS <- apply(coef.K5SRS, 2, var)
  var.coef.K7SRS <- apply(coef.K7SRS, 2, var)
  var.coef.K3approx <- matrix(apply(coef.K3approx, 2, var), ncol = 4)
  var.coef.K5approx <- matrix(apply(coef.K5approx, 2, var), ncol = 4)
  var.coef.K7approx <- matrix(apply(coef.K7approx, 2, var), ncol = 4)
  var.coef.K3SRSapprox <- apply(coef.K3SRSapprox, 2, var)
  var.coef.K5SRSapprox <- apply(coef.K5SRSapprox, 2, var)
  var.coef.K7SRSapprox <- apply(coef.K7SRSapprox, 2, var)
  var.coef.K3SRSinduced <- matrix(apply(coef.K3SRSinduced, 2, var), ncol = 4)
  var.coef.K5SRSinduced <- matrix(apply(coef.K5SRSinduced, 2, var), ncol = 4)
  var.coef.K7SRSinduced <- matrix(apply(coef.K7SRSinduced, 2, var), ncol = 4)
  var.coef.K3QR <- matrix(apply(coef.K3QR, 2, var), ncol = 4)
  var.coef.K5QR <- matrix(apply(coef.K5QR, 2, var), ncol = 4)
  var.coef.K7QR <- matrix(apply(coef.K7QR, 2, var), ncol = 4)
  
  var.coef = rbind(var.coef.mean, var.coef.QR, 
                   t(var.coef.K3), t(var.coef.K5), t(var.coef.K7), 
                   var.coef.K3SRS, var.coef.K5SRS, var.coef.K7SRS,
                   t(var.coef.K3approx), t(var.coef.K5approx), t(var.coef.K7approx), 
                   var.coef.K3SRSapprox, var.coef.K5SRSapprox, var.coef.K7SRSapprox,
                   t(var.coef.K3SRSinduced), t(var.coef.K5SRSinduced), t(var.coef.K7SRSinduced),
                   t(var.coef.K3QR), t(var.coef.K5QR), t(var.coef.K7QR))
  
  return(var.coef)
}

Bias.calc <- function(Coefficients, tau){
  
  # Coefficients for each method
  coef.mean <- Coefficients[[match(tau, taus)]][,1:d]
  coef.QR <- Coefficients[[match(tau, taus)]][,(d+1):(2*d)]
  coef.K3 <- Coefficients[[match(tau, taus)]][,(2*d+1):(6*d)]
  coef.K5 <- Coefficients[[match(tau, taus)]][,(6*d+1):(10*d)]
  coef.K7 <- Coefficients[[match(tau, taus)]][,(10*d+1):(14*d)]
  coef.K3SRS <- Coefficients[[match(tau, taus)]][,(14*d+1):(15*d)]
  coef.K5SRS <- Coefficients[[match(tau, taus)]][,(15*d+1):(16*d)]
  coef.K7SRS <- Coefficients[[match(tau, taus)]][,(16*d+1):(17*d)]
  coef.K3approx <- Coefficients[[match(tau, taus)]][,(17*d+1):(21*d)]
  coef.K5approx <- Coefficients[[match(tau, taus)]][,(21*d+1):(25*d)]
  coef.K7approx <- Coefficients[[match(tau, taus)]][,(25*d+1):(29*d)]
  coef.K3SRSapprox <- Coefficients[[match(tau, taus)]][,(29*d+1):(30*d)]
  coef.K5SRSapprox <- Coefficients[[match(tau, taus)]][,(30*d+1):(31*d)]
  coef.K7SRSapprox <- Coefficients[[match(tau, taus)]][,(31*d+1):(32*d)]
  coef.K3SRSinduced <- Coefficients[[match(tau, taus)]][,(32*d+1):(36*d)]
  coef.K5SRSinduced <- Coefficients[[match(tau, taus)]][,(36*d+1):(40*d)]
  coef.K7SRSinduced <- Coefficients[[match(tau, taus)]][,(40*d+1):(44*d)]

  coef.K3QR <- Coefficients[[match(tau, taus)]][,(44*d+1):(48*d)]
  coef.K5QR <- Coefficients[[match(tau, taus)]][,(48*d+1):(52*d)]
  coef.K7QR <- Coefficients[[match(tau, taus)]][,(52*d+1):(56*d)]
  
  # Average coefficients
  avg.coef.mean <- apply(coef.mean, 2, mean)
  avg.coef.QR <- apply(coef.QR, 2, mean)
  avg.coef.K3 <- matrix(apply(coef.K3, 2, mean), ncol = 4)
  avg.coef.K5 <- matrix(apply(coef.K5, 2, mean), ncol = 4)
  avg.coef.K7 <- matrix(apply(coef.K7, 2, mean), ncol = 4)
  avg.coef.K3SRS <- apply(coef.K3SRS, 2, mean)
  avg.coef.K5SRS <- apply(coef.K5SRS, 2, mean)
  avg.coef.K7SRS <- apply(coef.K7SRS, 2, mean)
  avg.coef.K3approx <- matrix(apply(coef.K3approx, 2, mean), ncol = 4)
  avg.coef.K5approx <- matrix(apply(coef.K5approx, 2, mean), ncol = 4)
  avg.coef.K7approx <- matrix(apply(coef.K7approx, 2, mean), ncol = 4)
  avg.coef.K3SRSapprox <- apply(coef.K3SRSapprox, 2, mean)
  avg.coef.K5SRSapprox <- apply(coef.K5SRSapprox, 2, mean)
  avg.coef.K7SRSapprox <- apply(coef.K7SRSapprox, 2, mean)
  avg.coef.K3SRSinduced <- matrix(apply(coef.K3SRSinduced, 2, mean), ncol = 4)
  avg.coef.K5SRSinduced <- matrix(apply(coef.K5SRSinduced, 2, mean), ncol = 4)
  avg.coef.K7SRSinduced <- matrix(apply(coef.K7SRSinduced, 2, mean), ncol = 4)
  avg.coef.K3QR <- matrix(apply(coef.K3QR, 2, mean), ncol = 4)
  avg.coef.K5QR <- matrix(apply(coef.K5QR, 2, mean), ncol = 4)
  avg.coef.K7QR <- matrix(apply(coef.K7QR, 2, mean), ncol = 4)
  
  # bias
  bias.coef.mean <- avg.coef.mean - beta_tau
  bias.coef.QR <- avg.coef.QR - beta_tau
  bias.coef.K3 <- avg.coef.K3 - beta_tau
  bias.coef.K5 <- avg.coef.K5 - beta_tau
  bias.coef.K7 <- avg.coef.K7 - beta_tau
  bias.coef.K3SRS <- avg.coef.K3SRS - beta_tau
  bias.coef.K5SRS <- avg.coef.K5SRS - beta_tau
  bias.coef.K7SRS <- avg.coef.K7SRS - beta_tau
  bias.coef.K3approx <- avg.coef.K3approx - beta_tau
  bias.coef.K5approx <- avg.coef.K5approx - beta_tau
  bias.coef.K7approx <- avg.coef.K7approx - beta_tau
  bias.coef.K3SRSapprox <- avg.coef.K3SRSapprox - beta_tau
  bias.coef.K5SRSapprox <- avg.coef.K5SRSapprox - beta_tau
  bias.coef.K7SRSapprox <- avg.coef.K7SRSapprox - beta_tau
  bias.coef.K3SRSinduced <- avg.coef.K3SRSinduced - beta_tau
  bias.coef.K5SRSinduced <- avg.coef.K5SRSinduced - beta_tau
  bias.coef.K7SRSinduced <- avg.coef.K7SRSinduced - beta_tau
  bias.coef.K3QR <- avg.coef.K3QR - beta_tau
  bias.coef.K5QR <- avg.coef.K5QR - beta_tau
  bias.coef.K7QR <- avg.coef.K7QR - beta_tau
  
  bias.coef = rbind(bias.coef.mean, bias.coef.QR, 
                    t(bias.coef.K3), t(bias.coef.K5), t(bias.coef.K7), 
                    bias.coef.K3SRS, bias.coef.K5SRS, bias.coef.K7SRS,
                    t(bias.coef.K3approx), t(bias.coef.K5approx), t(bias.coef.K7approx), 
                    bias.coef.K3SRSapprox, bias.coef.K5SRSapprox, bias.coef.K7SRSapprox,
                    t(bias.coef.K3SRSinduced), t(bias.coef.K5SRSinduced), t(bias.coef.K7SRSinduced),
                    t(bias.coef.K3QR), t(bias.coef.K5QR), t(bias.coef.K7QR))
  return(bias.coef)
}

RB.calc <- function(Coefficients, tau){
  
  # Coefficients for each method
  coef.mean <- Coefficients[[match(tau, taus)]][,1:d]
  coef.QR <- Coefficients[[match(tau, taus)]][,(d+1):(2*d)]
  coef.K3 <- Coefficients[[match(tau, taus)]][,(2*d+1):(6*d)]
  coef.K5 <- Coefficients[[match(tau, taus)]][,(6*d+1):(10*d)]
  coef.K7 <- Coefficients[[match(tau, taus)]][,(10*d+1):(14*d)]
  coef.K3SRS <- Coefficients[[match(tau, taus)]][,(14*d+1):(15*d)]
  coef.K5SRS <- Coefficients[[match(tau, taus)]][,(15*d+1):(16*d)]
  coef.K7SRS <- Coefficients[[match(tau, taus)]][,(16*d+1):(17*d)]
  coef.K3approx <- Coefficients[[match(tau, taus)]][,(17*d+1):(21*d)]
  coef.K5approx <- Coefficients[[match(tau, taus)]][,(21*d+1):(25*d)]
  coef.K7approx <- Coefficients[[match(tau, taus)]][,(25*d+1):(29*d)]
  coef.K3SRSapprox <- Coefficients[[match(tau, taus)]][,(29*d+1):(30*d)]
  coef.K5SRSapprox <- Coefficients[[match(tau, taus)]][,(30*d+1):(31*d)]
  coef.K7SRSapprox <- Coefficients[[match(tau, taus)]][,(31*d+1):(32*d)]
  coef.K3SRSinduced <- Coefficients[[match(tau, taus)]][,(32*d+1):(36*d)]
  coef.K5SRSinduced <- Coefficients[[match(tau, taus)]][,(36*d+1):(40*d)]
  coef.K7SRSinduced <- Coefficients[[match(tau, taus)]][,(40*d+1):(44*d)]
  
  coef.K3QR <- Coefficients[[match(tau, taus)]][,(44*d+1):(48*d)]
  coef.K5QR <- Coefficients[[match(tau, taus)]][,(48*d+1):(52*d)]
  coef.K7QR <- Coefficients[[match(tau, taus)]][,(52*d+1):(56*d)]
  
  # MEDIANS
  med.coef.mean <- apply(coef.mean, 2, median)
  med.coef.QR <- apply(coef.QR, 2, median)
  med.coef.K3 <- matrix(apply(coef.K3, 2, median), ncol = 4)
  med.coef.K5 <- matrix(apply(coef.K5, 2, median), ncol = 4)
  med.coef.K7 <- matrix(apply(coef.K7, 2, median), ncol = 4)
  med.coef.K3SRS <- apply(coef.K3SRS, 2, median)
  med.coef.K5SRS <- apply(coef.K5SRS, 2, median)
  med.coef.K7SRS <- apply(coef.K7SRS, 2, median)
  med.coef.K3approx <- matrix(apply(coef.K3approx, 2, median), ncol = 4)
  med.coef.K5approx <- matrix(apply(coef.K5approx, 2, median), ncol = 4)
  med.coef.K7approx <- matrix(apply(coef.K7approx, 2, median), ncol = 4)
  med.coef.K3SRSapprox <- apply(coef.K3SRSapprox, 2, median)
  med.coef.K5SRSapprox <- apply(coef.K5SRSapprox, 2, median)
  med.coef.K7SRSapprox <- apply(coef.K7SRSapprox, 2, median)
  med.coef.K3SRSinduced <- matrix(apply(coef.K3SRSinduced, 2, median), ncol = 4)
  med.coef.K5SRSinduced <- matrix(apply(coef.K5SRSinduced, 2, median), ncol = 4)
  med.coef.K7SRSinduced <- matrix(apply(coef.K7SRSinduced, 2, median), ncol = 4)
  med.coef.K3QR <- matrix(apply(coef.K3QR, 2, median), ncol = 4)
  med.coef.K5QR <- matrix(apply(coef.K5QR, 2, median), ncol = 4)
  med.coef.K7QR <- matrix(apply(coef.K7QR, 2, median), ncol = 4)
  
  ### ROBUST BIAS
  Medians <- rbind(med.coef.mean, med.coef.QR, 
                   t(med.coef.K3), t(med.coef.K5), t(med.coef.K7),
                   med.coef.K3SRS, med.coef.K5SRS, med.coef.K7SRS,
                   t(med.coef.K3approx), t(med.coef.K5approx), t(med.coef.K7approx),
                   med.coef.K3SRSapprox, med.coef.K5SRSapprox, med.coef.K7SRSapprox,
                   t(med.coef.K3SRSinduced), t(med.coef.K5SRSinduced), t(med.coef.K7SRSinduced),
                   t(med.coef.K3QR), t(med.coef.K5QR), t(med.coef.K7QR))
  Robust.Bias <- sweep(Medians, 2, beta_tau)
  
  return(Robust.Bias)
}
