##### CHAPTER 3: SIMULATION 1 PLOTS

# load packages
library(MASS)
library(ald)
library(quantreg)
library(xtable)
source("M-estimation Functions.R")
source("Other Functions.R")
source("RSS Functions.R")


case = 1 # choose 1, 5, 6 for different error setting in paper
palette(grey.colors(10, start=0, end = 0.8))

#################### START OF RESULTS CALCULATION ####################
{
  beta = c(0, 1, 1, 1) # parameter values
  gamma = c(1, 0, 0, 0)
  d = length(beta)
  
  taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
  model.names <- c("MeanReg", "QR", 
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
  
  ### Read in data
  Coefficients <- Coefficients.biased <- IMSE <- list()
  IMSE.avg <- matrix(NA, nrow = length(taus), ncol = length(model.names)-1)
  MSE <- BIAS <- RB <- VAR <- list()
  
  for(tau in taus){
    beta_tau = beta_calc(tau)
    
    # IMSE results from simulation
    IMSE[[match(tau, taus)]] <- read.csv(here::here("Ch3.1-5000-05162024", paste("IMSE-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
    #IMSE[[match(tau, taus)]] <- read.csv(here::here("Ch3.1-5000-05162024", "Case 6 (5% outliers)", paste("IMSE-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
    
    IMSE.avg[match(tau, taus),] <- colMeans(IMSE[[match(tau, taus)]])
    
    
    # Coefficient results from simulation
    Coefficients[[match(tau, taus)]] <- read.csv(here::here("Ch3.1-5000-05162024", paste("Coef-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
    #Coefficients[[match(tau, taus)]] <- read.csv(here::here("Ch3.1-5000-05162024", "Case 6 (5% outliers)", paste("Coef-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]

    # calculate variance, bias, MSE, robust bias
    VAR[[match(tau, taus)]] <- Variance.calc(Coefficients, tau)
    BIAS[[match(tau, taus)]] <- Bias.calc(Coefficients, tau)
    MSE[[match(tau, taus)]] <- BIAS[[match(tau, taus)]]^2 + VAR[[match(tau, taus)]]
    RB[[match(tau, taus)]] <- RB.calc(Coefficients, tau)
    
    rownames(VAR[[match(tau, taus)]]) <- rownames(BIAS[[match(tau, taus)]]) <-
      rownames(MSE[[match(tau, taus)]]) <- rownames(RB[[match(tau, taus)]]) <- model.names
    
    ### BIAS ADJUSTMENT FOR INTERCEPT
    Coefficients.biased[[match(tau, taus)]] <- Coefficients[[match(tau, taus)]] # save pre-correction coefficients
    
    adj <- c(1-2*tau, rep(0, d-1)) # adjustment for each beta value
    no.adj <- rep(0, d) # no adjustment for existing methods
    rep.adj <- c(rep(no.adj, times = 2), rep(adj, times = 30), rep(no.adj, times = 24))
    Coefficients[[match(tau, taus)]] <- sweep(Coefficients[[match(tau, taus)]], 2, rep.adj, "+")
  }
}
#################### END OF RESULTS CALCULATION ####################



#################### START OF ONLY RE PLOTS ####################
{
  png(file=paste0("Ch3.1Sim-Case", case, "REsummary.png"), width = 900, height = 600, res = 100)
  
  colnames(IMSE.avg) <- model.names[-1]
  RE.IMSE <- IMSE.avg[,1]/IMSE.avg
  
  leg.pos = "top"
  
  par(mfrow=c(3,3), mar=c(3, 3, 2, 0.5), mgp=c(2,1,0), family="Times", pty = "m") # plot settings
  
  ## ORIGINAL MRSS LOSS
  
  for(k in c(3, 5, 7)){
    taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
    if(k==3){selection=2:5}
    if(k==5){selection=6:9}
    if(k==7){
      selection=10:13
      taus = c(0.45, 0.50, 0.55)}
    
    ylim = range(RE.IMSE[,1:13], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(i) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
    if(k==3){legend(leg.pos, legend = paste0("R", 1:4), 
                    lty = 2:5, lwd = 2, col =2:5, bty = "n", seg.len = 3)}
    
  }
  
  ## APPROXIMATION MRSS LOSS
  
  for(k in c(3, 5, 7)){
    taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
    if(k==3){selection=17:20}
    if(k==5){selection=21:24}
    if(k==7){selection=25:28
    taus = c(0.45, 0.50, 0.55)}
    
    ylim = range(RE.IMSE[,c(1, 17:28)], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(ii) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
  }
  
  
  ## SRS FROM INDUCED POPULATION (CHANGING TAU)
  
  for(k in c(3, 5, 7)){
    taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
    if(k==3){selection=32:35}
    if(k==5){selection=36:39}
    if(k==7){selection=40:43}
    
    ylim = range(RE.IMSE[,c(1, 32:43)], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(iii) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
  }
  
  if(F){
  ## MEDNS IN CHECK LOSS (NOT CHANGING TAU)
  
  for(k in c(3, 5, 7)){
    taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
    if(k==3){selection=44:47}
    if(k==5){selection=48:51}
    if(k==7){selection=52:55}
    
    ylim = range(RE.IMSE[,c(1, 44:55)], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(iv) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
  }
  }
  
  
  dev.off()
}
#################### END OF ONLY RE PLOTS ####################



#################### START OF COEF TABLES ####################
tau = 0.5

avg.coefs = colMeans(Coefficients[[match(tau, taus)]])
sd.coefs<- apply(Coefficients[[match(tau, taus)]], 2, sd)

avg.coefs.mat = matrix(avg.coefs, nrow = d)
sd.coefs.mat = matrix(sd.coefs, nrow = d)
table.coefs.mat = matrix(paste(format(round(avg.coefs.mat, 2), nsmall = 2), 
                             paste0("(", format(round(sd.coefs.mat, 2), nsmall = 2), ")")), nrow = d) # combining mean and sd

col.selection = c(c(1, 2), c(3, 7, 11, 18, 22, 26, 33, 37, 41) + 3) # add 1 for R2, 2 for R3, 3 for R4
#model.names[col.selection] # check which ones are selected

(coef.table = as.table(cbind(beta, table.coefs.mat[,col.selection]))) # long version
coef.table.stacked = rbind(coef.table[,1:6], coef.table[,7:12]) # stacked version

colnames(coef.table) = c("Beta", "OLS", "QR", 
                         paste("(i)", c("k=3", "k=5", "k=7")), 
                         paste("(ii)", c("k=3", "k=5", "k=7")), 
                         paste("(iii)", c("k=3", "k=5", "k=7")))
rownames(coef.table) = c("Intercept", "X1", "X2", "X3")
rownames(coef.table.stacked) = rep(c("Intercept", "X1", "X2", "X3"), 2)

print(xtable(coef.table, type = "latex"), file = paste0("Ch3-1Sim_Case", case, "CoefTable.tex"))
print(xtable(coef.table.stacked, type = "latex"), file = paste0("Ch3-1Sim_Case", case, "CoefTable.tex"))

##################### END OF COEF TABLES #####################


##################### START OF POPULATION PLOT #####################
{
  #png(file="Ch3.1Sim-Population.png", width = 900, height = 1200, res = 100)
  
  par(mfrow=c(1, 3), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 1))
  for(c in 1:3){
    cases = c(1, 5, 6)
    case = cases[c]
    
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
      } else if(case == 6){ # 10% x and y outliers
        # Add Normal error to response
        y <- as.vector( X %*% beta + rnorm(N) )
        # Add outliers
        X[1:(0.05*N),2:4] <- 10
        y[1:(0.05*N)] <- 50
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
    
    plot(X[,2], y, xlab = "X", ylab = "Y", main = paste("Population", c, "(N=1000)"),
         pch = 16, col = "lightgrey")
    
  }
  
  #dev.off()
}
##################### END OF POPULATION PLOT #####################

