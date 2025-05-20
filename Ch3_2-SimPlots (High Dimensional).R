##### CHAPTER 3: SIMULATION 2 PLOTS & TABLES
##### HIGH DIMENSIONAL DATA

# load packages
library(MASS)
library(ald)
library(quantreg)
library(xtable)
source("/Users/neveloewen/Dropbox/Neve - Research Materials/Research - R Code/Functions/M-estimation Functions.R")
source("/Users/neveloewen/Dropbox/Neve - Research Materials/Research - R Code/Functions/Other Functions.R")
source("/Users/neveloewen/Dropbox/Neve - Research Materials/Research - R Code/Functions/RSS Functions.R")

################################################## 
#################### SETTINGS #################### 
##################################################

case = 3
tau = 0.5 # same for all cases
model.names <- c("MeanReg", "QR", 
                 "k=3 R1", "k=3 R2", "k=3 R3", "k=3 R4",
                 "k=5 R1", "k=5 R2", "k=5 R3", "k=5 R4",
                 "k=7 R1", "k=7 R2", "k=7 R3", "k=7 R4",
                 "k=3 R1 SRS", "k=3 R2 SRS", "k=3 R3 SRS", "k=3 R4 SRS",
                 "k=5 R1 SRS", "k=5 R2 SRS", "k=5 R3 SRS", "k=5 R4 SRS",
                 "k=7 R1 SRS", "k=7 R2 SRS", "k=7 R3 SRS", "k=7 R4 SRS")

if(case == 1){
  beta = c(0, rgamma(10, shape = 1, rate = 10), rep(1, 5), 
           rgamma(10, shape = 1, rate = 10), rep(4, 5))
}else if(case == 2){
  beta = c(0, rep(3, 5), 
           rgamma(10, shape = 1, rate = 10), rep(1, 10), 
           rgamma(10, shape = 1, rate = 10), rep(2, 5))
}else if(case == 3){
  beta = c(0, rgamma(10, shape = 1, rate = 10), rep(5, 3), 
           rgamma(10, shape = 1, rate = 10), rep(1, 5),
           rgamma(10, shape = 1, rate = 10), rep(10, 2))
}

# this is the same for all cases
d = length(beta)
gamma = c(1, rep(0, d-1))
lambdas = exp(seq(0, 5, length.out = 6))



################################################## 
###################### DATA ######################
##################################################
IMSE <- read.csv(here::here("Ch3.2-1000-06302024", paste("1000HighDim-IMSE-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
Coefficients <- read.csv(here::here("Ch3.2-1000-06302024", paste("1000HighDim-Coef-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
PATH <- read.csv(here::here("Ch3.2-1000-06302024", paste("1000HighDim-PATH-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]


################################################## 
#################### IMSE / RE ################### 
##################################################

IMSE.avg <- colMeans(IMSE)
RE.IMSE <- IMSE.avg[1]/IMSE.avg
RE.table <- matrix(RE.IMSE[-1], ncol = 4, byrow = TRUE)
colnames(RE.table) <- c("R1", "R2", "R3", "R4")
rownames(RE.table) <- c("k=3", "k=5", "k=7", "k=3 SRS", "k=5 SRS", "k=7 SRS")
print(xtable(RE.table, type = "latex"), file = paste0("Ch3-2Sim_Case", case, "IMSETable.tex"))


################################################## 
################## COEFFICIENTS ################## 
##################################################

coef.avg <- colMeans(Coefficients)
coef.avg.mat <- matrix(coef.avg, ncol = d, byrow = TRUE)
rownames(coef.avg.mat) = model.names
#round(coef.avg.mat, 3)


################################################## 
################### PATH PLOTS ################### 
##################################################

L = length(lambdas)

avg.path.meanLASSO = PATH[,1:L]
avg.path.QRLASSO = PATH[,(L+1):(2*L)]
avg.path.K3LASSO = PATH[,(2*L+1):(6*L)]
avg.path.K5LASSO = PATH[,(6*L+1):(10*L)]
avg.path.K7LASSO = PATH[,(10*L+1):(14*L)]
avg.path.K3SRSinducedLASSO = PATH[,(14*L+1):(18*L)]
avg.path.K5SRSinducedLASSO = PATH[,(18*L+1):(22*L)]
avg.path.K7SRSinducedLASSO = PATH[,(22*L+1):(26*L)]

{
  png(file=paste0("Ch3.2Sim-Case", case, "PathMedNS.png"), width = 900, height = 600, res = 100)
  
  par(mfrow = c(2,2), family="Times", pty = "m", mar=c(4, 4, 2, 1), mgp=c(2.5,1,0)) # plot settings
  
  matplot(t(avg.path.QRLASSO[-1,]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "QR")
  matplot(t(avg.path.K3LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=3") # R2
  matplot(t(avg.path.K5LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=5") # R2
  matplot(t(avg.path.K7LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=7") # R2
  
  dev.off()
}

# SRS INDUCED
{
  png(file=paste0("Ch3.2Sim-Case", case, "PathSRSinduced.png"), width = 900, height = 600, res = 100)
  
  par(mfrow = c(2,2), family="Times", pty = "m", mar=c(4, 4, 2, 1), mgp=c(2.5,1,0)) # plot settings
  
  matplot(t(avg.path.QRLASSO[-1,]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "QR")
  matplot(t(avg.path.K3SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=3")
  matplot(t(avg.path.K5SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=5")
  matplot(t(avg.path.K7SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "k=7")
  
  dev.off()
}


# ALL IN ONE PLOT
{
  png(file=paste0("Ch3.2Sim-Case", case, "Path.png"), width = 900, height = 1200, res = 100)
  
  par(mfrow = c(4,2), family="Times", pty = "m", mar=c(4, 4, 2, 1), mgp=c(2.5,1,0)) # plot settings
  
  matplot(t(avg.path.meanLASSO[-1,]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "OLS")
  matplot(t(avg.path.QRLASSO[-1,]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "QR")
  matplot(t(avg.path.K3LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (i) k=3") # R2
  matplot(t(avg.path.K3SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (iii) k=3")
  matplot(t(avg.path.K5LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (i) k=5") # R2
  matplot(t(avg.path.K5SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (iii) k=5")
  matplot(t(avg.path.K7LASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (i) k=7") # R2
  matplot(t(avg.path.K7SRSinducedLASSO[-1,(L+1):(2*L)]), type = "l",
          xlab = expression("log("*lambda*")"), ylab = "Coefficients",
          main = "Method (iii) k=7")
  
  dev.off()
}






