##### SIMULATION 4 PLOTS & TABLES
##### VIOLATION OF LINEAR POPULATION DATA


# load packages
library(MASS)
library(ald)
library(quantreg)
source("M-estimation Functions.R")
source("Other Functions.R")
source("RSS Functions.R")




{
  png(file="ViolationSim-Population.png", width = 900, height = 600, res = 100)
  
  #par(mfrow=c(1, 1), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 2)) # plot settings

  par(mfrow=c(1, 2), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 1))
  
  for(case in c(2, 3)){
    
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
  
  all.taus = seq(0.3, 0.7, by = 0.2)
  
  plot(x, y, xlab = "X", ylab = "Y", main = paste("Population", case, "(N=10,000)"),
       pch = 16, col = "lightgrey")
  
  
  lty.vec <- c(2, 1, 3)
  for(t in all.taus){
    curve(2 + 3*x^2  + qnorm(t, 0, 1), 
          from = -4, to = 4, add=TRUE,
          lty = lty.vec[match(t, all.taus)], lwd = 2)
  }
  
  }

  legend("bottomright", title=expression(tau), legend = all.taus, 
         lwd = 2, lty = lty.vec, seg.len = 3, bty = "n")
  
  dev.off()
}



case = 3
palette(grey.colors(10, start=0, end = 0.8))
#################### START OF RESULTS CALCULATION ####################
{
  taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
  model.names <- c("OLS", "QR", 
                   "k=3 R1", "k=3 R2", "k=3 R3", "k=3 R4",
                   "k=5 R1", "k=5 R2", "k=5 R3", "k=5 R4",
                   "k=7 R1", "k=7 R2", "k=7 R3", "k=7 R4",
                   "k=3 R1 SRS", "k=3 R2 SRS", "k=3 R3 SRS", "k=3 R4 SRS",
                   "k=5 R1 SRS", "k=5 R2 SRS", "k=5 R3 SRS", "k=5 R4 SRS",
                   "k=7 R1 SRS", "k=7 R2 SRS", "k=7 R3 SRS", "k=7 R4 SRS")
  
  ### Read in data
  IMSE <- list()
  IMSE.avg <- matrix(NA, nrow = length(taus), ncol = length(model.names))
  
  for(tau in taus){
    # IMSE results from simulation
    IMSE[[match(tau, taus)]] <- read.csv(here::here("Violation-5000-05172025", paste("Violation-ISE-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
    IMSE.avg[match(tau, taus),] <- colMeans(IMSE[[match(tau, taus)]])
  }
  colnames(IMSE.avg) <- model.names
}
#################### END OF RESULTS CALCULATION ####################




#################### START OF ONLY RE PLOTS ####################
{
  png(file=paste0("ViolationSim-Case", case, "RE.png"), width = 900, height = 600, res = 100)
  
  RE.IMSE <- IMSE.avg[,2]/IMSE.avg
  
  taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
  leg.pos = "bottomleft"
  
  par(mfrow=c(2,3), mar=c(3, 3, 2, 0.5), mgp=c(2,1,0), family="Times", pty = "m") # plot settings
  
  ## ORIGINAL MRSS LOSS
  
  for(k in c(3, 5, 7)){
    
    if(k==3){selection=3:6}
    if(k==5){selection=7:10}
    if(k==7){selection=11:14      
    taus = c(0.45, 0.50, 0.55)}
    
    ylim = range(RE.IMSE[,2:14], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(i) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
    if(k==3){legend(leg.pos, legend = paste0("R", 1:4), 
                    lty = 2:5, lwd = 2, col =2:5, bty = "n", seg.len = 3)}
    
  }
  
  
  ## SRS FROM INDUCED POPULATION (CHANGING TAU)
  
  taus = c(0.4, 0.45, 0.50, 0.55, 0.6)
  
  for(k in c(3, 5, 7)){
    if(k==3){selection=15:18}
    if(k==5){selection=19:22}
    if(k==7){selection=23:26}
    
    ylim = range(RE.IMSE[, c(2, 15:26)], na.rm = TRUE)
    matplot(taus, na.omit(RE.IMSE[,selection]), main = paste0("(iii) k=",k),
            xlab = expression(tau), ylab = "RE",
            type = "l", lty = 2:5, lwd = 2, col = 2:5,
            ylim = ylim)
    abline(h=1, col = "red")
  }
  
  dev.off()
}
#################### END OF ONLY RE PLOTS ####################

