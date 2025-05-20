##### CHAPTER 3: SIMULATION 3 PLOTS & TABLES
##### NON-LINEAR DATA


# load packages
library(MASS)
library(ald)
library(quantreg)
#source("M-estimation Functions.R")
#source("Other Functions.R")
#source("RSS Functions.R")




{
  png(file="Ch3_3Sim-Population.png", width = 900, height = 800, res = 100)
  
  par(mfrow=c(2, 1), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 1))
  for(case in 1:2){
    
    if(T){
      set.seed(8)
      N = 10000 # population size
      
      if(case == 1){
        x <- rnorm(N, 0, 1) # covariate
        Z <- rnorm(N, 0, 1) # error
        y <- 2.5 + sin(2*x) + 2*exp(-16*x^2) + 0.5*Z # response
      } else if(case == 2){
        x <- rnorm(N, 0, 1) # covariate
        E <- rexp(N, 1) # error
        y <- 2 + 2*cos(x) + exp(-4*x^2) + E # response
      }
    }
    
    
    all.taus = seq(0.3, 0.7, by = 0.1)
    
    plot(x, y, xlab = "X", ylab = "Y", main = paste("Population", case, "(N=10,000)"),
         pch = 16, col = "lightgrey")
    
    
    lty.vec <- c(2, 3, 1, 4, 5)
    for(t in all.taus){
      
      if(case==1){
      curve(2.5 + sin(2*x) + 2*exp(-16*x^2) + qnorm(t, 0, 1), 
            from = -4, to = 4, add=TRUE,
            lty = lty.vec[match(t, all.taus)], lwd = 2)
      }
      if(case == 2){
        curve(2 + 2*cos(x) + exp(-4*x^2) + qexp(t, 1), 
              from = -4, to = 4, add=TRUE,
              lty = lty.vec[match(t, all.taus)], lwd = 2)
      }
    }
  }
  
  legend("topright", title=expression(tau), legend = all.taus, 
         lwd = 2, lty = lty.vec, seg.len = 3, bty = "n")
  
  dev.off()
}



case = 1
palette(grey.colors(10, start=0, end = 0.8))
#################### START OF RESULTS CALCULATION ####################
{
  taus = c(0.45, 0.50, 0.55)
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
    IMSE[[match(tau, taus)]] <- read.csv(here::here("Ch3.3-5000-07082024", paste("ISE-Case", case, "-tau", tau, ".csv", sep = "")))[,-1]
    IMSE.avg[match(tau, taus),] <- colMeans(IMSE[[match(tau, taus)]])
  }
  colnames(IMSE.avg) <- model.names
}
#################### END OF RESULTS CALCULATION ####################




#################### START OF ONLY RE PLOTS ####################
{
  png(file=paste0("Ch3.3Sim-Case", case, "RE.png"), width = 900, height = 600, res = 100)
  
  RE.IMSE <- IMSE.avg[,2]/IMSE.avg
  
  taus = c(0.45, 0.50, 0.55)
  leg.pos = "top"
  
  par(mfrow=c(2,3), mar=c(3, 3, 2, 0.5), mgp=c(2,1,0), family="Times", pty = "m") # plot settings
  
  ## ORIGINAL MRSS LOSS
  
  for(k in c(3, 5, 7)){

    if(k==3){selection=3:6}
    if(k==5){selection=7:10}
    if(k==7){selection=11:14}
    
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


# Presentation population plot

{
  png(file="Presentation-Ch3.3Sim-Population.png", width = 450, height = 600, res = 100)  
  
  par(mfrow = c(2,1), family="Times", pty = "m", mar=c(4, 4, 2, 1), mgp=c(2.5,1,0)) # plot settings
  for(case in 1:2){
    
    if(T){
      set.seed(8)
      N = 10000 # population size
      
      if(case == 1){
        x <- rnorm(N, 0, 1) # covariate
        Z <- rnorm(N, 0, 1) # error
        y <- 2.5 + sin(2*x) + 2*exp(-16*x^2) + 0.5*Z # response
      } else if(case == 2){
        x <- rnorm(N, 0, 1) # covariate
        E <- rexp(N, 1) # error
        y <- 2 + 2*cos(x) + exp(-4*x^2) + E # response
      }
    }
    
    
    all.taus = seq(0.3, 0.7, by = 0.1)
    
    plot(x, y, xlab = "X", ylab = "Y", main = paste("Population", case, "(N=10,000)"),
         pch = 16, col = "lightgrey")
    
    
    lty.vec <- c(2, 3, 1, 4, 5)
    for(t in all.taus){
      
      if(case==1){
        curve(2.5 + sin(2*x) + 2*exp(-16*x^2) + qnorm(t, 0, 1), 
              from = -4, to = 4, add=TRUE,
              lty = lty.vec[match(t, all.taus)], lwd = 2)
      }
      if(case == 2){
        curve(2 + 2*cos(x) + exp(-4*x^2) + qexp(t, 1), 
              from = -4, to = 4, add=TRUE,
              lty = lty.vec[match(t, all.taus)], lwd = 2)
      }
    }
  }
  
  legend("topright", title=expression(tau), legend = all.taus, 
         lwd = 2, lty = lty.vec, seg.len = 3, bty = "n")
  
  dev.off()
}
