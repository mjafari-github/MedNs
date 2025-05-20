##### REAL DATA: CI RESULTS #####
##### March 12, 2024

library(plotrix)

#################### START OF POPULATION PLOTS ####################

# loading and cleaning population
{
  # load data
  BodyFat <- read.csv(here::here("Data", "Body Fat Data", "BodyFat - Extended.csv"))
  
  # clean data
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

cor(R1, BodyFat$BodyFat)
cor(R1[-42], BodyFat$BodyFat[-42])

{
  # histogram of response
  
  png(file="BodyFatHist.png", width = 900, height = 400, res = 100)
  par(mfrow=c(1, 2), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 1)) 
  
  hist(BodyFat$BodyFat[which(BodyFat$Sex == 1)], xlab = "Body Fat (%)", main = "Male")# , col = 4)
  hist(BodyFat$BodyFat[which(BodyFat$Sex == 0)], xlab = "Body Fat (%)", main = "Female")#, col = 2)
  
  dev.off()
}

{
  png(file="BodyFatScatter.png", width = 900, height = 400, res = 100)
  par(mfrow = c(1, 1), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 2)) # plot settings
  
  # plot of ranker against response
  plot(y, R1, pch = 16, col = "dark grey",
       ylab = "Ranker", xlab = "Body Fat (%)")
  #plot(y[-42], R1[-42], pch = 16, col = "dark grey",
  #     ylab = "Ranker", xlab = "Body Fat (%)")
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  dev.off()
}

{
  # COOKS DISTANCE
  model = lm(y ~ features)
  cooksd = cooks.distance(model)
  
  cols = ifelse(cooksd>4*mean(cooksd, na.rm=T), 2, 1)
  
  png(file="BodyFatCooks.png", width = 900, height = 400, res = 100)
  par(mfrow = c(1, 1), family="Times", pty = "m", mgp=c(2.5,1,0), mar=c(4, 4, 2, 2)) # plot settings
  
  
  plot(cooksd, pch=4, main=" ", ylab = "Cook's Distance", col = cols)  # plot cook's distance
  abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
 
  dev.off()
  
  # outliers according to cooks distance
  BodyFat[cooksd>4*mean(cooksd, na.rm=T),]
}
{
  # MAHALANOBIS
  data <- BodyFat[,6:15]
  distances <- mahalanobis(x = data , center = colMeans(data) , cov = cov(data))
  # Cutoff value for ditances from Chi-Sqaure Dist. 
  # with p = 0.95 df = 2 which in ncol(air)
  cutoff <- qchisq(p = 0.95 , df = ncol(data))
  ## Display observation whose distance greater than cutoff value
  data[distances > cutoff ,]
  dim(BodyFat[distances > cutoff ,])
}

#################### END OF POPULATION PLOTS ####################










#################### START OF CI PLOTS ####################
{
  
ns = c(40, 60, 80, 100) # sample size 
R = 1 # number of rankers
method = 3 # choose 1: exact loss, 2: approx loss, 3: srs induced

png(file=paste0("BodyFatCI-Method", method, ".png"), width = 900, height = 600, res = 100)  
par(mfrow = c(2,2), family="Times", pty = "m", mar=c(4, 4, 2, 1), mgp=c(2.5,1,0)) # plot settings

for(n in ns){

coef.all <- read.csv(here::here("Data", "Body Fat Data", paste0("CoefResults-BodyFat", n, ".csv")))[,-1] 

# using complete cases coefficients
if(T){
  
  coef.complete <- coef.all[complete.cases(coef.all),]
  
  # Coefficients for each method
  coef.mean <- coef.complete[,1:d]
  coef.QR <- coef.complete[,(d+1):(2*d)]
  coef.K3 <- coef.complete[,(2*d+1):((R+2)*d)]
  coef.K5 <- coef.complete[,((R+2)*d+1):((2*R+2)*d)]
  coef.K7 <- coef.complete[,((2*R+2)*d+1):((3*R+2)*d)]
  
  coef.K3approx <- coef.complete[,((3*R+2)*d+1):((3*R+3)*d)]
  coef.K5approx <- coef.complete[,((3*R+3)*d+1):((3*R+4)*d)]
  coef.K7approx <- coef.complete[,((3*R+4)*d+1):((3*R+5)*d)]

  coef.K3SRSinduced <- coef.complete[,((3*R+5)*d+1):((4*R+5)*d)]
  coef.K5SRSinduced <- coef.complete[,((4*R+5)*d+1):((5*R+5)*d)]
  coef.K7SRSinduced <- coef.complete[,((5*R+5)*d+1):((6*R+5)*d)]
  
  coef.K3SRS <- coef.complete[,((6*R+5)*d+1):((6*R+6)*d)]
  coef.K5SRS <- coef.complete[,((6*R+6)*d+1):((6*R+7)*d)]
  coef.K7SRS <- coef.complete[,((6*R+7)*d+1):((6*R+8)*d)]
  
  coef.K3SRSapprox <- coef.complete[,((6*R+8)*d+1):((7*R+8)*d)]
  coef.K5SRSapprox <- coef.complete[,((7*R+8)*d+1):((8*R+8)*d)]
  coef.K7SRSapprox <- coef.complete[,((8*R+8)*d+1):((9*R+8)*d)]
  
  # Average coefficients
  avg.coef.mean <- apply(coef.mean, 2, mean)
  avg.coef.QR <- apply(coef.QR, 2, mean)
  avg.coef.K3 <- matrix(apply(coef.K3, 2, mean), ncol = R)
  avg.coef.K5 <- matrix(apply(coef.K5, 2, mean), ncol = R)
  avg.coef.K7 <- matrix(apply(coef.K7, 2, mean), ncol = R)
  avg.coef.K3SRS <- apply(coef.K3SRS, 2, mean)
  avg.coef.K5SRS <- apply(coef.K5SRS, 2, mean)
  avg.coef.K7SRS <- apply(coef.K7SRS, 2, mean)
  avg.coef.K3SRSinduced <- matrix(apply(coef.K3SRSinduced, 2, mean), ncol = R)
  avg.coef.K5SRSinduced <- matrix(apply(coef.K5SRSinduced, 2, mean), ncol = R)
  avg.coef.K7SRSinduced <- matrix(apply(coef.K7SRSinduced, 2, mean), ncol = R)
  avg.coef.K3approx <- matrix(apply(coef.K3approx, 2, mean), ncol = R)
  avg.coef.K5approx <- matrix(apply(coef.K5approx, 2, mean), ncol = R)
  avg.coef.K7approx <- matrix(apply(coef.K7approx, 2, mean), ncol = R)
  avg.coef.K3SRSapprox <- apply(coef.K3SRSapprox, 2, mean)
  avg.coef.K5SRSapprox <- apply(coef.K5SRSapprox, 2, mean)
  avg.coef.K7SRSapprox <- apply(coef.K7SRSapprox, 2, mean)
  
  # variance 
  sd.coef.mean <- apply(coef.mean, 2, sd)
  sd.coef.QR <- apply(coef.QR, 2, sd)
  sd.coef.K3 <- matrix(apply(coef.K3, 2, sd), ncol = R)
  sd.coef.K5 <- matrix(apply(coef.K5, 2, sd), ncol = R)
  sd.coef.K7 <- matrix(apply(coef.K7, 2, sd), ncol = R)
  sd.coef.K3SRS <- apply(coef.K3SRS, 2, sd)
  sd.coef.K5SRS <- apply(coef.K5SRS, 2, sd)
  sd.coef.K7SRS <- apply(coef.K7SRS, 2, sd)
  sd.coef.K3SRSinduced <- matrix(apply(coef.K3SRSinduced, 2, sd), ncol = R)
  sd.coef.K5SRSinduced <- matrix(apply(coef.K5SRSinduced, 2, sd), ncol = R)
  sd.coef.K7SRSinduced <- matrix(apply(coef.K7SRSinduced, 2, sd), ncol = R)
  sd.coef.K3approx <- matrix(apply(coef.K3approx, 2, sd), ncol = R)
  sd.coef.K5approx <- matrix(apply(coef.K5approx, 2, sd), ncol = R)
  sd.coef.K7approx <- matrix(apply(coef.K7approx, 2, sd), ncol = R)
  sd.coef.K3SRSapprox <- apply(coef.K3SRSapprox, 2, sd)
  sd.coef.K5SRSapprox <- apply(coef.K5SRSapprox, 2, sd)
  sd.coef.K7SRSapprox <- apply(coef.K7SRSapprox, 2, sd)
}

# confidence bounds: normal
if(T){
  
  ### CONFIDENCE INTERVAL BOUNDS
  z = qnorm(0.975)
  
  # Lower bounds
  lb.coef.mean <- avg.coef.mean - z*sd.coef.mean  
  lb.coef.QR <- avg.coef.QR - z*sd.coef.QR   
  lb.coef.K3 <- avg.coef.K3 - z*sd.coef.K3   
  lb.coef.K5 <- avg.coef.K5 - z*sd.coef.K5   
  lb.coef.K7 <- avg.coef.K7 - z*sd.coef.K7   
  lb.coef.K3SRS <- avg.coef.K3SRS - z*sd.coef.K3SRS   
  lb.coef.K5SRS <- avg.coef.K5SRS - z*sd.coef.K5SRS   
  lb.coef.K7SRS <- avg.coef.K7SRS - z*sd.coef.K7SRS   
  lb.coef.K3SRSinduced <- avg.coef.K3SRSinduced - z*sd.coef.K3SRSinduced   
  lb.coef.K5SRSinduced <- avg.coef.K5SRSinduced - z*sd.coef.K5SRSinduced   
  lb.coef.K7SRSinduced <- avg.coef.K7SRSinduced - z*sd.coef.K7SRSinduced   
  lb.coef.K3approx <- avg.coef.K3approx - z*sd.coef.K3approx   
  lb.coef.K5approx <- avg.coef.K5approx - z*sd.coef.K5approx   
  lb.coef.K7approx <- avg.coef.K7approx - z*sd.coef.K7approx   
  lb.coef.K3SRSapprox <- avg.coef.K3SRSapprox - z*sd.coef.K3SRSapprox   
  lb.coef.K5SRSapprox <- avg.coef.K5SRSapprox - z*sd.coef.K5SRSapprox   
  lb.coef.K7SRSapprox <- avg.coef.K7SRSapprox - z*sd.coef.K7SRSapprox   
  
  # Upper bounds
  ub.coef.mean <- avg.coef.mean + z*sd.coef.mean   
  ub.coef.QR <- avg.coef.QR + z*sd.coef.QR   
  ub.coef.K3 <- avg.coef.K3 + z*sd.coef.K3   
  ub.coef.K5 <- avg.coef.K5 + z*sd.coef.K5   
  ub.coef.K7 <- avg.coef.K7 + z*sd.coef.K7   
  ub.coef.K3SRS <- avg.coef.K3SRS + z*sd.coef.K3SRS   
  ub.coef.K5SRS <- avg.coef.K5SRS + z*sd.coef.K5SRS   
  ub.coef.K7SRS <- avg.coef.K7SRS + z*sd.coef.K7SRS   
  ub.coef.K3SRSinduced <- avg.coef.K3SRSinduced + z*sd.coef.K3SRSinduced   
  ub.coef.K5SRSinduced <- avg.coef.K5SRSinduced + z*sd.coef.K5SRSinduced   
  ub.coef.K7SRSinduced <- avg.coef.K7SRSinduced + z*sd.coef.K7SRSinduced   
  ub.coef.K3approx <- avg.coef.K3approx + z*sd.coef.K3approx   
  ub.coef.K5approx <- avg.coef.K5approx + z*sd.coef.K5approx   
  ub.coef.K7approx <- avg.coef.K7approx + z*sd.coef.K7approx   
  ub.coef.K3SRSapprox <- avg.coef.K3SRSapprox + z*sd.coef.K3SRSapprox   
  ub.coef.K5SRSapprox <- avg.coef.K5SRSapprox + z*sd.coef.K5SRSapprox   
  ub.coef.K7SRSapprox <- avg.coef.K7SRSapprox + z*sd.coef.K7SRSapprox   
  
}

r = 1 # which ranker

if(method == 1){
## EXACT LOSS FUNCTION
lower.bounds = as.vector(matrix(c(lb.coef.QR, lb.coef.K3[,r], lb.coef.K5[,r], lb.coef.K7[,r]), 
                                ncol = d, byrow = TRUE))
upper.bounds = as.vector(matrix(c(ub.coef.QR, ub.coef.K3[,r], ub.coef.K5[,r], ub.coef.K7[,r]), 
                                ncol = d, byrow = TRUE))
means = as.vector(matrix(c(avg.coef.QR, avg.coef.K3[,r], avg.coef.K5[,r], avg.coef.K7[,r]), 
                         ncol = d, byrow = TRUE))
}
if(method == 2){
# APPROX
lower.bounds = as.vector(matrix(c(lb.coef.QR, lb.coef.K3approx[,r], lb.coef.K5approx[,r], lb.coef.K7approx[,r]), 
                                ncol = d, byrow = TRUE))
upper.bounds = as.vector(matrix(c(ub.coef.QR, ub.coef.K3approx[,r], ub.coef.K5approx[,r], ub.coef.K7approx[,r]), 
                                ncol = d, byrow = TRUE))
means = as.vector(matrix(c(avg.coef.QR, avg.coef.K3approx[,r], avg.coef.K5approx[,r], avg.coef.K7approx[,r]), 
                         ncol = d, byrow = TRUE))
}
if(method == 3){
# SRS INDUCES
lower.bounds = as.vector(matrix(c(lb.coef.QR, lb.coef.K3SRSinduced[,r], lb.coef.K5SRSinduced[,r], lb.coef.K7SRSinduced[,r]), 
                                ncol = d, byrow = TRUE))
upper.bounds = as.vector(matrix(c(ub.coef.QR, ub.coef.K3SRSinduced[,r], ub.coef.K5SRSinduced[,r], ub.coef.K7SRSinduced[,r]), 
                                ncol = d, byrow = TRUE))
means = as.vector(matrix(c(avg.coef.QR, avg.coef.K3SRSinduced[,r], avg.coef.K5SRSinduced[,r], avg.coef.K7SRSinduced[,r]), 
                         ncol = d, byrow = TRUE))
}

# indicates CIs not containing 0
colours = ifelse(sign(upper.bounds) == sign(lower.bounds), 2, 1)
# beta labels on x axis
x.labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]),
             expression(beta[4]), expression(beta[5]), expression(beta[6]),
             expression(beta[7]), expression(beta[8]), expression(beta[9]), expression(beta[10]))

# exclude intercept from plots
plot.selection = 5:(4*d)


plotCI(y = means[plot.selection], x = 1:length(plot.selection), #x = 1:(4*d), 
       li = lower.bounds[plot.selection], ui = upper.bounds[plot.selection], 
       err = "y", 
       col = colours[plot.selection], lwd = colours[plot.selection],
       main = paste0("n = ", n), ylab = "Mean Estimate", xlab = "",
       xaxt = "n")
axis(1, at = seq(from = 2.5, to = 40.5, by = 4),
     labels = x.labels)
abline(v = seq(4, 4*(d-2), by = 4)+0.5, col = "lightgrey") # parameter dividing lines
abline(h = 0, col = 2, lwd=2) # line at mean = 0

}

dev.off()

}

