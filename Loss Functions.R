# Median Ranked Set Sampling (MRSS) Loss

# When set size is odd:

med_rho1 <- function(x, m, tau){
  (tau - 1)*x - ((m-1)/(2*m-1))*log((exp((tau-1)*x)-tau)/(1-tau))
}

med_rho2 <- function(x, m, tau){
  tau*x - ((m-1)/(2*m-1))*log((exp(tau*x)-(1-tau))/tau)
}

MRSSloss <- function(x, m, tau){
  ifelse(x <= 0, med_rho1(x, m, tau), med_rho2(x, m, tau))
}


# regret loss function
MedNSRegret <- function(x, m, tau){
  c1 = log((2*m-1)*tau/m) - ((m-1)/(2*m-1))*log(((m-1)/m)*(tau/(1-tau)))
  c2 = log((2*m-1)*(1-tau)/m) - ((m-1)/(2*m-1))*log(((m-1)/m)*((1-tau)/tau))
  ifelse(tau <= 0.5, MRSSloss(x, m, tau) - c1, MRSSloss(x, m, tau) - c2)
}

# MRSS Loss Approximation
MRSSlossApprox <- function(x, m, tau){
  ifelse(x<= 0, (tau-1)*x, tau*x) - ((m-1)/(2*m-1))*abs(x) + ((m-1)/(2*(2*m-1)))*x^2
}
