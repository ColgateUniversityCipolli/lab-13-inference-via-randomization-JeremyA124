library(tidyverse)
library(VGAM)
set.seed(7272)

# Parameters
R <- 1000  # Number of simulations
n = 30 # Sample size

# Set up the results storage
simulation.dat <- tibble(
  sample.s.reject = rep(NA, 500),
  sample.sr.reject = rep(NA, 500)
)

for(k in 1:500) {
  sim.dat <- rlaplace(n = n, location = 0, scale = 4) #simulate rlapace data for n = 30
  sim.dat.null <- sim.dat - mean(sim.dat) #Shift the data to null
  s <- sd(sim.dat)
  resample.dat <- tibble(
    sample.s.ts = rep(NA, R),
    sample.sr.ts = rep(NA, R)
  )
  
  # Resampling loop
  for(i in 1:R) {
    resamples <- sample(sim.dat.null, size = n, replace = TRUE) # Take resamples
    xbar <- mean(resamples) 
    
    # Calculate test statistics
    resample.dat$sample.s.ts[i] <- xbar / (s / sqrt(n))  # Using original standard deviation
    resample.dat$sample.sr.ts[i] <- xbar / (sd(resamples) / sqrt(n))  # Using resampled standard deviation
  }
  
  # Observed t-statistic based on centered data
  observed.t <- mean(sim.dat) / (sd(sim.dat) / sqrt(n))  # Corrected t-statistic
  
  # Two-tailed p-value calculation
  pval.s <- mean(resample.dat$sample.s.ts >= observed.t)  # For original sample SD
  pval.sr <- mean(resample.dat$sample.sr.ts >= observed.t)  # For resampled SD
  
  # Rejection decision based on two-tailed test
  if(pval.s < 0.05) {
    simulation.dat$sample.s.reject[k] <- 1
  } else {
    simulation.dat$sample.s.reject[k] <- 0
  }
  
  if(pval.sr < 0.05) {
    simulation.dat$sample.sr.reject[k] <- 1
  } else {
    simulation.dat$sample.sr.reject[k] <- 0
  }
}

# Summarize Type I error rates
type_I_error_rate_s <- mean(simulation.dat$sample.s.reject)
type_I_error_rate_sr <- mean(simulation.dat$sample.sr.reject)

#Store capture results
simulation.mean.dat <- tibble(
  sample.capture = rep(NA, R)
)

for(i in 1:R) {
  sim.dat <- rlaplace(n = n, location = 0, scale = 4) #simulate rlapace data for n = 30
  mu0 <- mean(sim.dat) #original mean
  sim.dat.null <- sim.dat - mean(sim.dat) #Shift the data to null
  resample.dats <- tibble(
    xbars = rep(NA, R)
  )
  
  # Resampling loop
  for(k in 1:R) {
    resamples <- sample(sim.dat.null, size = n, replace = TRUE) # Take resamples
    resample.dats$xbars[k] <- mean(resamples) 
  }
  
  #Compute confidence interval 
  ci <- quantile(resample.dats$xbars, c(0.025, 0.975))
  
  #Check for coverage
  if((mu0 >= ci[1] && mu0 <= ci[2])){
    simulation.mean.dat$sample.capture[i] <- 1
  } else{
    simulation.mean.dat$sample.capture[i] <- 0
  }
}

#Calculate the capture rate
capture.porportion <- mean(simulation.mean.dat$sample.capture)
