library(tidyverse)  # For data manipulation and plotting
library(e1071)      # For skewness calculation

set.seed(7272)

# Read in the data from a CSV file containing the Zebra Finch data
finches <- read_csv("zebrafinches.csv")

# Initialize the number of data points (n) and calculate the t-statistic for the "further" data
n <- length(finches$further)
t <- t.test(finches$further, mu = 0, alternative = "less")$statistic

# Error approximation based on skewness of the data using normal distribution
error <- pnorm(t) + (skewness(finches$further) / sqrt(n)) * ((2 * (t ** 2) + 1) / 6) * dnorm(t)

# Create a tibble to store t-values and the associated error for visualization
error.stats <- tibble(
  t = rep(NA, 201),      # Placeholder for t values
  error = rep(NA, 201)   # Placeholder for calculated error values
)

# Loop over a range of t-values to compute error and store the results
k = 1
for(i in seq(-10, 10, 0.1)){
  error <- (skewness(finches$further) / sqrt(n)) * ((2 * (i ** 2) + 1) / 6) * dnorm(i)
  error.stats$t[k] <- i
  error.stats$error[k] <- error
  k = k + 1
}

# Plot the error values against t-values using ggplot2
ggplot() +
  geom_line(data = error.stats, aes(x = t, y = error)) +
  theme_bw()  # Using black and white theme for clarity

# Use quantile to get a specific value of t at the 5th percentile of the "further" data
t <- quantile(finches$further, 0.05)

# Calculate the sample size n based on skewness, standard deviation, and a chosen error margin (0.1) at the 5th percentile
n <- ((skewness(finches$further) / (6 * 0.10 * 0.05)) * (2 * (t ** 2) + 1) * dnorm(t)) ** 2

# Define the number of bootstrap iterations (R)
R <- 1000

# Function for bootstrapping and shifting resamples to approximate the sampling distribution
bootstrap.shift <- function(data, s, n, R) {
  data <- data - mean(data)
  resampled.data <- tibble(ts = rep(NA, R))  # Initialize a tibble to store the bootstrap t-statistics
  
  # Perform resampling R times
  for(i in 1:R) {
    resamples <- sample(data, size = n, replace = TRUE)  # Resample with replacement
    resampled.data$ts[i] <- mean(resamples) / (s / sqrt(n))  # Calculate t-statistic for each resample
  }
  
  resampled.shifted.data <- tibble(ts = resampled.data$ts) #store data in tibble
  
  return(resampled.shifted.data)  # Return the shifted resample data
}

# Generate shifted resamples for the 'closer', 'further', and 'difference' data
resamples.null.closer <- bootstrap.shift(finches$closer, s = sd(finches$closer), n, R)
resamples.null.farther <- bootstrap.shift(finches$further, s = sd(finches$further), n, R)
resamples.null.diff <- bootstrap.shift(finches$diff, s = sd(finches$diff), n, R)

# Function for calculating bootstrap p-values
bootstrap.pval <- function(shifted.data, observed.t, method) {
  if(method == "less") {
    boot.pval <- mean(shifted.data$ts <= observed.t)  # P-value for a "less" alternative
  } else if(method == "greater") {
    boot.pval <- mean(shifted.data$ts >= observed.t)  # P-value for a "greater" alternative
  } else if(method == "two.sided") {
    boot.pval <- mean(abs(shifted.data$ts) >= abs(observed.t))  # Two-sided test p-value
  } else {
    stop("Enter a valid method.")  # Handle invalid method input
  }
  
  return(boot.pval)  # Return the computed p-value
}

# Compute bootstrap p-values for each group: 'closer', 'further', and 'difference'
closer.pval.boot <- bootstrap.pval(resamples.null.closer, observed.t = mean(finches$closer) / (sd(finches$closer) / sqrt(n)), method = "greater")
farther.pval.boot <- bootstrap.pval(resamples.null.farther, observed.t = mean(finches$further) / (sd(finches$further) / sqrt(n)), method = "less")
diff.pval.boot <- bootstrap.pval(resamples.null.closer, observed.t = mean(finches$diff) / (sd(finches$diff) / sqrt(n)), method = "two.sided")

# Compute t-test p-values for comparison
closer.pval.ttest <- t.test(finches$closer, mu = 0, alternative = "greater")$p.value
farther.pval.ttest <- t.test(finches$further, mu = 0, alternative = "less")$p.value
diff.pval.ttest <- t.test(finches$diff, mu = 0, alternative = "two.sided")$p.value

# Compute the 5th percentile from the shifted resamples (bootstrap) and compare to the theoretical 5th percentile t-value
closer.5th <- quantile(resamples.null.closer$ts, probs = 0.05)
farther.5th <- quantile(resamples.null.farther$ts, probs = 0.05)
diff.5th <- quantile(resamples.null.diff$ts, probs = 0.05)
t.5th <- qt(0.05, df = n - 1)  # Theoretical t-value at the 5th percentile

# Function to compute bootstrap confidence intervals
bootstrap.CI <- function(data, s, n, R) {
  resampled.data <- tibble(ts = rep(NA, R))  # Initialize tibble to store means
  
  # Perform resampling R times
  for(i in 1:R) {
    resamples <- sample(data, size = n, replace = TRUE)  # Resample with replacement
    resampled.data$xbars[i] <- mean(resamples)  # Store the mean of each resample
  }
  
  return(quantile(resampled.data$xbars, c(0.025, 0.975)))  # Return 95% confidence interval using percentiles
}

# Compute bootstrap confidence intervals for the three datasets
boot.CI.closer <- bootstrap.CI(finches$closer, s = sd(finches$closer), n, R)
boot.CI.farther <- bootstrap.CI(finches$further, s = sd(finches$further), n, R)
boot.CI.diff <- bootstrap.CI(finches$diff, s = sd(finches$diff), n, R)

# Compute t-test confidence intervals for comparison
ttest.CI.closer <- t.test(finches$closer, mu = 0, alternative = "two.sided")$conf.int
ttest.CI.farther <- t.test(finches$further, mu = 0, alternative = "two.sided")$conf.int
ttest.CI.diff <- t.test(finches$diff, mu = 0, alternative = "two.sided")$conf.int

# Randomization procedure to compute p-values and confidence intervals for a random hypothesis test
random.pval <- function(data, mu0, method = "two.sided", R = 1000) {
  rand <- tibble(means = rep(NA, R))  # Initialize tibble to store means
  
  x.shift <- data - mu0  # Shift the data by subtracting the hypothesized mean
  
  for(i in 1:R) {
    curr.rand <- x.shift * sample(x = c(-1, 1), size = length(x.shift), replace = TRUE)  # Randomize the signs
    rand$means[i] <- mean(curr.rand)  # Store the mean of each randomization
  }
  
  delta <- abs(mean(data))  # Calculate the difference from the hypothesized mean
  low <- 0 - delta  # Lower bound for the CI
  high <- 0 + delta  # Upper bound for the CI
  
  # Compute p-values based on the randomization
  if(method == "less") {
    return(mean(rand$means <= low))
  } else if(method == "greater") {
    return(mean(rand$means >= high))
  } else if(method == "two.sided") {
    return(mean(rand$means <= low) + mean(rand$means >= high))
  } else {
    stop("Enter a valid method.")  # Handle invalid method input
  }
}

# Compute random p-values for each group
closer.rand.pval <- random.pval(finches$closer, 0, "greater")
farther.rand.pval <- random.pval(finches$further, 0, "less")
diff.rand.pval <- random.pval(finches$diff, 0)

# Randomized Confidence Interval procedure
random.CI <- function(data,
                      R = 1000){
  mu0.iterate <- 0.01 #Iteration Variable
  starting.point <- mean(data) #Starting Point
  mu.lower <- starting.point # Initialize lower bound for CI
  
  repeat{
    rand <- tibble(means = rep(NA, R)) # Initialize tibble for storing means
    x.shift <- data - mu.lower # Shift data
    
    for(i in 1:R){
      curr.rand <- x.shift * # Randomize the data
        sample(x = c(-1,1),
               size = length(x.shift),
               replace = TRUE)
      
      rand$means[i] <- mean(curr.rand) # Store mean of each randomization
    }
    
    delta <- abs(mean(data)) # Calculate delta from the me
    low <- 0 - delta # Lower bound for CI
    high<- 0 + delta # Upper bound for CI
    
    rand <- rand |>
      mutate(means = means + mu.lower) # Compute p-value
    delta <- abs(mean(data) - mu.lower)
    low <- mu.lower - delta 
    high<- mu.lower + delta   
    p.val <- mean(rand$means <= low) +
      mean(rand$means >= high)
    
    if(p.val < 0.05){ # Stop if p-value is less than 0.05
      break
    }else{
      mu.lower <- mu.lower - mu0.iterate # Otherwise, shift lower bound further
    }
  }
  
  mu0.iterate <- 0.01  # Repeat for upper bound
  starting.point <- mean(data)
  mu.higher <- starting.point
  
  repeat{
    rand <- tibble(means = rep(NA, R))
    x.shift <- data - mu.higher
    
    for(i in 1:R){
      curr.rand <- x.shift *
        sample(x = c(-1,1),
               size = length(x.shift),
               replace = TRUE)
      
      rand$means[i] <- mean(curr.rand)
    }
    
    delta <- abs(mean(data))
    high <- 0 + delta
    
    rand <- rand |>
      mutate(means = means + mu.higher)
    delta <- abs(mean(data) - mu.higher)
    low <- mu.higher - delta 
    high<- mu.higher + delta   
    p.val <- mean(rand$means <= low) +
      mean(rand$means >= high)
    
    if(p.val < 0.05){
      break
    }else{
      mu.higher <- mu.higher + mu0.iterate # Shift upper bound higher
    }
  }
  
  return(c(mu.lower, mu.higher))
}

# Compute random confidence intervals for each group
closer.rand.CI <- random.CI(finches$closer)
farther.rand.CI <- random.CI(finches$further)
diff.rand.CI <- random.CI(finches$diff)
