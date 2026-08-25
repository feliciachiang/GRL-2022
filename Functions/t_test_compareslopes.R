# Function to conduct two-tailed t-test to test whether slopes of two datasets are equal or not
# Assumes error variances are approximately equal, pools estimates of the error variances, weighing each by their degrees of freedom
# 
# Website reference: https://www.real-statistics.com/regression/hypothesis-testing-significance-regression-line-slope/comparing-slopes-two-independent-samples/
# Academic reference: Lothar Sachs (1982), "Applied Statistics: A Handbook of Techniques", Springer-Verlag, pp. 413-442.
# Input: 2 2-dimensional datasets to use
# Output: p-value of t-test

t_test_compareslopes <-  function(data_b1, data_b2) {
  #find b1 and b2 sample sizes
  n_b1 <- length(data_b1[,1])
  n_b2 <- length(data_b2[,1])
  
  #degrees of freedom
  df = n_b1 + n_b2 - 4
  
  #fit linear models to b1 and b2 data
  lm_b1 <- lm(data_b1[,2]~data_b1[,1])
  lm_b2 <- lm(data_b2[,2]~data_b2[,1])
  
  #find b1 and b2 slopes
  slope_b1 <- lm_b1$coefficients[[2]]
  slope_b2 <- lm_b2$coefficients[[2]]
  
  #find standard error of Y given X for the linear regression line for b1 and b2
  steyx_b1 <- steyx(data_b1)
  steyx_b2 <- steyx(data_b2)
  
  #find standard deviation of X for b1 and b2
  sd_b1x <- sd(data_b1[,1])
  sd_b2x <- sd(data_b2[,1])
  
  #pooled error variance approach
  s_res_squared <- ((n_b1-2)*steyx_b1^2 + (n_b2-2)*steyx_b2^2)/(n_b1+n_b2-4)
  
  s_b1b2 <- sqrt(s_res_squared)*sqrt(1/((sd_b1x^2)*(n_b1-1)) + 1/((sd_b2x^2)*(n_b2-1)))
  
  #find t-statistic
  t_stat <- (slope_b1-slope_b2)/(s_b1b2)
  
  #find p-value
  #https://stackoverflow.com/questions/58686432/is-there-an-r-function-similar-to-tdistx-degrees-of-freedom-tails-in-excel
  #Run two-sided t-test 
  #p-value = probability less than -t_stat or greater than t_stat
  #Find the probability less than -t_stat and multiply by 2 to find two-sided p-value
  #pt(x, df) = P(x<t_stat) where T~t(df)
  pval <- 2*pt(-abs(t_stat), df)
  
  return(pval)
}

steyx <-  function(data) {
  # find standard error of Y given X for the 2dimensional data
  # https://support.microsoft.com/en-us/office/steyx-function-6ce74b2c-449d-4a6e-b9ac-f9cef5ba48ab
  x <- data[,1]
  y <- data[,2]
  
  n <- length(x)
  
  #https://rdrr.io/github/Thonnard/THON/src/R/steyx.R 
  output <- sqrt((1/(n-2))*(sum((y-mean(y))^2)-((sum((x-mean(x))*(y-mean(y))))^2)/(sum((x-mean(x))^2))))
  
  #return
  return(output)
}

# #test with real statistics sample dataset
# 
# #change working directory
# setwd('~/Desktop/')
# 
# #import data
# data_b1 <- read.csv('men_lifeexpvscig.csv', header = FALSE)
# data_b2 <- read.csv('women_lifeexpvscig.csv', header = FALSE)