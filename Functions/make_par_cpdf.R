#Function to create conditional PDFs based on:
#   1) fitted copula and copula parameters
#   2) raw x and y vectors

#Return: function to find conditional PDFs based on y value

make_par_cpdf <-  function(y, fittedcopula, data_x, data_y) {
  library(copula)
  library(fitdistrplus)
  
  #find marginal distributions for x and positive y
  fittedx <- fitdist(data_x, 'norm')
  fittedy <- fitdist(data_y, 'gamma')
  
  #check fit
  
  return(function(x) dCopula(cbind(pnorm(x, mean = fittedx$estimate[[1]], sd = fittedx$estimate[[2]]), pgamma(y, shape = fittedy$estimate[[1]], rate = fittedy$estimate[[2]])), fittedcopula)*dnorm(x, mean = fittedx$estimate[[1]], sd = fittedx$estimate[[2]]))
}
