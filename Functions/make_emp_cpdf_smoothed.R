#Function to create conditional PDFs based on:
#   1) fitted copula and copula parameters
#   2) raw x and y vectors

#Return: function to find conditional PDFs based on y value

make_emp_cpdf_smoothed <-  function(y, fittedcopula, data_x, data_y) {
  library(copula) #import function dCopula to calculate conditional PDF
  library(spatstat.core) 
  library(kdensity) #import function kdensity to generate kernel density (PDF) estimate for univariate data
  library(EnvStats)
  
  #establish bounds for x variable
  xminfrom = min(data_x)-(range(data_x)[2]-range(data_x)[1])*0.25
  xmaxfrom = max(data_x)+(range(data_x)[2]-range(data_x)[1])*0.25
  
  #create empirical cdf function
  #density function finds the density of the x variable from the min and max x values of the variable
  bx <- density(data_x, from = xminfrom, to = xmaxfrom, n=1e5)
  #CDF function creates the cumulative distribution function whose probability density is estimated as b_x 
  fx <- CDF(bx)
  
  #establish bounds for y variable
  yminfrom = min(data_y)-(range(data_y)[2]-range(data_y)[1])*0.25
  ymaxfrom = max(data_y)+(range(data_y)[2]-range(data_y)[1])*0.25
  
  #create empirical cdf function
  #density function finds the density of the y variable from the min and max x values of the variable
  by <- density(data_y, from = yminfrom, to = ymaxfrom, n=1e5)
  #CDF function creates the cumulative distribution function whose probability density is estimated as b_y
  fy <- CDF(by)  
  
  #define density function for x variable
  #kdex = kdensity(data_x)
  
  #return(function(x) dCopula(cbind(fx(x), fy(y)), fittedcopula)*kdex(x))
  return(function(x) dCopula(cbind(fx(x), fy(y)), fittedcopula)*demp(x, data_x))
}

