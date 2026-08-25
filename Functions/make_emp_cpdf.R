#Function to create conditional PDFs based on:
#   1) fitted copula and copula parameters
#   2) raw x and y vectors

#Return: function to find conditional PDFs based on y value

make_emp_cpdf <-  function(y, fittedcopula, data_x, data_y) {
  library(copula)
  library(EnvStats)
  
  return(function(x) dCopula(cbind(pemp(x, data_x), pemp(y, data_y)), fittedcopula)*demp(x, data_x))
}

