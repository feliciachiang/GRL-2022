# Function to find the best-fit copula family
# Prints "Independence!" if independence copula is chosen

# Written by: Felicia Chiang, felicia.chiang@nasa.gov
# Last updated: 28 Jan 2022

findbestfitcopula <- function(data_x, data_y, familynums = c(1,3:6), familysetind = c(1,3:6,13,14,16,23,24,26,33,34,36), selectionmethod) {
  #import necessary libraries
  library(copula)
  library(VineCopula)
  library(spcopula)
  
  #transform data into uniform marginals
  var_u <- pobs(data_x) 
  var_v <- pobs(data_y) 
  
  #run independence test AND
  #select best fit copula family based on the user defined selectionmethod out of the specified families
  selectedCopula <- BiCopSelect(var_u, var_v, familyset = familynums, selectioncrit = selectionmethod, indeptest = TRUE)

  #pvalue should be greater than 0.05 to be considered a good fit
  #pval will be NaN if independence copula is chosen
  pval <- BiCopGofTest(var_u, var_v, selectedCopula$family, selectedCopula$par, selectedCopula$par2)$p.value

  copfitlist <- list(selectedCopula, pval)

  #evaluate and replace copfitlist if other family has higher pval
  if (is.finite(copfitlist[[2]])) {
    #print("True")
    library(TeachingDemos) #for digits function
    #familysetind = c(1,3:6,13,14,16,23,24,26,33,34,36)
    familysetmat = digits(familysetind, simplify = TRUE)
    
    while(copfitlist[[2]] < 0.05 && length(familysetind) > 1) {
      familysetind = familysetind[familysetmat[2,] != digits(copfitlist[[1]]$family, 2)[2]]
      familysetmat = digits(familysetind, simplify = TRUE)    
      
      #try to look for family with p.value greater than 0.05
      #copfitlist <- findbestfitcopula(data_x, data_y, familysetind, selectionmethod)
      
      #run independence test AND
      #select best fit copula family based on the user defined selectionmethod out of the specified families
      selectedCopula <- BiCopSelect(var_u, var_v, familyset = familysetind, selectioncrit = selectionmethod)
      
      #pvalue should be greater than 0.05 to be considered a good fit
      #pval will be NaN if independence copula is chosen
      pval <- BiCopGofTest(var_u, var_v, selectedCopula$family, selectedCopula$par, selectedCopula$par2)$p.value
      
      copfitlist <- list(selectedCopula, pval)      
      
    }
    #if no other families pass the goodness-of-fit test, revert back to original best fit copula family
    if (copfitlist[[2]] < 0.05) {
      #run independence test AND
      #select best fit copula family based on the user defined selectionmethod out of the specified families
      selectedCopula <- BiCopSelect(var_u, var_v, familyset = familynums, selectioncrit = selectionmethod)
      
      #pvalue should be greater than 0.05 to be considered a good fit
      #pval will be NaN if independence copula is chosen
      pval <- BiCopGofTest(var_u, var_v, selectedCopula$family, selectedCopula$par, selectedCopula$par2)$p.value
      
      copfitlist <- list(selectedCopula, pval)    
    }
  } else {
    print("Independence!")
  }
  
  return(copfitlist)
}