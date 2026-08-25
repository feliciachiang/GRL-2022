# Generate non-exceedance probabilities for each degree of warming for each model 
# Non-exceedance based on the Xth percentile of baseline tasmax or VPD
# 
# Input: Read in copula parameters and gof results, moisture and tasmax/VPD variables
# Output: Jun-Aug non-exceedance probabilities for each degree of warming for each model in csv format (one csv file for each pair of variables)

library(copula)
library(VineCopula)
library(spcopula)
library(ncdf4)
library(EnvStats)
library(kdensity)
library(tools)
#require(pracma)
require(MESS)

#working directory
wdir = "/Users/fchiang/GISS/"
scriptdir = paste0(wdir, "HSM 2022 scripts/")
resultsdir = paste0(wdir, "HSM project/results/")
figdir = paste0(wdir, "HSM project/figures")

#import user defined functions
#source(paste0(scriptdir, "functions/", "make_emp_cpdf.R"))
source(paste0(scriptdir, "functions/", "make_emp_cpdf_smoothed.R"))
source(paste0(scriptdir, "functions/", "findbestfitcopula.R"))

monstrlist <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

#--------------------------------specify---------------------------------------#

#3-month aggreegate
sc <- 3

#region name 
regionname <- 'CNA'

#read in model names
modelnames <- read.table(paste0(wdir, "HSM project/modelnames.txt"), header = FALSE, sep = ",")$V1

#select models with at least 10 ensemble members
#modelinds <- c(2,5,9,10,11)
modelinds <- seq(1,12)


#summer months 
#endmonthinds <- c(6,8,10)
#endmonthinds <- 8

endmonthpos <- 2
endmonth <- 8
savestr <- paste0(monstrlist[endmonth-sc+1], '-', monstrlist[endmonth])

#set percentile to sample 
percentile <- 0.75
#heatpercentile <- 0.75
#heatpercentile <- 0.8
heatpercentile <- 0.9

#set colors for regression lines for each degree of warming
colvector <- c('orange', 'red', 'purple')

yuplim <- 0.85

#------------------------------------------------------------------------------#

#import moisture and heat data: pr and tmax------------------------------------
mvartype <- 'pr'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
maggtype <- 'sum'

vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum Temperature (K)'
aggtype <- 'mean'

#import moisture and heat data: mrsos and tmax---------------------------------
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface soil moisture (kg m-2)'
maggtype <- 'sum'

vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum Temperature (K)'
aggtype <- 'mean'

#import moisture and heat data: pr and vpd------------------------------------
mvartype <- 'pr'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
maggtype <- 'sum'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#import moisture and heat data: mrsos and vpd---------------------------------
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface soil moisture (kg m-2)'
maggtype <- 'sum'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#-----------------------------------------------------------------------------

#import goodness of fit pvalues and copula parameters
#for each model and summer month period 

# copfamily <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4), c(length(endmonthinds), length(modelinds), 4))
# names(dim(copfamily)) <- c('endmonth', 'model', 'warmingdegree')
# gofpval <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4), c(length(endmonthinds), length(modelinds), 4))
# names(dim(gofpval)) <- c('endmonth', 'model', 'warmingdegree')
# copparams <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4*3), c(length(endmonthinds), length(modelinds), 4, 3))
# names(dim(copparams)) <- c('endmonth', 'model', 'warmingdegree', 'params')

#load goodness of fit pvalues
filename=paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-gofpval-JJA-allmodels.nc')

nc <- nc_open(filename)
gofpval <- ncvar_get(nc, nc$var[['gofpvals']])

dendingmonth <- nc$dim$endingmonth$vals
dmodels <- nc$dim$model$vals

nc_close(nc)

#load copparams

filename=paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-copparams-JJA-allmodels.nc')

nc <- nc_open(filename)
copparams <- ncvar_get(nc, nc$var[["copparams"]])  
#dim[[1]] = ending month c(6, 8, 10)
#dim[[2]] = model c(2, 5, 9, 10, 11)
#dim[[3]] = warmingdegree c(0, 1, 2, 3)
#dim[[4]] = parameters
# dendingmonth <- nc$dim[[1]]$vals
# dmodels <- nc$dim[[2]]$vals
# dwarmingdegree <- nc$dim[[3]]$vals
# dparameters <- nc$dim[[4]]$vals
nc_close(nc)

#create array to save non-exceedance values
nonexceeds <- array(rep(NaN, length(modelinds)*4), c(length(modelinds), 4))
names(dim(nonexceeds)) <- c('model', 'warmingdegree')

#for each model
for (modelpos in 1:length(modelinds)) {
  #print(modelpos)
  currentmodel <- modelnames[modelinds[modelpos]]
  
  print(currentmodel)
  
  #----------------------------import data------------------------------------
  
  #import moisture variable data
  nc <- nc_open(paste0(resultsdir, currentmodel, '_', mvartype, '_', maggtype, '_', savestr, '_datafromdegreesofwarming_first5ens.nc'))
  mvardata <- ncvar_get(nc, nc$var[[mvartype]])  
  degrees <- nc$dim$degree$vals
  nc_close(nc)
  
  #import baseline 1850-1899 data
  nc <- nc_open(paste0(resultsdir, currentmodel, '_', mvartype, '_', maggtype, '_', savestr, '_1850-1899baseline_first5ens.nc'))
  mvarbaseline <- ncvar_get(nc, nc$var[[mvartype]])
  nc_close(nc)
  
  #import heat variable data
  nc <- nc_open(paste0(resultsdir, currentmodel, '_', vartype, '_', aggtype, '_', savestr, '_datafromdegreesofwarming_first5ens.nc'))
  vardata <- ncvar_get(nc, nc$var[[vartype]])
  nc_close(nc)
  
  #import baseline 1850-1899 data for the relevant season
  nc <- nc_open(paste0(resultsdir, currentmodel, '_', vartype, '_', aggtype, '_', savestr, '_1850-1899baseline_first5ens.nc'))
  varbaseline <- ncvar_get(nc, nc$var[[vartype]])
  nc_close(nc)
  
  
  #------------------------------------------------------------------------------
  #find seasonal averages from the baseline (1850-1899)
  baseline_seasonalmean <- mean(varbaseline)
  
  #anomalize data according to seasonal average from the baseline from all ensemble members
  vardata <- vardata - baseline_seasonalmean
  
  varbaseline <- varbaseline-baseline_seasonalmean
  
  #find 75th percentile line from the vardata (temperature or VPD) baseline time period of 1850-1899
  varbaselinethreshold <- quantile(varbaseline, heatpercentile)
  
  #use negative precipitation in order to have positive correlation for x and y
  data_baseline <- data.frame(as.vector(varbaseline), as.vector(-mvarbaseline))
  data_1deg <- data.frame(as.vector(vardata[,,1]), as.vector(-mvardata[,,1]))
  data_2deg <- data.frame(as.vector(vardata[,,2]), as.vector(-mvardata[,,2]))
  data_3deg <- data.frame(as.vector(vardata[,,3]), as.vector(-mvardata[,,3]))
  #data_4deg <- data.frame(as.vector(vardata[,,4]), as.vector(-mvardata[,,4]))
  
  #combine all data to create comprehensive data frame
  var_all = c(as.vector(varbaseline), as.vector(vardata))
  #find the x and y limits given all data
  xlim = c(min(var_all), max(var_all))
  
  #import baseline 1850-1899 data for the relevant season
  baselinedata <- mvarbaseline
  
  #75th percentile based on data from 1850-1899 (also negative)
  baseline75p <- quantile(c(baselinedata), percentile, na.rm = TRUE)
  baseline75p <- -baseline75p
  
  #baseline
  fittedcopulabaseline <- copulaFromFamilyIndex(copparams[modelpos, 1, 1], copparams[modelpos, 1, 2], copparams[modelpos, 1, 3])
  #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
  cond_dens_empbase <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopulabaseline, data_baseline[,1], data_baseline[,2]))

  
  #1 degree of warming
  fittedcopula1 <- copulaFromFamilyIndex(copparams[modelpos, 2, 1], copparams[modelpos, 2, 2], copparams[modelpos, 2, 3])
  #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
  cond_dens_emp1 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula1, data_1deg[,1], data_1deg[,2]))
  plot(cond_dens_emp1, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))
  
  #create empirical estimation of marginal distribution of temperature
  #make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
  
  #2 degrees of warming
  fittedcopula2 <- copulaFromFamilyIndex(copparams[modelpos, 3, 1], copparams[modelpos, 3, 2], copparams[modelpos, 3, 3])
  #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
  cond_dens_emp2 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula2, data_2deg[,1], data_2deg[,2]))
  plot(cond_dens_emp2, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))
  
  
  #create empirical estimation of marginal distribution of temperature
  #make_marg_emp2 <-  function(x) demp(x, as.vector(data_2deg[,1]))
  
  #3 degrees of warming
  fittedcopula3 <- copulaFromFamilyIndex(copparams[modelpos, 4, 1], copparams[modelpos, 4, 2], copparams[modelpos, 4, 3])
  #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
  cond_dens_emp3 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula3, data_3deg[,1], data_3deg[,2]))
  plot(cond_dens_emp3, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))

  #visualize the conditional PDFs
  plot(cond_dens_empbase, from=xlim[1], to=xlim[2], ylim = c(0,yuplim), xlab = paste0('Conditional ', varstring, ' distribution'), ylab = 'Density', main = paste0(currentmodel, ' ', savestr), col = 'black')
  plot(cond_dens_emp1, from=xlim[1], to=xlim[2], col = 'orange', ylim = c(0,yuplim), add=TRUE)
  plot(cond_dens_emp2, from=xlim[1], to=xlim[2], col = 'red', ylim = c(0,yuplim), add=TRUE)
  plot(cond_dens_emp3, from=xlim[1], to=xlim[2], col = 'purple', ylim = c(0,yuplim), add=TRUE)
  #add vline for var baseline threshold
  abline(v=varbaselinethreshold, col = 'black', lty=2)
  
  #find the x and y coordinates associated with each cPDF to use MESS::auc to estimate the area under the curve
  cond_dens_empbase_x <- seq(from = xlim[1], to = xlim[2], length.out = 1000)
  cond_dens_empbase_y <- cond_dens_empbase(cond_dens_empbase_x)
  
  cond_dens_emp1_x <- seq(from = xlim[1], to = xlim[2], length.out = 1000)
  cond_dens_emp1_y <- cond_dens_emp1(cond_dens_emp1_x)
  
  cond_dens_emp2_x <- seq(from = xlim[1], to = xlim[2], length.out = 1000)
  cond_dens_emp2_y <- cond_dens_emp2(cond_dens_emp2_x)
  
  cond_dens_emp3_x <- seq(from = xlim[1], to = xlim[2], length.out = 1000)
  cond_dens_emp3_y <- cond_dens_emp3(cond_dens_emp3_x)  
  
  #find area under each conditional curve
  
  nonexceeds[modelpos, 1] <- MESS::auc(cond_dens_empbase_x[cond_dens_empbase_x < varbaselinethreshold], cond_dens_empbase_y[cond_dens_empbase_x < varbaselinethreshold], type = 'spline')/MESS::auc(cond_dens_empbase_x, cond_dens_empbase_y, type = 'spline')
  nonexceeds[modelpos, 2] <- MESS::auc(cond_dens_emp1_x[cond_dens_emp1_x < varbaselinethreshold], cond_dens_emp1_y[cond_dens_emp1_x < varbaselinethreshold], type = 'spline')/MESS::auc(cond_dens_emp1_x, cond_dens_emp1_y, type = 'spline')
  nonexceeds[modelpos, 3] <- MESS::auc(cond_dens_emp2_x[cond_dens_emp2_x < varbaselinethreshold], cond_dens_emp2_y[cond_dens_emp2_x < varbaselinethreshold], type = 'spline')/MESS::auc(cond_dens_emp2_x, cond_dens_emp2_y, type = 'spline')
  nonexceeds[modelpos, 4] <- MESS::auc(cond_dens_emp3_x[cond_dens_emp3_x < varbaselinethreshold], cond_dens_emp3_y[cond_dens_emp3_x < varbaselinethreshold], type = 'spline')/MESS::auc(cond_dens_emp3_x, cond_dens_emp3_y, type = 'spline')
  
}

#save nonexceeds as csv
write.csv(nonexceeds, paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-nonexceeding-baseline75th-JJA-allmodels_wbaseline_', as.character(heatpercentile*100), 'th.csv'), row.names= FALSE)
    
