# Bivariate data visualization
#
# Import: 3-month warm season data corresponding to each degree (1, 2, 3, 4) of warming
#         -Apr-June, June-Aug, Aug-Oct
# Output: Plots of bivariate data corresponding to each degree of warming

#import libraries
library(copula)
library(VineCopula)
library(spcopula)
library(ncdf4)
library(EnvStats)
library(kdensity)
library(tools)

#libraries
library(ggplot2)
library(ggridges)

#working directory
wdir = "/Users/fchiang/GISS/"
scriptdir = paste0(wdir, "HSM 2022 scripts/")
resultsdir = paste0(wdir, "HSM project/results/")
figdir = paste0(wdir, "HSM project/figures")

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
endmonthinds <- 8

#set colors for regression lines for each degree of warming
colvector <- c('orange', 'red', 'purple')

#------------------------------------------------------------------------------#

#import moisture data: pr------------------------------------------
# mvartype <- 'pr'
# mvarstring <- 'Pr'
# mvaraxlab <- 'Precipitation (kg m-2 s-1)'
# maggtype <- 'sum'

#import mrsos-------------------------------------------------------
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface soil moisture (kg m-2)'
maggtype <- 'sum'

#import placeholder heat data: tasmax-------------------------------------------------------
vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum Temperature (K)'
aggtype <- 'mean'

#import heat data: vpd----------------------------------------------------------
# vartype <- 'vpd'
# varstring <- 'VPD'
# varaxlab <- 'Vapor pressure deficit (mb)'
# aggtype <- 'mean'

#combine data from each model
#dataarray <- array(rep(NaN, 3*length(modelinds)*4*100), c(length(modelinds)*4*100, 3))
#colnames(dataarray) <- c('model', 'degree', 'value')

dataarray <- array(rep(NaN, length(modelinds)*550*3), c(length(modelinds)*550, 3))
colnames(dataarray) <- c('model', 'degree', 'value')
dataframe <- data.frame(dataarray)

indexval <- 1

for (modelpos in 1:length(modelinds)) {
  #print(modelpos)
  currentmodel <- modelnames[modelinds[modelpos]]
  
  print(currentmodel)
  
  #for each summer period
  for (endmonthpos in 1:length(endmonthinds)) {
    endmonth <- endmonthinds[endmonthpos]
    savestr <- paste0(monstrlist[endmonth-sc+1], '-', monstrlist[endmonth])
    
    print(savestr)
    
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
    baseline_seasonalmean <- mean(mvarbaseline)
    
    #anomalize data according to seasonal average from the baseline from all ensemble members
    mvardata <- mvardata - baseline_seasonalmean
    
    mvarbaseline <- mvarbaseline-baseline_seasonalmean
    
    #find 75th percentile line from the vardata (temperature or VPD) baseline time period of 1850-1899
    varbaselinethreshold <- quantile(varbaseline, percentile)
    
    #use negative precipitation in order to have positive correlation for x and y
    data_baseline <- data.frame(as.vector(varbaseline), as.vector(mvarbaseline))
    data_1deg <- data.frame(as.vector(vardata[,,1]), as.vector(mvardata[,,1]))
    data_2deg <- data.frame(as.vector(vardata[,,2]), as.vector(mvardata[,,2]))
    data_3deg <- data.frame(as.vector(vardata[,,3]), as.vector(mvardata[,,3]))
    #data_4deg <- data.frame(as.vector(vardata[,,4]), as.vector(-mvardata[,,4]))
    
    dataframe[indexval:(indexval+549), 1] <- currentmodel
    
    dataframe[indexval:(indexval+249),2] <- '0'
    dataframe[(indexval+250):(indexval+250+99),2] <- '1'
    dataframe[(indexval+250+100):(indexval+250+199),2] <- '2'
    dataframe[(indexval+250+200):(indexval+250+299),2] <- '3'
    
    dataframe[indexval:(indexval+249),3] <- data_baseline[,2]
    dataframe[(indexval+250):(indexval+250+99),3] <- data_1deg[,2]
    dataframe[(indexval+250+100):(indexval+250+199),3] <- data_2deg[,2]
    dataframe[(indexval+250+200):(indexval+250+299),3] <- data_3deg[,2]
    
    indexval <- indexval + 550
  }
}

dataframe$model = factor(dataframe$model, levels = modelnames)
dataframe$model = with(dataframe, factor(model, levels = rev(levels(model))))

jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_ridge.jpeg'), width = 4, height = 8, units = 'in', res = 300)
ggplot(dataframe, aes(y=as.factor(model),
                      x=value,
                      fill=degree)) +
  geom_density_ridges(alpha=0.25, scale = 0.95) +
  scale_y_discrete(expand = c(0.01, 0)) +  
  scale_x_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c("green", "orange", "red", "purple"))
dev.off()



