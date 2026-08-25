# Find "signal" to "noise" ratio for heat marginals 
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

library(tidyr)
library(ggplot2)
library(plotly)

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

#import placeholder moisture data: pr------------------------------------------
mvartype <- 'pr'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
maggtype <- 'sum'

#import heat data: tasmax-------------------------------------------------------
# vartype <- 'tasmax'
# varstring <- 'temperature'
# varaxlab <- 'Maximum temperature (K)'
# aggtype <- 'mean'

#import heat data: vpd----------------------------------------------------------
vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#import heat data: tas----------------------------------------------------------
# vartype <- 'tas'
# varstring <- 'temperature'
# varaxlab <- 'Temperature (K)'
# aggtype <- 'mean'

#additional variables----------------------------------------

# #latent heat flux (hfls), Bowen ratio (bo), relative humidity (hurs), specific humidity (huss), actual vapor pressure (ea)
# vartype <- 'hfls'
# varaxlab <- 'Latent heat flux (W m-2)'
# aggtype <- 'mean'
# 
# vartype <- 'bo'
# varaxlab <- 'Bowen ratio (1)'
# aggtype <- 'mean'
# 
# vartype <- 'hurs'
# varaxlab <- 'Relative humidity (%)'
# aggtype <- 'mean'
# 
# vartype <- 'huss'
# varaxlab <- 'Specific humidity (1)'
# aggtype <- 'mean'
# 
# vartype <- 'ea'
# varaxlab <- 'Vapor pressure (mb)'
# aggtype <- 'mean'

# vartype <- 'es'
# varaxlab <- 'Saturated vapor pressure (mb)'
# aggtype <- 'mean'


#create grid to plot on:
#jpeg(paste0(figdir, '/CMIP6multimodel_', vartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_demp_col_legend.jpeg'), width = 10/1.4, height = 12/1.4, units = 'in', res = 300)

# jpeg(paste0(figdir, '/CMIP6multimodel_', vartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_demp_col.jpeg'), width = 3.5/1.4, height = 12/1.4, units = 'in', res = 300)
# par(mfcol = c(12, 1), mar = numeric(4), oma = c(4,4,.5,.5), mgp = c(2, .6, 0))

snr <- array(rep(NaN, length(modelinds)), c(length(modelinds)))

#for each model
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
    baseline_seasonalmean <- mean(varbaseline)
    
    #anomalize data according to seasonal average from the baseline from all ensemble members
    vardata <- vardata - baseline_seasonalmean
    
    varbaseline <- varbaseline-baseline_seasonalmean
    
    #find 75th percentile line from the vardata (temperature or VPD) baseline time period of 1850-1899
    varbaselinethreshold <- quantile(varbaseline, percentile)
    
    #use negative precipitation in order to have positive correlation for x and y
    data_baseline <- data.frame(as.vector(varbaseline), as.vector(-mvarbaseline))
    data_1deg <- data.frame(as.vector(vardata[,,1]), as.vector(-mvardata[,,1]))
    data_2deg <- data.frame(as.vector(vardata[,,2]), as.vector(-mvardata[,,2]))
    data_3deg <- data.frame(as.vector(vardata[,,3]), as.vector(-mvardata[,,3]))
    #data_4deg <- data.frame(as.vector(vardata[,,4]), as.vector(-mvardata[,,4]))
    
    #find snr using mean shift from baseline to 2nd degree of warming and sd from baseline
    snr[modelpos] <- (mean(data_2deg[,1])-mean(data_baseline[,1]))/(sd(data_baseline[,1]))
    #snr[modelpos] <- (median(data_2deg[,1])-median(data_baseline[,1]))/(sd(data_baseline[,1]))
    #snr[modelpos] <- (median(data_2deg[,1])-median(data_baseline[,1]))/(IQR(data_baseline[,1]))
    #snr[modelpos] <- (quantile(data_2deg[,1], 0.1)-quantile(data_baseline[,1], 0.1))/sd(data_baseline[,1])
  }
}


#convert csv table into dataframe
df <- data.frame(snr, row.names = modelnames)
colnames(df) <- c('SNR')
df_models <- cbind(Model = rownames(df), df)
rownames(df_models) <- NULL

dat.g <- gather(df_models, type, value, -Model)

#create nested barplot
jpeg(paste0(figdir, '/CMIP6multimodel_',  vartype, '-snr-iqr-JJA-allmodels-barplot.jpeg'), width = 8, height = 4, units = 'in', res = 300)
g <- ggplot(dat.g, aes(type, value)) + geom_bar(aes(fill = Model), stat = 'identity', position = 'dodge') + labs(y = 'Signal to noise ratio') #+ ylim(0, 5)
print(g)
dev.off()

#export equivalent percentiles information
write.csv(snr, file = paste0(resultsdir, '/CMIP6multimodel_', vartype, '-snr-iqr-JJA-allmodels.csv'), row.names = modelnames)


