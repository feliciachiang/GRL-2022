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
library(dplyr)

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
vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum temperature (K)'
aggtype <- 'mean'

#import heat data: vpd----------------------------------------------------------
# vartype <- 'vpd'
# varstring <- 'VPD'
# varaxlab <- 'Vapor pressure deficit (mb)'
# aggtype <- 'mean'

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
# vartype <- 'hfss'
# varaxlab <- 'Sensible heat flux (W m-2)'
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
vartype <- 'ea'
varaxlab <- 'Vapor pressure (mb)'
aggtype <- 'mean'
# 
vartype <- 'es'
varaxlab <- 'Saturated vapor pressure (mb)'
aggtype <- 'mean'


#create grid to plot on:
#jpeg(paste0(figdir, '/CMIP6multimodel_', vartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_demp_col_legend.jpeg'), width = 10/1.4, height = 12/1.4, units = 'in', res = 300)

jpeg(paste0(figdir, '/CMIP6multimodel_', vartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_demp_col_withyaxis.jpeg'), width = 3.5/1.4, height = 12/1.4, units = 'in', res = 300)
par(mfcol = c(12, 1), mar = numeric(4), oma = c(4,4,.5,.5), mgp = c(2, .6, 0))

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
    
    #import baseline 1850-1899 data for the relevant season
    baselinedata <- mvarbaseline
    
    #75th percentile based on data from 1850-1899 (also negative)
    baseline75p <- quantile(c(baselinedata), percentile, na.rm = TRUE)
    baseline75p <- -baseline75p
    
    #combine all data to create comprehensive data frame
    var_all = c(as.vector(varbaseline), as.vector(vardata))
    #find the x and y limits given all data
    xlim = c(min(var_all), max(var_all))

    #create empirical estimation of marginal distribution of variable
    make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,1]))
    #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
    
    
    
    #1 degree of warming
    make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
    #make_marg_emp1 <- suppressWarnings(kdensity(as.vector(data_1deg[,1])))
    
    #2 degrees of warming
    make_marg_emp2 <-  function(x) demp(x, as.vector(data_2deg[,1]))
    #make_marg_emp2 <- suppressWarnings(kdensity(as.vector(data_2deg[,1])))
    
    #3 degrees of warming
    make_marg_emp3 <-  function(x) demp(x, as.vector(data_3deg[,1]))
    #make_marg_emp3 <- suppressWarnings(kdensity(as.vector(data_3deg[,1])))
    
    #--------------------------------PLOTS!----------------------------------------#
    
    #establish x lim (standard across models)
    xlowlim <- -4   #for tasmax
    xhilim <- 10    #for tasmax
    
    #xlowlim <- -35   #for hfls
    #xhilim <- 30     #for hfls
    
    #xlowlim <- -25   #for hfss
    #xhilim <- 35     #for hfss
    
    #xlowlim <- -10   #for VPD
    #xhilim <- 20     #for VPD
    
    #xlowlim <- -4     #for ea
    #xhilim <- 8       #for ea
    
    #xlowlim <- -10     #for es
    #xhilim <- 20      #for es
    
    #set upper limit for y axis
    #yuplim <- 0.65
    yuplim <- 0.65 #-> tasmax, vpd, tas
    #yuplim <- 0.2 #-> hfls
    #yuplim <- 10 #-> bo
    #yuplim <- 0.25 #-> hurs
    #yuplim <- 1200 #-> huss
    #yuplim <- 1 #-> ea
    
    #establish y lim (not standard across models)
    #account for highest value in all curves
    marg_base_ymax <- max(make_marg_empbase(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp1_ymax <- max(make_marg_emp1(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp2_ymax <- max(make_marg_emp2(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp3_ymax <- max(make_marg_emp3(seq(xlowlim, xhilim, length.out = 1000)))
    
    ylowlim <- 0
    yhilim <- max(marg_base_ymax, marg_emp1_ymax, marg_emp2_ymax, marg_emp3_ymax)
    yhilim <- yhilim*1.15
    
    
    #plot empirical marginal temperature/VPD/alt distributions for baseline, 1-3 degrees of warming
    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', yaxt = 'n')
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7", border = NA)
    #for tasmax
    if (yhilim > 0.7) {
      axis(2, at = c(0, 0.25, 0.5, 1), labels = c(0, 0.25, 0.5, 1))
    } else {
      print('confirm')
      axis(2, at = c(0, 0.25, 1), labels = c(0, 0.25, 1))
    }

    #for vpd
    #axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1))
    
    # #for ea
    # if (between(yhilim, 0.45, 0.55)) {
    #   axis(2, at = c(0, 0.25, 1), labels = c(0, 0.25, 1))
    # } else {
    #   axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1))
    # }
    
    #for es
    # if (between(yhilim, 0.25, 0.35)) {
    #   axis(2, at = c(0, 0.15, 1), labels = c(0, 0.15, 1))
    # } else {
    #   axis(2, at = c(0, 0.15, 0.3, 0.45, 1), labels = c(0, 0.15, 0.3, 0.45, 1))
    # }


    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', add = TRUE)
    title(main = paste0(currentmodel), line = 0.5)
    plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
    plot(make_marg_emp2, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'red', add=TRUE)
    plot(make_marg_emp3, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'purple', add=TRUE)
    
    # if (modelpos < length(modelinds)) {
    #   axis(side = 1, pos = 0, tck = 0, labels = F)
    # } else {
    #   axis(side = 1, pos = 0)
    # }
    
    axis(side = 1, pos = 0, tck = 0, labels = F)
    
    # plot(make_marg_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, xlim = c(min(data_baseline[,1]), max(data_3deg[,1])), ylim = c(0,yuplim), xlab = 'Maximum Temperature (K)', ylab = 'Density',  main = paste0(currentmodel, ' ', savestr), col = 'black')
    # plot(make_marg_emp1, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'orange', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp2, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'red', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp3, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'purple', ylim = c(0,yuplim), add=TRUE)
    
  }
}
#legend(x = 'bottom', inset = c(0, -0.13), legend = c('Baseline', '1 Degree', '2 Degrees', '3 Degrees'), col = c('green', 'orange', 'red', 'purple'), ncol = 4, lty = 1, xpd = TRUE)
#save plot
title(xlab = varaxlab,
      ylab = 'Density',
      outer = TRUE, line = 2.8, cex.lab = 1.5)
dev.off()

#for EF
vartype <- 'EF'
varaxlab <- 'Evaporative Fraction'
aggtype <- 'mean'

jpeg(paste0(figdir, '/CMIP6multimodel_', vartype, '-marginalPDFs_withbaseline_anom-JJA-allmodels_demp_col.jpeg'), width = 3.5/1.4, height = 12/1.4, units = 'in', res = 300)
par(mfcol = c(12, 1), mar = numeric(4), oma = c(4,4,.5,.5), mgp = c(2, .6, 0))

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
    
    lhvartype <- 'hfls'
    
    #import heat variable data
    nc <- nc_open(paste0(resultsdir, currentmodel, '_', lhvartype, '_', aggtype, '_', savestr, '_datafromdegreesofwarming_first5ens.nc'))
    lhvardata <- ncvar_get(nc, nc$var[[lhvartype]])
    nc_close(nc)
    
    #import baseline 1850-1899 data for the relevant season
    nc <- nc_open(paste0(resultsdir, currentmodel, '_', lhvartype, '_', aggtype, '_', savestr, '_1850-1899baseline_first5ens.nc'))
    lhvarbaseline <- ncvar_get(nc, nc$var[[lhvartype]])
    nc_close(nc)
    
    shvartype <- 'hfss'
    
    #import heat variable data
    nc <- nc_open(paste0(resultsdir, currentmodel, '_', shvartype, '_', aggtype, '_', savestr, '_datafromdegreesofwarming_first5ens.nc'))
    shvardata <- ncvar_get(nc, nc$var[[shvartype]])
    nc_close(nc)    
    
    #import baseline 1850-1899 data for the relevant season
    nc <- nc_open(paste0(resultsdir, currentmodel, '_', shvartype, '_', aggtype, '_', savestr, '_1850-1899baseline_first5ens.nc'))
    shvarbaseline <- ncvar_get(nc, nc$var[[shvartype]])
    nc_close(nc)    
    
    #calculate EF data, where EF = lh/(lh + sh)
    vardata <- lhvardata/(lhvardata + shvardata)
    varbaseline <- lhvarbaseline/(lhvarbaseline + shvarbaseline)
    
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
    
    #import baseline 1850-1899 data for the relevant season
    baselinedata <- mvarbaseline
    
    #75th percentile based on data from 1850-1899 (also negative)
    baseline75p <- quantile(c(baselinedata), percentile, na.rm = TRUE)
    baseline75p <- -baseline75p
    
    #combine all data to create comprehensive data frame
    var_all = c(as.vector(varbaseline), as.vector(vardata))
    #find the x and y limits given all data
    xlim = c(min(var_all), max(var_all))
    
    #create empirical estimation of marginal distribution of variable
    make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,1]))
    #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
    
    
    
    #1 degree of warming
    make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
    #make_marg_emp1 <- suppressWarnings(kdensity(as.vector(data_1deg[,1])))
    
    #2 degrees of warming
    make_marg_emp2 <-  function(x) demp(x, as.vector(data_2deg[,1]))
    #make_marg_emp2 <- suppressWarnings(kdensity(as.vector(data_2deg[,1])))
    
    #3 degrees of warming
    make_marg_emp3 <-  function(x) demp(x, as.vector(data_3deg[,1]))
    #make_marg_emp3 <- suppressWarnings(kdensity(as.vector(data_3deg[,1])))
    
    #--------------------------------PLOTS!----------------------------------------#
    
    #establish x lim (standard across models)
    xlowlim <- -0.25   #for EF
    xhilim <- 0.15      
    
    #set upper limit for y axis
    yuplim <- 10       #-> EF
    
    #establish y lim (not standard across models)
    #account for highest value in all curves
    marg_base_ymax <- max(make_marg_empbase(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp1_ymax <- max(make_marg_emp1(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp2_ymax <- max(make_marg_emp2(seq(xlowlim, xhilim, length.out = 1000)))
    marg_emp3_ymax <- max(make_marg_emp3(seq(xlowlim, xhilim, length.out = 1000)))
    
    ylowlim <- 0
    yhilim <- max(marg_base_ymax, marg_emp1_ymax, marg_emp2_ymax, marg_emp3_ymax)
    yhilim <- yhilim*1.15
    
    
    #plot empirical marginal temperature/VPD/alt distributions for baseline, 1-3 degrees of warming
    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', axes = FALSE)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7", border = NA)
    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', add = TRUE)
    title(main = paste0(currentmodel), line = 0.5)
    plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
    plot(make_marg_emp2, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'red', add=TRUE)
    plot(make_marg_emp3, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'purple', add=TRUE)
    
    if (modelpos < length(modelinds)) {
      axis(side = 1, pos = 0, tck = 0, labels = F)
    } else {
      axis(side = 1, pos = 0)
    }
    
    # plot(make_marg_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, xlim = c(min(data_baseline[,1]), max(data_3deg[,1])), ylim = c(0,yuplim), xlab = 'Maximum Temperature (K)', ylab = 'Density',  main = paste0(currentmodel, ' ', savestr), col = 'black')
    # plot(make_marg_emp1, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'orange', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp2, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'red', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp3, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'purple', ylim = c(0,yuplim), add=TRUE)
    
  }
}
#legend(x = 'bottom', inset = c(0, -0.13), legend = c('Baseline', '1 Degree', '2 Degrees', '3 Degrees'), col = c('green', 'orange', 'red', 'purple'), ncol = 4, lty = 1, xpd = TRUE)
#save plot
title(xlab = varaxlab,
      ylab = 'Density',
      outer = TRUE, line = 2.8, cex.lab = 1.5)
dev.off()

