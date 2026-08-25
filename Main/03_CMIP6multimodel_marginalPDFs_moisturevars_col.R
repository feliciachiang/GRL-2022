# Moisture data marginal distribution visualization
#
# Import: 3-month warm season data corresponding to each degree (1, 2, 3, 4) of warming
#         -Apr-June, June-Aug, Aug-Oct
# Output: Plots of moisture data corresponding to each degree of warming

#import libraries
library(copula)
library(VineCopula)
library(spcopula)
library(ncdf4)
library(EnvStats)
library(kdensity)
library(tools)

#working directory
wdir = "/Users/fchiang/GISS/"
scriptdir = paste0(wdir, "HSM 2022 scripts/")
resultsdir = paste0(wdir, "HSM project/results/")
figdir = paste0(wdir, "HSM project/figures")

monstrlist <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

#--------------------------------specify---------------------------------------#

#percentile
percentile <- 0.75

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
# mvaraxlab <- 'Precipitation (mm/day)'
# maggtype <- 'sum'

#import mrsos 
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface soil moisture (mm)'
maggtype <- 'sum'


#import placeholder heat data: tasmax-------------------------------------------------------
vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum Temperature (K)'
aggtype <- 'mean'

# equivpercent <- array(rep(NaN, length(modelinds)*length(endmonthinds)*3), c(length(endmonthinds), length(modelinds), 3))
# names(dim(copfamily)) <- c('endmonth', 'model', 'warmingdegree')

equivpercent <- array(rep(NaN, length(endmonthinds)*3), c(length(modelinds), 3))
names(dim(copfamily)) <- c('model', 'warmingdegree')

#create grid to plot on:
# jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, '-marginalPDFs_withbaseline_anom_demp.jpeg'), width = 8*2, height = 4*2, units = 'in', res = 300)
# par(mfcol = c(length(endmonthinds), length(modelinds)))

jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, '-marginalPDFs_withbaseline_anom_demp_JJA-allmodels_col.jpeg'), width = 3.5/1.4, height = 12/1.4, units = 'in', res = 300)
#par(mfcol = c(length(modelinds), 1))
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
    
    #for precipitation, convert kg/m2/s == mm/s to mm/day
    #mvardata <- mvardata*24*60*60
    #mvarbaseline <- mvarbaseline*24*60*60
    
    #for soil moisture, kg/m2 == mm
    
    #------------------------------------------------------------------------------
    #find seasonal averages from the baseline (1850-1899)
    baseline_seasonalmean <- mean(mvarbaseline)
    
    #anomalize data according to seasonal average from the baseline from all ensemble members
    mvardata <- mvardata - baseline_seasonalmean
    
    mvarbaseline <- mvarbaseline-baseline_seasonalmean
    
    #find 75th percentile line from the mvardata (precipitation or soil moisture) baseline time period of 1850-1899
    mvarbaselinethreshold <- quantile(mvarbaseline, percentile, na.rm = TRUE)
    
    #use negative precipitation in order to have positive correlation for x and y
    data_baseline <- data.frame(as.vector(varbaseline), as.vector(mvarbaseline))
    data_1deg <- data.frame(as.vector(vardata[,,1]), as.vector(mvardata[,,1]))
    data_2deg <- data.frame(as.vector(vardata[,,2]), as.vector(mvardata[,,2]))
    data_3deg <- data.frame(as.vector(vardata[,,3]), as.vector(mvardata[,,3]))
    #data_4deg <- data.frame(as.vector(vardata[,,4]), as.vector(-mvardata[,,4]))
    
    #print out what percentile the baseline threshold is at for each degree of warming
    #print(mvarbaselinethreshold*1.05)
    #print(max(data_1deg[,2]))
    #print(max(data_2deg[,2]))
    #print(max(data_3deg[,2]))
    
    #save the ecdf value for each degree of warming
    equivpercent[modelpos, 1] <- ecdf(data_1deg[,2])(mvarbaselinethreshold)
    equivpercent[modelpos, 2] <- ecdf(data_2deg[,2])(mvarbaselinethreshold)
    equivpercent[modelpos, 3] <- ecdf(data_3deg[,2])(mvarbaselinethreshold)
    
    #combine all data to create comprehensive data frame
    var_all = c(as.vector(mvarbaseline), as.vector(mvardata))
    #find the x and y limits given all data
    xlim = c(min(var_all), max(var_all))

    #create empirical estimation of marginal distribution of variable
    make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,2]))
    #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
    
    
    
    #1 degree of warming
    make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,2]))
    #make_marg_emp1 <- suppressWarnings(kdensity(as.vector(data_1deg[,1])))
    
    #2 degrees of warming
    make_marg_emp2 <-  function(x) demp(x, as.vector(data_2deg[,2]))
    #make_marg_emp2 <- suppressWarnings(kdensity(as.vector(data_2deg[,1])))
    
    #3 degrees of warming
    make_marg_emp3 <-  function(x) demp(x, as.vector(data_3deg[,2]))
    #make_marg_emp3 <- suppressWarnings(kdensity(as.vector(data_3deg[,1])))
    
    #--------------------------------PLOTS!----------------------------------------#
    
    #establish x lim (standard across models)
    #xlowlim <- -6 #for pr
    #xhilim <- 6   #for pr
    
    xlowlim <- -30 #for mrsos
    xhilim <- 20 #for mrsos
    
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
    #plot(make_marg_empbase, from=xlim[1], to=xlim[2], xlim = c(xlim[1], xlim[2]), ylim = c(0,yuplim), xlab = mvaraxlab, ylab = 'Density',  main = paste0(currentmodel, ' ', savestr), col = 'black')
    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', axes = FALSE)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7", border = NA)
    plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', add = TRUE)
    
    title(main = paste0(currentmodel), line = 0.5)
    plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
    plot(make_marg_emp2, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'red', add=TRUE)
    plot(make_marg_emp3, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'purple', add=TRUE)
    abline(v=mvarbaselinethreshold, col = 'black', lty=2)
    #points(xpoints, ypoints, pch=8, col= c("orange", "red", "purple"), ylim = c(0,yuplim), add = TRUE)
    
    # plot(make_marg_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, xlim = c(min(data_baseline[,1]), max(data_3deg[,1])), ylim = c(0,yuplim), xlab = 'Maximum Temperature (K)', ylab = 'Density',  main = paste0(currentmodel, ' ', savestr), col = 'black')
    # plot(make_marg_emp1, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'orange', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp2, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'red', ylim = c(0,yuplim), add=TRUE)
    # plot(make_marg_emp3, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, col = 'purple', ylim = c(0,yuplim), add=TRUE)
    
    if (modelpos < length(modelinds)) {
      axis(side = 1, pos = 0, tck = 0, labels = F)
    } else {
      axis(side = 1, pos = 0)
    }
  }
}
#save plot
title(xlab = mvaraxlab,
      ylab = 'Density',
      outer = TRUE, line = 2.8, cex.lab = 1.5)
dev.off()

#export equivalent percentiles information
#write.csv(equivpercent, file = paste0(resultsdir, '/CMIP6multimodel_', mvartype, '-equivalentpercent-JJA-allmodels.csv'), row.names = modelnames)
