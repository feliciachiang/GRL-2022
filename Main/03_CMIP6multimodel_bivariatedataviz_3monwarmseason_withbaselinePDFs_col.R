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
endmonthinds <- 8

#set percentile to sample 
percentile <- 0.75

#set colors for regression lines for each degree of warming
colvector <- c('orange', 'red', 'purple')

#------------------------------------------------------------------------------#

#import moisture and heat data: pr and tmax------------------------------------
mvartype <- 'pr'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
mvaraxlabnounit <- 'Precipitation'
maggtype <- 'sum'
# 
vartype <- 'tasmax'
varstring <- 'temperature'
varaxlab <- 'Maximum temperature (K)'
aggtype <- 'mean'

#import moisture and heat data: mrsos and tmax---------------------------------
# mvartype <- 'mrsos'
# mvarstring <- 'Mrsos'
# mvaraxlab <- 'Surface soil moisture (kg m-2)'
# mvaraxlabnounit <- 'Surface Soil Moisture'
# maggtype <- 'sum'
# 
# vartype <- 'tasmax'
# varstring <- 'temperature'
# varaxlab <- 'Maximum temperature (K)'
# aggtype <- 'mean'

#import moisture and heat data: pr and vpd------------------------------------
# mvartype <- 'pr'
# mvarstring <- 'Pr'
# mvaraxlab <- 'Precipitation (kg m-2 s-1)'
# mvaraxlabnounit <- 'Precipitation'
# maggtype <- 'sum'
# 
# vartype <- 'vpd'
# varstring <- 'VPD'
# varaxlab <- 'Vapor pressure deficit (mb)'
# aggtype <- 'mean'

#import moisture and heat data: mrsos and vpd---------------------------------
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface Soil Moisture (kg m-2)'
mvaraxlabnounit <- 'Surface Soil Moisture'
maggtype <- 'sum'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#import moisture and heat data: pr and tas------------------------------------
# mvartype <- 'pr'
# mvarstring <- 'Pr'
# mvaraxlab <- 'Precipitation (kg m-2 s-1)'
# maggtype <- 'sum'
# 
# vartype <- 'tas'
# varstring <- 'temperature'
# varaxlab <- 'Temperature (K)'
# aggtype <- 'mean'

#import moisture and heat data: mrsos and tas---------------------------------
# mvartype <- 'mrsos'
# mvarstring <- 'Mrsos'
# mvaraxlab <- 'Surface soil moisture (kg m-2)'
# maggtype <- 'sum'
# 
# vartype <- 'tas'
# varstring <- 'temperature'
# varaxlab <- 'Temperature (K)'
# aggtype <- 'mean'

#precipitation x additional variables----------------------------------------
mvartype <- 'pr'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
maggtype <- 'sum'

#latent heat flux (hfls), Bowen ratio (bo), relative humidity (hurs), specific humidity (huss), actual vapor pressure (ea)
vartype <- 'hfls'
varaxlab <- 'Latent heat flux (W m-2)'
aggtype <- 'mean'
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

#soil moisture x additional variables----------------------------------------
# mvartype <- 'mrsos'
# mvarstring <- 'Mrsos'
# mvaraxlab <- 'Surface soil moisture (kg m-2)'
# maggtype <- 'sum'
# 
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

#import equivalent percentile information for the moisture variable
equivpercent <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, '-equivalentpercent-JJA-allmodels.csv'))


#save copula family, goodness of fit pvalue, and copula parameters
#for each pair of variables
#save for each model and summer month period 
copfamily <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4), c(length(endmonthinds), length(modelinds), 4))
names(dim(copfamily)) <- c('endmonth', 'model', 'warmingdegree')
gofpval <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4), c(length(endmonthinds), length(modelinds), 4))
names(dim(gofpval)) <- c('endmonth', 'model', 'warmingdegree')
copparams <- array(rep(NaN, length(modelinds)*length(endmonthinds)*4*3), c(length(endmonthinds), length(modelinds), 4, 3))
names(dim(copparams)) <- c('endmonth', 'model', 'warmingdegree', 'params')

#create grid to plot on:
#jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-conditionalPDFs_withbaseline_anom-allmonths.jpeg'), width = 8*2, height = 4*2, units = 'in', res = 300)
#par(mfrow = c(length(modelinds), length(endmonthinds)))

jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-conditionalPDFs_withbaseline_anom-JJA-allmodels_col.jpeg'), width = 3.5/1.4, height = 12/1.4, units = 'in', res = 300)
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
    
    #baseline data
    cor(data_baseline, method = 'kendall')[2]
    
    #use function to determine copula family, parameters, and goodness of fit pval
    copfitlist <- findbestfitcopula(data_x = data_baseline[,1], data_y = data_baseline[,2], familynums = c(1,3:6), selectionmethod = "BIC")
    print(copfitlist)
    
    #save copula family name
    copfamily[endmonthpos, modelpos, 1] <- copfitlist[[1]]$familyname
    
    #save gof pval
    gofpval[endmonthpos, modelpos, 1] <- copfitlist[[2]]
    
    #save copula params
    copparams[endmonthpos, modelpos, 1, 1] <- copfitlist[[1]]$family
    copparams[endmonthpos, modelpos, 1, 2] <- copfitlist[[1]]$par
    copparams[endmonthpos, modelpos, 1, 3] <- copfitlist[[1]]$par2
    
    fittedcopulabaseline <- copulaFromFamilyIndex(copfitlist[[1]]$family, copfitlist[[1]]$par, copfitlist[[1]]$par2)
    #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
    cond_dens_empbase <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopulabaseline, data_baseline[,1], data_baseline[,2]))
    #cond_dens_parbase <- make_par_cpdf(baseline75p, fittedcopulabaseline, data_baseline[,1], data_baseline[,2])
    #plot(cond_dens_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim), col = 'black')
    
    #create empirical estimation of marginal distribution of temperature
    #make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,1]))
    #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
    
    
    
    #1 degree of warming
    
    #look at the joint behavior - measure kendall's tau for each degree
    cor(data_1deg, method = 'kendall')[2]
    
    #use function to determine copula family, parameters, and goodness of fit pval
    copfitlist <- findbestfitcopula(data_x = data_1deg[,1], data_y = data_1deg[,2], familynums = c(1,3:6), selectionmethod = "BIC")
    print(copfitlist)
    
    #save copula family name
    copfamily[endmonthpos, modelpos, 2] <- copfitlist[[1]]$familyname
    
    #save gof pval
    gofpval[endmonthpos, modelpos, 2] <- copfitlist[[2]]
    
    #save copula params
    copparams[endmonthpos, modelpos, 2, 1] <- copfitlist[[1]]$family
    copparams[endmonthpos, modelpos, 2, 2] <- copfitlist[[1]]$par
    copparams[endmonthpos, modelpos, 2, 3] <- copfitlist[[1]]$par2
    
    #OR:
    # var_u <- pobs(data_1deg)[,1] #tasmax
    # var_v <- pobs(data_1deg)[,2] #pr
    # #choose from Gaussian, Clayton, Gumbel, Frank, Joe, and associated rotated copulas
    # selectedCopula <- BiCopSelect(var_u, var_v, familyset = c(1,3:6), selectioncrit = "BIC")
    # print(selectedCopula)
    # 
    # #pvalue should be greater than 0.05 to be considered a good fit
    # pval <- BiCopGofTest(var_u, var_v, selectedCopula$family, selectedCopula$par, selectedCopula$par2)$p.value
    # print(pval)
    
    fittedcopula1 <- copulaFromFamilyIndex(copfitlist[[1]]$family, copfitlist[[1]]$par, copfitlist[[1]]$par2)
    #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
    cond_dens_emp1 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula1, data_1deg[,1], data_1deg[,2]))
    #plot(cond_dens_emp1, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))
    
    #create empirical estimation of marginal distribution of temperature
    #make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
    
    #2 degrees of warming
    
    #look at the joint behavior - measure kendall's tau for each degree
    cor(data_2deg, method = 'kendall')[2]
    
    #use function to determine copula family, parameters, and goodness of fit pval
    copfitlist <- findbestfitcopula(data_x = data_2deg[,1], data_y = data_2deg[,2], familynums = c(1,3:6), selectionmethod = "BIC")
    print(copfitlist)
    
    #save copula family name
    copfamily[endmonthpos, modelpos, 3] <- copfitlist[[1]]$familyname
    
    #save gof pval
    gofpval[endmonthpos, modelpos, 3] <- copfitlist[[2]]
    
    #save copula params
    copparams[endmonthpos, modelpos, 3, 1] <- copfitlist[[1]]$family
    copparams[endmonthpos, modelpos, 3, 2] <- copfitlist[[1]]$par
    copparams[endmonthpos, modelpos, 3, 3] <- copfitlist[[1]]$par2
    
    fittedcopula2 <- copulaFromFamilyIndex(copfitlist[[1]]$family, copfitlist[[1]]$par, copfitlist[[1]]$par2)
    #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
    cond_dens_emp2 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula2, data_2deg[,1], data_2deg[,2]))
    #plot(cond_dens_emp2, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))
    
    
    #create empirical estimation of marginal distribution of temperature
    #make_marg_emp2 <-  function(x) demp(x, as.vector(data_2deg[,1]))
    
    #3 degrees of warming
    
    #look at the joint behavior - measure kendall's tau for each degree
    cor(data_3deg, method = 'kendall')[2]
    
    #use function to determine copula family, parameters, and goodness of fit pval
    copfitlist <- findbestfitcopula(data_x = data_3deg[,1], data_y = data_3deg[,2], familynums = c(1,3:6), selectionmethod = "BIC")
    print(copfitlist)
    
    #save copula family name
    copfamily[endmonthpos, modelpos, 4] <- copfitlist[[1]]$familyname
    
    #save gof pval
    gofpval[endmonthpos, modelpos, 4] <- copfitlist[[2]]
    
    #save copula params
    copparams[endmonthpos, modelpos, 4, 1] <- copfitlist[[1]]$family
    copparams[endmonthpos, modelpos, 4, 2] <- copfitlist[[1]]$par
    copparams[endmonthpos, modelpos, 4, 3] <- copfitlist[[1]]$par2
    
    fittedcopula3 <- copulaFromFamilyIndex(copfitlist[[1]]$family, copfitlist[[1]]$par, copfitlist[[1]]$par2)
    #use imported make_emp_cpdf function to estimate conditional PDF at the given percentile
    cond_dens_emp3 <- suppressWarnings(make_emp_cpdf_smoothed(baseline75p, fittedcopula3, data_3deg[,1], data_3deg[,2]))
    #plot(cond_dens_emp3, from=min(data_baseline[,1])-(max(data_4deg[,1])-min(data_baseline[,1]))*.1, to=max(data_4deg[,1])+(max(data_4deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim))
    
    
    #create empirical estimation of marginal distribution of temperature
    #make_marg_emp3 <-  function(x) demp(x, as.vector(data_3deg[,1]))
    
    #--------------------------------PLOTS!----------------------------------------#
    
    #establish x lim (from marginals)
    xlowlim <- -4   #for tasmax
    xhilim <- 10     #for tasmax
    
    #xlowlim <- -10   #for VPD
    #xhilim <- 20     #for VPD    
    
    #establish y lim (not standard, not from marginals)
    cond_base_ymax <- max(cond_dens_empbase(seq(xlowlim, xhilim, length.out = 1000)))
    if(equivpercent$V1[modelpos] < 0.99) {
      cond_emp1_ymax <- max(cond_dens_emp1(seq(xlowlim, xhilim, length.out = 1000)))
    } else {
      cond_emp1_ymax <- 0
    }
    if (equivpercent$V2[modelpos] < 0.99) {
      cond_emp2_ymax <- max(cond_dens_emp2(seq(xlowlim, xhilim, length.out = 1000)))
    } else {
      cond_emp2_ymax <- 0
    }
    if (equivpercent$V3[modelpos] < 0.99) {
      cond_emp3_ymax <- max(cond_dens_emp3(seq(xlowlim, xhilim, length.out = 1000)))
    } else  {
      cond_emp3_ymax <- 0
    }
    
    ylowlim <- 0
    yhilim <- max(cond_base_ymax, cond_emp1_ymax, cond_emp2_ymax, cond_emp3_ymax)
    yhilim <- yhilim*1.15
    
    #plot conditional PDFs for 1, 2, 3, and 4 degrees of warming
    #jpeg(paste0(figdir, '/', mvartype, vartype, '-', savestr, '-conditionalPDFs_withbaseline_anom.jpeg'), width = 9, height = 8, units = 'in', res = 300)
    #plot(cond_dens_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim), xlab = paste0('Conditional ', varstring, ' distribution'), ylab = 'Density', main = paste0(currentmodel, ' ', savestr), col = 'black')
    plot(cond_dens_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', axes = FALSE)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7", border = NA)
    plot(cond_dens_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', add = TRUE)
    title(main = paste0(currentmodel), line = 0.5)
    if (equivpercent$V1[modelpos] < 0.99) {
      plot(cond_dens_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
    }
    if (equivpercent$V2[modelpos] < 0.99) {
      plot(cond_dens_emp2, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'red', add=TRUE)
    }
    if (equivpercent$V3[modelpos] < 0.99) {
      plot(cond_dens_emp3, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'purple', add=TRUE)
    }
    #add vline for var baseline threshold
    abline(v=varbaselinethreshold, col = 'black', lty=2)
    #add legend
    #legend('topright', inset = 0.025, legend = c('1850-1899 baseline', paste0(toTitleCase(vartype), ' at 75th percentile of baseline'), '1 Degree', '2 Degrees', '3 Degrees', '4 Degrees'), col = c('black', 'black', 'green', 'orange', 'red', 'purple'), lty = c(1,2,1,1,1,1), cex = 1, box.lty = 1)
    #add title
    #title(paste0(savestr, ' PDFs given ', mvarstring, ' at 75th percentile of baseline'))
    #title('Conditional PDF of Tmax given Pr at the 75th percentile in 1850-1899')
    #export plot 
    #dev.off()
    
    if (modelpos < length(modelinds)) {
      axis(side = 1, pos = 0, tck = 0, labels = F)
    } else {
      axis(side = 1, pos = 0)
    }
    
  }
}
#save plot
#title(xlab = paste0('Conditional ', varstring, ' distribution (K)'),
title(xlab = paste0(varaxlab, ' conditioned on 75th percentile baseline ', mvaraxlabnounit),
      ylab = 'Density',
      outer = TRUE, line = 2.8, cex.lab = 1.5)
dev.off()

#save copfamily
# 
# filename=paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-copfamily-allmonthsalldegs.nc')
# 
# nx <- length(endmonthinds)
# ny <- length(modelinds)
# endmonth_def <- ncdim_def("endingmonth", "monthind", endmonthinds)
# model_def <- ncdim_def("model", "modelind", modelinds)
# 
# warmingdegree_def <- ncdim_def("warmingdegree","warmingdegreeind", 0:3, unlim=TRUE)
# mv <- -999 #missing value to use
# var_temp <- ncvar_def("copulafams", "familynames", list(endmonth_def, model_def, warmingdegree_def), longname="Copula_family_names", mv) 
# 
# ncnew <- nc_create(filename, list(var_temp))
# 
# print(paste("The file has", ncnew$nvars,"variables"))
# #[1] "The file has 1 variables"
# print(paste("The file has", ncnew$ndim,"dimensions"))
# #[1] "The file has 3 dimensions"
# 
# ncvar_put(ncnew, var_temp, copfamily, start=c(1,1,1), count=c(nx,ny,4), verbose = TRUE)
# 
# # Don't forget to close the file
# nc_close(ncnew)
# 
# #read in to check netcdf files
# #import moisture variable data
# nc <- nc_open(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-copfamily-allmonthsalldegs.nc'))
# data <- ncvar_get(nc, nc$var[["copulafams"]])  
# nc_close(nc)

#save gofpval

filename=paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-gofpval-JJA-allmodels.nc')

nx <- length(endmonthinds)
ny <- length(modelinds)
endmonth_def <- ncdim_def("endingmonth", "monthind", endmonthinds)
model_def <- ncdim_def("model", "modelind", modelinds)

warmingdegree_def <- ncdim_def("warmingdegree","warmingdegreeind", 0:3, unlim=TRUE)
mv <- -999 #missing value to use
var_temp <- ncvar_def("gofpvals", "pvals", list(endmonth_def, model_def, warmingdegree_def), longname="Goodness-of-fit-pvalues", mv) 

ncnew <- nc_create(filename, list(var_temp))

print(paste("The file has", ncnew$nvars,"variables"))
#[1] "The file has 1 variables"
print(paste("The file has", ncnew$ndim,"dimensions"))
#[1] "The file has 3 dimensions"

ncvar_put(ncnew, var_temp, gofpval, start=c(1,1,1), count=c(nx,ny,4))

# Don't forget to close the file
nc_close(ncnew)

#read in to check netcdf files
#import moisture variable data
nc <- nc_open(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-gofpval-JJA-allmodels.nc'))
data <- ncvar_get(nc, nc$var[["gofpvals"]])  
nc_close(nc)

#save copparams

filename=paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-copparams-JJA-allmodels.nc')

nx <- length(endmonthinds)
ny <- length(modelinds)
endmonth_def <- ncdim_def("endingmonth", "monthind", endmonthinds)
model_def <- ncdim_def("model", "modelind", modelinds)

warmingdegree_def <- ncdim_def("warmingdegree","warmingdegreeind", 0:3)

params_def <- ncdim_def("paramters", "parameterind", 1:3, unlim=TRUE)

mv <- -999 #missing value to use
var_temp <- ncvar_def("copparams", "parameters", list(endmonth_def, model_def, warmingdegree_def, params_def), longname="Copula_parameters", mv) 

ncnew <- nc_create(filename, list(var_temp))

print(paste("The file has", ncnew$nvars,"variables"))
#[1] "The file has 1 variables"
print(paste("The file has", ncnew$ndim,"dimensions"))
#[1] "The file has 4 dimensions"

ncvar_put(ncnew, var_temp, copparams, start=c(1,1,1,1), count=c(nx,ny,4, 3))

# Don't forget to close the file
nc_close(ncnew)

#read in to check netcdf files
#import moisture variable data
nc <- nc_open(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-copparams-JJA-allmodels.nc'))
data <- ncvar_get(nc, nc$var[["copparams"]])  
nc_close(nc)
