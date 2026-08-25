# EM-Earth analysis 
#
# Import: 3-month warm season (June-Aug) data corresponding to 1950-1969 and 2000-2019 years from EM-Earth data
# Output: Scatterplots (pr, tmax), (pr, VPD); conditional PDFs for 1950-1969 and 2000-2019

# Written by: Felicia Chiang, felicia.chiang@nasa.gov

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
source(paste0(scriptdir, "functions/", "t_test_compareslopes.R"))

monstrlist <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

#--------------------------------specify---------------------------------------#

#3-month aggreegate
sc <- 3

#region name 
regionname <- 'CNA'

#summer months 
endmonthinds <- 8

#set percentile to sample 
percentile <- 0.75

#set colors for regression lines for each degree of warming
colvector <- c('orange')

#------------------------------------------------------------------------------#

#import moisture and heat data: pr and tmax------------------------------------
mvartype <- 'prcp'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (mm day-1)'
mvaraxlabnounit <- 'Precipitation'
maggtype <- 'sum'

vartype <- 'tmax'
varstring <- 'temperature'
varaxlab <- 'Maximum temperature (K)'
aggtype <- 'mean'

#import moisture and heat data: pr and vpd------------------------------------
mvartype <- 'prcp'
mvarstring <- 'Pr'
mvaraxlab <- 'Precipitation (kg m-2 s-1)'
mvaraxlabnounit <- 'Precipitation'
maggtype <- 'sum'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#import equivalent percentile information for the moisture variable
#equivpercent <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, '-equivalentpercent-JJA-allmodels.csv'))


#save copula family, goodness of fit pvalue, and copula parameters
#for each pair of variables
#save for each model and summer month period 
copfamily <- array(rep(NaN, length(endmonthinds)*2), c(length(endmonthinds), 2))
names(dim(copfamily)) <- c('endmonth', 'warmingdegree')
gofpval <- array(rep(NaN, length(endmonthinds)*2), c(length(endmonthinds), 2))
names(dim(gofpval)) <- c('endmonth', 'warmingdegree')
copparams <- array(rep(NaN, length(endmonthinds)*2*3), c(length(endmonthinds), 2, 3))
names(dim(copparams)) <- c('endmonth', 'warmingdegree', 'params')

#create grid to plot on:
#jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-conditionalPDFs_withbaseline_anom-allmonths.jpeg'), width = 8*2, height = 4*2, units = 'in', res = 300)
#par(mfrow = c(length(modelinds), length(endmonthinds)))

#jpeg(paste0(figdir, '/EMEarth_', mvartype, vartype, '-conditionalPDFs_anom-JJA.jpeg'), width = 8, height = 4, units = 'in', res = 300)
#par(mfcol = c(12, 1), mar = numeric(4), oma = c(4,4,.5,.5), mgp = c(2, .6, 0))

#specify season
endmonth <- endmonthinds[1]
savestr <- paste0(monstrlist[endmonth-sc+1], '-', monstrlist[endmonth])

print(savestr)

#----------------------------import data------------------------------------

#import moisture variable data
nc <- nc_open(paste0(resultsdir, 'EM_Earth_probabilistic_monthly_CNAregion_', mvartype, '_', maggtype, '_3mon_', monstrlist[endmonth], 'end_aggregated195001-196912timeseries_alt.nc'))
mvarp1data <- ncvar_get(nc, nc$var[[mvartype]])  
nc_close(nc)

nc <- nc_open(paste0(resultsdir, 'EM_Earth_probabilistic_monthly_CNAregion_', mvartype, '_', maggtype, '_3mon_', monstrlist[endmonth], 'end_aggregated200001-201912timeseries_alt.nc'))
mvarp2data <- ncvar_get(nc, nc$var[[mvartype]])  
nc_close(nc)

#import heat variable data
nc <- nc_open(paste0(resultsdir, 'EM_Earth_probabilistic_monthly_CNAregion_', vartype, '_', aggtype, '_3mon_', monstrlist[endmonth], 'end_aggregated195001-196912timeseries_alt.nc'))
varp1data <- ncvar_get(nc, nc$var[[vartype]])  
nc_close(nc)

nc <- nc_open(paste0(resultsdir, 'EM_Earth_probabilistic_monthly_CNAregion_', vartype, '_', aggtype, '_3mon_', monstrlist[endmonth], 'end_aggregated200001-201912timeseries_alt.nc'))
varp2data <- ncvar_get(nc, nc$var[[vartype]])  
nc_close(nc)


#------------------------------------------------------------------------------
#find seasonal averages from the baseline (1950-1969)
baseline_seasonalmean <- mean(varp1data)

#anomalize heat data according to seasonal average from the baseline from all ensemble members
varp2data <- varp2data - baseline_seasonalmean

varp1data <- varp1data-baseline_seasonalmean

#find 75th percentile line from the anomalized vardata (temperature or VPD) baseline time period of 1950-1969
varbaselinethreshold <- quantile(varp1data, percentile)

#check that 75th percentile from 1950-1969 is still within the range of the 2000-2019 data

#use negative precipitation in order to have positive correlation for x and y
data_baseline <- data.frame(as.vector(varp1data), as.vector(-mvarp1data))
data_1deg <- data.frame(as.vector(varp2data), as.vector(-mvarp2data))

#import baseline 1850-1899 data for the relevant season
baselinedata <- mvarp1data

#75th percentile based on data from 1850-1899 (also negative)
baseline75p <- quantile(c(baselinedata), percentile, na.rm = TRUE)
baseline75p <- -baseline75p

#baseline data
cor(data_baseline, method = 'kendall')[2]

# #plot 1950-1969 against 2000-2019 for pair of variables
# xlim = c(min(varp1data) - 1, max(varp2data + 1))
# ylim = c(min(-mvarp1data)*1.05, max(-mvarp2data)*0.95)
# 
# #fit linear models to baseline and degree data
# data_base_lm <- lm(data_baseline[,2]~data_baseline[,1])
# data_deg_lm <- lm(data_1deg[,2]~data_1deg[,1])
# 
# #find slopes for baseline and degree data
# data_base_slope <- data_base_lm$coefficients[[2]]
# #print(data_base_slope)
# data_deg_slope <- data_deg_lm$coefficients[[2]]
# #print(data_deg_slope)
# 
# #find t-test p-value (comparing slopes for baseline and degree data)
# pval <- t_test_compareslopes(data_base, data_deg)
# #print(pval)
# 
# #plot bivariate data and include best fit lines for baseline and for degree data---------------------------------------------
# #print t-test results for slope
# jpeg(paste0(figdir, '/EM-Earth_', mvartype, vartype, '-', savestr, '-scatter.jpeg'), width = 4.5*1.25, height = 4*1.25, units = 'in', res = 300)
# 
# #par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
# if (pval < 0.05) {
#   maintext = paste0('2000-2019', '*')
# } else {
#   maintext = paste0('2000-2019')
# }
# 
# plot(data_baseline[,1], data_baseline[,2], xlim = xlim, ylim = ylim, xlab = varaxlab, ylab = mvaraxlab, main = maintext, font.main = 1, col= 'green')
# 
# points(data_1deg[,1], data_1deg[,2], xlim = xlim, ylim = ylim, xlab = varaxlab, ylab = mvaraxlab, main = maintext, font.main = 1, col= 'orange')
# #add Kendall Tau's correlation value
# text(xlim[1], ylim[2]- 0.000002, bquote(tau == .(round(cor(data_deg, method = 'kendall')[2], 3))), col = 'black', pos = 4)
# 
# abline(data_base_lm, col = 'green')
# #text(xlim[1], ylim[2] -0.000009, paste0("y ~ ",round(data_base_lm$coefficients[1],4)," + ", round(data_base_lm$coefficients[2],8),"x"), col = 'black', pos = 4)
# abline(data_deg_lm, col = 'orange')
# #text(xlim[1], ylim[2] -0.000016, paste0("y ~ ",round(data_deg_lm$coefficients[1],4)," + ", round(data_deg_lm$coefficients[2],8),"x"), col = 'black', pos = 4)  
# 
# #save plot
# dev.off()
# 
# #plot marginals for current variable (tmax, VPD)------------------------------------------------------------------
# 
# #create empirical estimation of marginal distribution of variable
# #make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,1]))
# #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
# 
# #1 degree of warming
# #make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
# #make_marg_emp1 <- suppressWarnings(kdensity(as.vector(data_1deg[,1])))
# 
# 
# #--------------------------------PLOTS!----------------------------------------#
# 
# #establish x lim (standard across models)
# xlowlim <- -4   #for tasmax
# xhilim <- 4     #for tasmax
# 
# #xlowlim <- -4   #for VPD
# #xhilim <- 5     #for VPD
# 
# #set upper limit for y axis
# #yuplim <- 0.65 #-> tasmax, vpd, tas
# #yuplim <- 0.2 #-> hfls
# #yuplim <- 10 #-> bo
# #yuplim <- 0.25 #-> hurs
# #yuplim <- 1200 #-> huss
# #yuplim <- 1 #-> ea
# 
# #establish y lim (not standard across models)
# #account for highest value in all curves
# marg_base_ymax <- max(make_marg_empbase(seq(xlowlim, xhilim, length.out = 1000)))
# marg_emp1_ymax <- max(make_marg_emp1(seq(xlowlim, xhilim, length.out = 1000)))
# 
# ylowlim <- 0
# yhilim <- max(marg_base_ymax, marg_emp1_ymax)
# yhilim <- yhilim*1.15
# 
# jpeg(paste0(figdir, '/EM-Earth_', vartype, '-', savestr, '-marginalPDF.jpeg'), width = 4.5*1.25, height = 4*1.25, units = 'in', res = 300)
# 
# #plot empirical marginal temperature/VPD/alt distributions for baseline, 1-3 degrees of warming
# plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'Maximum temperature anomaly (K)', ylab = 'Density')
# title(main = 'Maximum temperature distribution', line = 0.5)
# plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
# legend('topright', inset = 0.025, legend = c('1950-1969', '2000-2019'), col = c('green', 'orange'), lty = c(1,1), cex = 1, box.lty = 1)
# dev.off()
# 
# # plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'VPD anomaly (mb)', ylab = 'Density')
# # title(main = 'Vapor pressure deficit distribution', line = 0.5)
# # plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
# # legend('topright', inset = 0.025, legend = c('1950-1969', '2000-2019'), col = c('green', 'orange'), lty = c(1,1), cex = 1, box.lty = 1)
# # dev.off()
# 
# #------------------------------------------------
# 
# #plot marginals for current variable (tmax, VPD) for each ensemble member-------#
# #ens = 5
# 
# #create empirical estimation of marginal distribution of variable
# #make_marg_empbase <-  function(x) demp(x, as.vector(data_baseline[,1]))
# #make_marg_empbase <- suppressWarnings(kdensity(as.vector(data_baseline[,1])))
# #make_marg_empbase <- function(x) demp(x, as.vector(varp1data[,ens]))
# 
# #1 degree of warming
# #make_marg_emp1 <-  function(x) demp(x, as.vector(data_1deg[,1]))
# #make_marg_emp1 <- suppressWarnings(kdensity(as.vector(data_1deg[,1])))
# #make_marg_emp1 <- function(x) demp(x, as.vector(varp2data[,ens]))
# 
# #--------------------------------PLOTS!----------------------------------------#
# 
# #establish x lim (standard across models)
# xlowlim <- -4   #for tasmax
# xhilim <- 4     #for tasmax
# 
# #xlowlim <- -4   #for VPD
# #xhilim <- 5     #for VPD
# 
# #set upper limit for y axis
# #yuplim <- 0.65 #-> tasmax, vpd, tas
# #yuplim <- 0.2 #-> hfls
# #yuplim <- 10 #-> bo
# #yuplim <- 0.25 #-> hurs
# #yuplim <- 1200 #-> huss
# #yuplim <- 1 #-> ea
# 
# #establish y lim (not standard across models)
# #account for highest value in all curves
# marg_base_ymax <- max(make_marg_empbase(seq(xlowlim, xhilim, length.out = 1000)))
# marg_emp1_ymax <- max(make_marg_emp1(seq(xlowlim, xhilim, length.out = 1000)))
# 
# ylowlim <- 0
# yhilim <- max(marg_base_ymax, marg_emp1_ymax)
# yhilim <- yhilim*1.15
# 
# jpeg(paste0(figdir, '/EM-Earth_ensmember', as.character(ens), '_', vartype, '-', savestr, '-marginalPDF.jpeg'), width = 4.5*1.25, height = 4*1.25, units = 'in', res = 300)
# 
# #plot empirical marginal temperature/VPD/alt distributions for baseline, 1-3 degrees of warming
# plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'Maximum temperature anomaly (K)', ylab = 'Density')
# title(main = paste0('Ensemble member ', as.character(ens)), line = 0.5)
# plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
# legend('topright', inset = 0.025, legend = c('1950-1969', '2000-2019'), col = c('green', 'orange'), lty = c(1,1), cex = 1, box.lty = 1)
# dev.off()
# 
# #plot CDF
# jpeg(paste0(figdir, '/EM-Earth_ensmember', as.character(ens), '_', vartype, '-', savestr, '-cdf.jpeg'), width = 4.5*1.25, height = 4*1.25, units = 'in', res = 300)
# plot(ecdf(varp1data[,ens]), col = 'green', xlab = 'Maximum temperature anomaly (K)')
# plot(ecdf(varp2data[,ens]), col = 'orange', add = TRUE)
# title(main = paste0('Ensemble member ', as.character(ens)), line = 0.5)
# 
# dev.off()
# 
# # plot(make_marg_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'VPD anomaly (mb)', ylab = 'Density')
# # title(main = 'Vapor pressure deficit distribution', line = 0.5)
# # plot(make_marg_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)
# # legend('topright', inset = 0.025, legend = c('1950-1969', '2000-2019'), col = c('green', 'orange'), lty = c(1,1), cex = 1, box.lty = 1)
# # dev.off()

#------------------------------------------------

#use function to determine copula family, parameters, and goodness of fit pval
copfitlist <- findbestfitcopula(data_x = data_baseline[,1], data_y = data_baseline[,2], familynums = c(1,3:6), selectionmethod = "BIC")
print(copfitlist)

#save copula family name
copfamily[endmonthpos, 1] <- copfitlist[[1]]$familyname

#save gof pval
gofpval[endmonthpos, 1] <- copfitlist[[2]]

#save copula params
copparams[endmonthpos, 1, 1] <- copfitlist[[1]]$family
copparams[endmonthpos, 1, 2] <- copfitlist[[1]]$par
copparams[endmonthpos, 1, 3] <- copfitlist[[1]]$par2

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
copfamily[endmonthpos, 2] <- copfitlist[[1]]$familyname

#save gof pval
gofpval[endmonthpos, 2] <- copfitlist[[2]]

#save copula params
copparams[endmonthpos, 2, 1] <- copfitlist[[1]]$family
copparams[endmonthpos, 2, 2] <- copfitlist[[1]]$par
copparams[endmonthpos, 2, 3] <- copfitlist[[1]]$par2

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

#--------------------------------PLOTS!----------------------------------------#

#establish x lim (from marginals)
xlowlim <- -3   #for tasmax
xhilim <- 3     #for tasmax

#xlowlim <- -4   #for VPD
#xhilim <- 4     #for VPD    

#establish y lim (not standard, not from marginals)
cond_base_ymax <- max(cond_dens_empbase(seq(xlowlim, xhilim, length.out = 1000)))
cond_emp1_ymax <- max(cond_dens_emp1(seq(xlowlim, xhilim, length.out = 1000)))

ylowlim <- 0
yhilim <- max(cond_base_ymax, cond_emp1_ymax, cond_emp2_ymax, cond_emp3_ymax)
yhilim <- yhilim*1.15

#plot conditional PDFs for 1, 2, 3, and 4 degrees of warming
#jpeg(paste0(figdir, '/', mvartype, vartype, '-', savestr, '-conditionalPDFs_withbaseline_anom.jpeg'), width = 9, height = 8, units = 'in', res = 300)

jpeg(paste0(figdir, '/EMEarth_', mvartype, vartype, '-conditionalPDFs_anom-JJA.jpeg'), width = 4.5*1.2, height = 4*1.2, units = 'in', res = 300)
#plot(cond_dens_empbase, from=min(data_baseline[,1])-(max(data_3deg[,1])-min(data_baseline[,1]))*.1, to=max(data_3deg[,1])+(max(data_3deg[,1])-min(data_baseline[,1]))*.1, ylim = c(0,yuplim), xlab = paste0('Conditional ', varstring, ' distribution'), ylab = 'Density', main = paste0(currentmodel, ' ', savestr), col = 'black')
#plot(cond_dens_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'Maximum temperature (K)', ylab = 'Density')
plot(cond_dens_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', xlab = 'Vapor pressure deficit (mb)', ylab = 'Density')
#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7", border = NA)
plot(cond_dens_empbase, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'green', add = TRUE)
title(paste0(varaxlab, ' conditioned on 75th percentile of ', mvaraxlabnounit, ' (1950-1969)'))
plot(cond_dens_emp1, from=xlowlim, to=xhilim, xlim = c(xlowlim, xhilim), ylim = c(ylowlim, yhilim), col = 'orange', add=TRUE)

#add vline for var baseline threshold
abline(v=varbaselinethreshold, col = 'black', lty=2)
#legend('topleft', inset = 0.025, legend = c(paste0('75th percentile of ', toTitleCase(vartype), ' (1950-1969)'), '1950-1969', '2000-2019'), col = c('black', 'green', 'orange'), lty = c(2,1,1), cex = 1, box.lty = 1)
#legend('topleft', inset = 0.025, legend = c('1950-1969', '2000-2019'), col = c('green', 'orange'), lty = 1, cex = 1, box.lty = 1)

dev.off()
#add legend
#add title
#title(paste0(savestr, ' PDFs given ', mvarstring, ' at 75th percentile of baseline'))
#title('Conditional PDF of Tmax given Pr at the 75th percentile in 1850-1899')
#export plot 
#dev.off()
  
#save plot
#title(xlab = paste0('Conditional ', varstring, ' distribution (K)'),
#dev.off()

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

filename=paste0(resultsdir, '/EM-Earth_', mvartype, vartype, '-gofpval-JJA.nc')

nx <- length(endmonthinds)
endmonth_def <- ncdim_def("endingmonth", "monthind", endmonthinds)

warmingdegree_def <- ncdim_def("warmingdegree","warmingdegreeind", 0:1, unlim=TRUE)
mv <- -999 #missing value to use
var_temp <- ncvar_def("gofpvals", "pvals", list(endmonth_def, warmingdegree_def), longname="Goodness-of-fit-pvalues", mv) 

ncnew <- nc_create(filename, list(var_temp))

print(paste("The file has", ncnew$nvars,"variables"))
#[1] "The file has 1 variables"
print(paste("The file has", ncnew$ndim,"dimensions"))
#[1] "The file has 2 dimensions"

ncvar_put(ncnew, var_temp, gofpval, start=c(1,1), count=c(nx,2))

# Don't forget to close the file
nc_close(ncnew)

#read in to check netcdf files
#import moisture variable data
nc <- nc_open(paste0(resultsdir, '/EM-Earth_', mvartype, vartype, '-gofpval-JJA.nc'))
data <- ncvar_get(nc, nc$var[["gofpvals"]])  
nc_close(nc)

#save copparams

filename=paste0(resultsdir, '/EM-Earth_', mvartype, vartype, '-copparams-JJA.nc')

nx <- length(endmonthinds)
endmonth_def <- ncdim_def("endingmonth", "monthind", endmonthinds)

warmingdegree_def <- ncdim_def("warmingdegree","warmingdegreeind", 0:1)

params_def <- ncdim_def("paramters", "parameterind", 1:3, unlim=TRUE)

mv <- -999 #missing value to use
var_temp <- ncvar_def("copparams", "parameters", list(endmonth_def, warmingdegree_def, params_def), longname="Copula_parameters", mv) 

ncnew <- nc_create(filename, list(var_temp))

print(paste("The file has", ncnew$nvars,"variables"))
#[1] "The file has 1 variables"
print(paste("The file has", ncnew$ndim,"dimensions"))
#[1] "The file has 4 dimensions"

ncvar_put(ncnew, var_temp, copparams, start=c(1,1,1), count=c(nx,2, 3))

# Don't forget to close the file
nc_close(ncnew)

#read in to check netcdf files
#import moisture variable data
nc <- nc_open(paste0(resultsdir, '/EM-Earth_', mvartype, vartype, '-copparams-JJA.nc'))
data <- ncvar_get(nc, nc$var[["copparams"]])  
nc_close(nc)
