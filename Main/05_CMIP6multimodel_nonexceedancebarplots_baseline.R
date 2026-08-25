# Generate nested barplots for non-exceedance probabilities for each degree of warming for each model 
# Non-exceedance based on the 75th percentile of the marginal
# 
# Input: Jun-Aug non-exceedance probabilities for each degree of warming for each model in csv format
# Output: Nested bar plots (one for each pair of variables)

#libraries
library(tidyr)
library(ggplot2)
library(plotly)

#working directory
wdir = "/Users/fchiang/GISS/"
scriptdir = paste0(wdir, "HSM 2022 scripts/")
resultsdir = paste0(wdir, "HSM project/results/")
figdir = paste0(wdir, "HSM project/figures")

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

# #import moisture and heat data: mrsos and vpd---------------------------------
mvartype <- 'mrsos'
mvarstring <- 'Mrsos'
mvaraxlab <- 'Surface soil moisture (kg m-2)'
maggtype <- 'sum'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#-----------------------------------------------------------------------------

#read in csv file
#nonexceeds <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-nonexceeding-baseline75th-JJA-allmodels_wbaseline.csv'))
nonexceeds <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-nonexceeding-baseline75th-JJA-allmodels_wbaseline_', as.character(heatpercentile*100), 'th.csv'))
nonexceeds <- nonexceeds*100
nonexceeds[nonexceeds > 100] = 100

equivpercent <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, '-equivalentpercent-JJA-allmodels.csv'))
equivpercent$V0 <- rep(0.75, 12)
equivpercent <- equivpercent[, c('X', 'V0', 'V1', 'V2', 'V3')]
nonexceeds[equivpercent[,2:5] > 0.99] = NaN

#replace values where conditional PDFs are not drawn with NAN

#convert csv table into dataframe
nonexceedsdf <- data.frame(nonexceeds, row.names = modelnames)
colnames(nonexceedsdf) <- c('Baseline', '1 Degree', '2 Degrees', '3 Degrees')
nonexceedsdf_models <- cbind(Model = rownames(nonexceedsdf), nonexceedsdf)
rownames(nonexceedsdf_models) <- NULL

dat.g <- gather(nonexceedsdf_models, type, value, -Model)
dat.g$type <- factor(dat.g$type, levels = c('Baseline', '1 Degree', '2 Degrees', '3 Degrees'))

#create nested barplot
jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-nonexceeding-baseline75th-JJA-allmodels-barplot-wbaseline_', as.character(heatpercentile*100), 'th.jpeg'), width = 8, height = 4, units = 'in', res = 300)
g <- ggplot(dat.g, aes(type, value)) + geom_bar(aes(fill = Model), stat = 'identity', position = 'dodge') + labs(y = 'Non-exceedance percent', x = 'Degree of warming') + ylim(0, 100)
print(g)
dev.off()
