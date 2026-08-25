# Generate nested barplots for quantile associated with 75th percentile from baseline for each degree of warming for each model 
# Generate for tasmax and for VPD
# 
# Input: Jun-Aug quantiles for each degree of warming for each model in csv format
# Output: Nested bar plots (one for tasmax and one for VPD)

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

#set colors for regression lines for each degree of warming
colvector <- c('orange', 'red', 'purple')

#------------------------------------------------------------------------------#

#import variable------------------------------------

# vartype <- 'tasmax'
# varstring <- 'temperature'
# varaxlab <- 'Maximum Temperature (K)'
# aggtype <- 'mean'

vartype <- 'vpd'
varstring <- 'VPD'
varaxlab <- 'Vapor pressure deficit (mb)'
aggtype <- 'mean'

#-----------------------------------------------------------------------------

#read in csv file
quantiles <- read.csv(paste0(resultsdir, '/CMIP6multimodel_', vartype, '-equivalentpercent-JJA-allmodels.csv'))
finquantiles <- quantiles[,2:4]*100

#convert csv table into dataframe
df <- data.frame(finquantiles, row.names = quantiles$X)
colnames(df) <- c('1 Degree', '2 Degrees', '3 Degrees')
df_models <- cbind(Model = rownames(df), df)
rownames(df_models) <- NULL

dat.g <- gather(df_models, type, value, -Model)

#create nested barplot
jpeg(paste0(figdir, '/CMIP6multimodel_',  vartype, '-equivpercent-baseline75th-JJA-allmodels-barplot.jpeg'), width = 8, height = 4, units = 'in', res = 300)
g <- ggplot(dat.g, aes(type, value)) + geom_bar(aes(fill = Model), stat = 'identity', position = 'dodge') + labs(y = 'Equivalent percentile', x = 'Degree of warming') + ylim(0, 100)
print(g)
dev.off()
