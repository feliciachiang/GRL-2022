# Create scatterplot of SNR & non-exceedance values for VPD 
#
# Import: 1) SNR for VPD
#         2) Non-exceedance values for VPD conditioned on pr, conditioned on mrsos
# Output: Scatterplots of SNR & non-exceedance values

#------------------------------------------------------------------------------

#import libraries
library(tidyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(reshape2)

#working directory
wdir = "/Users/fchiang/GISS/"
scriptdir = paste0(wdir, "HSM 2022 scripts/")
resultsdir = paste0(wdir, "HSM project/results/")
figdir = paste0(wdir, "HSM project/figures")

vartype = 'vpd'

#import snr
snr = read.csv(paste0(resultsdir, '/CMIP6multimodel_', vartype, '-snr-JJA-allmodels.csv'))
colnames(snr) <- c('Model', 'SNR')

#import non-exceedance values for VPD
#conditioned on pr
mvartype = 'mrsos'
nonexceeds = read.csv(paste0(resultsdir, '/CMIP6multimodel_', mvartype, vartype, '-nonexceeding-baseline75th-JJA-allmodels_wbaseline_75th.csv'))
nonexceeds = nonexceeds*100

highexceedmodel = c(0,1,1,1,1,0,0,0,0,0,0,0)

df <- data.frame(snr, nonexceeds$V3, highexceedmodel)
colnames(df) <- c('Model', 'SNR', 'NE', 'High')

#remove MIROC-ES2L 2nd degree of warming non-exceedance value for mrsos
df <- df[-c(9),]

jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-snrxne-JJA-allmodels-scatter.jpeg'), width = 6, height = 4, units = 'in', res = 300)
g <- ggplot(df, aes(x = SNR, y = NE)) + geom_point(aes(color = High), size = 2) + geom_text_repel(label = df$Model, hjust = 0, nudge_x = 0.13, segment.color = NA) + xlim(0,4.8) + ylab('Non-exceedance percent') + xlab('Signal to noise ratio')
print(g)
dev.off()

jpeg(paste0(figdir, '/CMIP6multimodel_', mvartype, vartype, '-snrxne-JJA-allmodels-scatter_regline.jpeg'), width = 6, height = 4, units = 'in', res = 300)
g <- ggplot(df, aes(x = SNR, y = NE)) + geom_smooth(method = 'lm', se = FALSE) + geom_point(size = 2) + geom_text_repel(label = df$Model, hjust = 0, nudge_x = 0.13) + xlim(0,4.8) + ylab('Non-exceedance percent') + xlab('Signal to noise ratio')
print(g)
dev.off()