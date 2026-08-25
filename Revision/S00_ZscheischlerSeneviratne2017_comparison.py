#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:06:50 2022

Find warm season correlation coefficient for:
    1) Maximum temperature and precipitation
    2) Maximum temperature and soil moisture
    3) VPD and precipitation
    4) VPD and soil moisture
for 1870-1969 and for 2001-2100 time periods.

Find the difference between the warm season correlation coefficients for the two periods.
    
1) Find the warm season (yearly averaged JJA) values for the region for each model.
2) Linearly detrend the warm season time series.
3) Find the original and detrended correlation coefficient (Pearson's, Kendall's) values

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
#import fnmatch
#import re
import xarray as xr
import scipy.stats
#from natsort import natsorted
#import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/CMIP6/"

#import functions to be used
def seasonalmean(xd, sc, endmonth):
    #xd is the time series starting in january
    #sc is the number of months to aggregate
    #endmonth is the last month of the aggregate
    
    A1 = np.empty((len(xd)-sc+1, sc))
    A1[:] = np.nan
    for i in range(0, sc):
        A1[:,i] = xd[i:len(xd)-sc+1+i]

    X = np.nanmean(A1, 1)
        
    #add nan values to time series for sc > 1    
    nanvector = np.empty(sc-1)
    nanvector[:] = np.nan
    X = np.concatenate((nanvector, X))       
    
    #extract appropriate aggregate from monthly time series
    finalts = X[endmonth::12]
    
    #return annual time series 
    return finalts


#months
monstrlist = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

varlist = ['pr', 'mrsos', 'tasmax', 'vpd']

#region
regionname = 'CNA'

#length of aggregated data
sc = 3

#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#focus on models with 10 ensemble members
#modelinds = np.array([1,4,8,9,10])
#modelinds = np.array([0,2,3,5,6,7,11])

#Part 1: Find the warm season values for the regional time series for each model

for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)

    #import regional time series
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  

    members = ts['member'].values
    members = members[0:5]
    
    years = np.unique(ts.time.dt.year)
    varlist = list(ts.keys())
    
    #for early, mid, and late summer periods
    #end in June, so if 3-month aggregate: Apr-June
    #endmonth = 5
    #end in August, so if 3-month aggregate: June-Aug
    endmonth = 7
    #end in October, so if 3-month aggregate: Aug-Oct
    #endmonth = 9
    #for endmonth in [5,7,9]:
    
    #for each variable
    for varnum in range(0, len(varlist)):
        currentvariable = varlist[varnum]
        print(currentvariable)
    
        aggtype = "mean"
        print(aggtype)

        #create empty array to save to
        saveds = np.empty((len(members), len(years)))
        saveds[:] = np.nan
        
        for membernum in range(0, len(members)):
            currentmember = members[membernum]
            print(currentmember)
            
            saveds[membernum,:] = seasonalmean(ts[currentvariable].sel(member = currentmember).values, sc, endmonth)
            
        #save 
        _arrsave = xr.Dataset({currentvariable: (['member', 'year'], saveds)},
                              coords = {'member': members, 'year': years}) 
        
        _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')                  
                 
#Part 2: Linearly detrend the warm season values for each model's regional time series
    #Find the linear regression line
    #Remove the differences from the regression line
    #Save the detrended data as a netcdf file
    
#For each model
for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)

    endmonth = 7

    aggtype = "mean"
    print(aggtype)
    
    #for each variable
    for varnum in range(0, len(varlist)):
        currentvariable = varlist[varnum]
        print(currentvariable)

        #import regional time series
        arr = xr.open_dataset(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')                  
                     
        members = arr['member'].values
        members = members[0:5]
        
        saveds = np.empty((len(members), len(years)))
        saveds[:] = np.nan
        
        #find the regression line for each ensemble member
        for membernum in range(0, len(members)):
            currentmember = members[membernum]
            print(currentmember)
            
            best_x = arr[currentvariable].year.values
            og_y = arr[currentvariable].sel(member = currentmember).values
            
            m, b = np.polyfit(best_x, og_y, 1)
            
            best_y = m*best_x+b
            
            #subtract best_y from original y
            saveds[membernum, :] = og_y-best_y 
            
        _arrsave = xr.Dataset({currentvariable: (['member', 'year'], saveds)},
                              coords = {'member': members, 'year': arr.year.values}) 
        
        _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens_lineardetrended.nc')                  

            
    
#Part 3: Find the original and detrended Pearson's and Kendall's correlation coefficients for 1870-1969 and 2001-2100
#Import the original regional time series and calculate the correlation coefficients
pearsoncorr = np.empty((len(modelnames)*5, 2))
pearsoncorr[:] = np.nan

kendallcorr = np.empty((len(modelnames)*5, 2))
kendallcorr[:] = np.nan

mvar = 'pr'
var = 'tasmax'

countnum = -1

#For each model
for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)

    endmonth = 7

    aggtype = "mean"
    print(aggtype)
    
    #import regional time series
    mvards = xr.open_dataset(resultsdir + currentmodel + '_' + mvar + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')                  
    vards = xr.open_dataset(resultsdir + currentmodel + '_' + var + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')                  

                 
    members = mvards['member'].values
    members = members[0:5]
    
    #find the correlation between precipitation and temperature
    for membernum in range(0, len(members)):
        currentmember = members[membernum]
        print(currentmember)
        
        countnum = countnum + 1
        
        current_x = mvards[mvar].sel(member = currentmember).values
        current_y = vards[var].sel(member = currentmember).values

        #1870-1969
        pearsoncorr[countnum, 0] = scipy.stats.pearsonr(current_x[20:120], current_y[20:120])[0]
        #2001-2100
        pearsoncorr[countnum, 1] = scipy.stats.pearsonr(current_x[151:251], current_y[151:251])[0]
        
        kendallcorr[countnum, 0] = scipy.stats.kendalltau(current_x[20:120], current_y[20:120])[0]
        kendallcorr[countnum, 1] = scipy.stats.kendalltau(current_x[151:251], current_y[151:251])[0]
    
#Import the linearly detrended regional time series and calculate the correlation coefficients
pearsoncorr = np.empty((len(modelnames)*5, 2))
pearsoncorr[:] = np.nan

kendallcorr = np.empty((len(modelnames)*5, 2))
kendallcorr[:] = np.nan

mvar = 'pr'
var = 'tasmax'

countnum = -1

#For each model
for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)

    endmonth = 7

    aggtype = "mean"
    print(aggtype)
    
    #import regional time series
    mvards = xr.open_dataset(resultsdir + currentmodel + '_' + mvar + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens_lineardetrended.nc')                  
    vards = xr.open_dataset(resultsdir + currentmodel + '_' + var + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens_lineardetrended.nc')                  

                 
    members = mvards['member'].values
    members = members[0:5]
    
    #find the correlation between precipitation and temperature
    for membernum in range(0, len(members)):
        currentmember = members[membernum]
        print(currentmember)
        
        countnum = countnum + 1
        
        current_x = mvards[mvar].sel(member = currentmember).values
        current_y = vards[var].sel(member = currentmember).values

        #1870-1969
        pearsoncorr[countnum, 0] = scipy.stats.pearsonr(current_x[20:120], current_y[20:120])[0]
        #2001-2100
        pearsoncorr[countnum, 1] = scipy.stats.pearsonr(current_x[151:251], current_y[151:251])[0]
        
        kendallcorr[countnum, 0] = scipy.stats.kendalltau(current_x[20:120], current_y[20:120])[0]
        kendallcorr[countnum, 1] = scipy.stats.kendalltau(current_x[151:251], current_y[151:251])[0]