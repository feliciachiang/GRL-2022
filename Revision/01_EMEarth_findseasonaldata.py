#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 2022

Find relevant seasonal (Jun-Aug) data for PR, TMAX, VPD from EM-Earth data

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
#import fnmatch
#import re
import xarray as xr
#from natsort import natsorted
#import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/EM-Earth/"

#timeperiodstr = '195001-196912'
timeperiodstr = '200001-201912'

ts = xr.open_dataset(datadir + 'EM_Earth_probabilistic_monthly_CNAregion_' + timeperiodstr + 'meantimeseries_allalt.nc')                  

#import functions to be used
def aggregate(xd, sc, endmonth, function):
    #xd is the time series starting in january
    #sc is the number of months to aggregate
    #endmonth is the last month of the aggregate
    
    A1 = np.empty((len(xd)-sc+1, sc))
    A1[:] = np.nan
    for i in range(0, sc):
        A1[:,i] = xd[i:len(xd)-sc+1+i]

    if function == "sum":
        X = np.nansum(A1, 1)
    elif function == "mean":
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

#region
regionname = 'CNA'

#length of aggregated data
sc = 3


members = ts['member'].values
members = members[0:5]

years = np.unique(ts.time.dt.year)
varlist = list(ts.keys())

#end in August, so if 3-month aggregate: June-Aug
endmonth = 7

#for each variable
for varnum in range(0, len(varlist)):
    currentvariable = varlist[varnum]
    print(currentvariable)

    #how to aggregate:
    if currentvariable == "prcp":
        aggtype = "sum"
        print(aggtype)
    else:
        aggtype = "mean"
        print(aggtype)

    #create empty array to save to
    saveds = np.empty((len(members), len(years)))
    saveds[:] = np.nan
    
    for membernum in range(0, len(members)):
        currentmember = members[membernum]
        print(currentmember)
        
        saveds[membernum,:] = aggregate(ts[currentvariable].sel(member = currentmember).values, sc, endmonth, aggtype)
        
    #save 
    _arrsave = xr.Dataset({currentvariable: (['member', 'year'], saveds)},
                          coords = {'member': members, 'year': years}) 
    
    _arrsave.to_netcdf(resultsdir + 'EM_Earth_probabilistic_monthly_CNAregion_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregated' + timeperiodstr + 'timeseries_alt.nc')                  
             
