#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 2022

Find relevant seasonal data for PR, SM (Sept-Mar, Jan-Mar)
Find relevant seasonal data for TMAX, VPD (Apr-Oct)

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

datadir = "/Volumes/GISS/CMIP6/"

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

#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#focus on models with 10 ensemble members
#modelinds = np.array([1,4,8,9,10])
#modelinds = np.array([0,2,3,5,6,7,11])

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
    
        #how to aggregate:
        if currentvariable == "pr" or currentvariable == "mrsos":
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
        
        _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')                  
                 
