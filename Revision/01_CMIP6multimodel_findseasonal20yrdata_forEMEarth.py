#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 2022

Find 20-year data (1950-1969, 2000-2019) for analysis equivalent to EM-Earth
Input: Seasonal data (01_findseasonaldata.py)
Output: Netcdf files with 20-year data for 1950-1969, 2000-2019 to match EM-Earth analysis

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

#months
monstrlist = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

#length of aggregated data
sc = 3

regionname = 'CNA'

#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  
    varlist = list(ts.keys())
    #varlist = ['es']

    window = 20   
    
    endmonth = 7

    #for endmonth in [5,7,9]:
        
    savestr = monstrlist[endmonth-sc+1] + '-' + monstrlist[endmonth]  
    
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

        ds = xr.open_dataset(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries_first5ens.nc')               
        #ds = xr.open_dataset(resultsdir + 'CanESM5_' + var + '_' + aggtype + '_aggregatedtimeseries.nc')                  
        members = ds['member'].values

        #empty array
        saveds = np.empty((2, len(members), 20))
        saveds[:] = np.nan
        
        #for 1950-1969
        selectinds = np.arange(100,120)
        
        #find 20-year data for each ensemble member
        for variantnum in range(0, len(members)):
            currentmember = members[variantnum]
            print(currentmember)
                        
            saveds[0, variantnum, :] = ds[currentvariable].sel(member = currentmember)[selectinds].values
        
        #for 2000-2019
        selectinds = np.arange(150,170)
        
        #find 20-year data for each ensemble member
        for variantnum in range(0, len(members)):
            currentmember = members[variantnum]
            print(currentmember)
                        
            saveds[1, variantnum, :] = ds[currentvariable].sel(member = currentmember)[selectinds].values
        
        
        #save as netcdf
        _arrsave = xr.Dataset({currentvariable: (['period', 'member', 'years'], saveds)},
                              coords = {'degree': np.arange(1,2), 'member': members, 'years': np.arange(0,20)}) 
        _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + savestr + '_EMEarthperiods_first5ens.nc')                  
                    
