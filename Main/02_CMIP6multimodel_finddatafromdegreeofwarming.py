#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 2022

Find 20-year data corresponding to each degree of warming
Input: Seasonal data (01_findseasonaldata.py), start year of 20-year period (00_findstartyear.py)
Output: Netcdf files with 20-year data for 1-4 degrees of warming

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

#focus on models with 10 ensemble members
#modelinds = np.array([1,4,8,9,10])

for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  
    varlist = list(ts.keys())
    #varlist = ['es']
    
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
        
        #use first 10 ensemble members
        members = ds['member'].values
        print(members)
        
        #empty array
        saveds = np.empty((len(members), 50))
        saveds[:] = np.nan
        
        for membernum in range(0, len(members)):
            currentmember = members[membernum]
            print(currentmember)

            saveds[membernum, :] = ds[currentvariable].sel(member = currentmember)[0:50].values
    
            #save as netcdf
            _arrsave = xr.Dataset({currentvariable: (['member', 'years'], saveds)},
                                  coords = {'member': members, 'years': np.arange(1850,1900)}) 
            
            #_arrsave.to_netcdf(resultsdir + 'CanESM5_' + var + '_1850-1899baseline.nc') 
            _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + savestr + '_1850-1899baseline_first5ens.nc')                  
                             
            #np.nanpercentile(_arrsave[currentvariable].values, 75)
            #end in June, so if 3-month aggregate: Apr-June
            #0.00010765772549348185
            
            #end in August, so if 3-month aggregate: June-Aug
            #0.00010294077856087824
            
            #end in October, so if 3-month aggregate: Aug-Oct
            #8.118182540783891e-05


#--------------------------------------------------------------------------

for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  
    varlist = list(ts.keys())
    #varlist = ['es']

    #import start year array
    startyear = xr.open_dataset(resultsdir + currentmodel + '_globalwarmingdegrees.nc')

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
        saveds = np.empty((len(np.arange(1,4)), len(members), 20))
        saveds[:] = np.nan
        
        #find 20-year data for each ensemble member and degree of warming
        for degreenum in np.arange(1,4):
            print(degreenum)    
            for variantnum in range(0, len(members)):
                currentmember = members[variantnum]
                print(currentmember)
                
                #find indices for years closest to degree of warming for this ensemble member 
                startingyear = startyear['startyear'].sel(degree = degreenum).sel(member = currentmember).values.astype(int)
                selectinds = np.arange(startingyear, startingyear+window)
                
                saveds[degreenum-1, variantnum, :] = ds[currentvariable].sel(member = currentmember)[selectinds].values
                
        #save as netcdf
        _arrsave = xr.Dataset({currentvariable: (['degree', 'member', 'years'], saveds)},
                              coords = {'degree': np.arange(1,4), 'member': members, 'years': np.arange(0,20)}) 
        _arrsave.to_netcdf(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + savestr + '_datafromdegreesofwarming_first5ens.nc')                  
                    
