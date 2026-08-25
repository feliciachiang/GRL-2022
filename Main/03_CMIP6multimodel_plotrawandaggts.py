#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 2022
Created on Mon Mar 21 2022

1) Plot sample time series showing original time series for IPCC regions
2) Plot sample time series showing aggregated time series

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
#import fnmatch
#import re
import xarray as xr
#from natsort import natsorted
import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"
figdir = wdir + "HSM project/figures/"

datadir = "/Volumes/GISS/CMIP6/"

#months
monstrlist = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

#3-month aggregate
sc = 3

#region
regionname = 'CNA'

#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#focus on models with 10 ensemble members
modelinds = np.array([1,4,8,9,10])

#for modelnum in range(0, len(modelnames)):
for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  
    varlist = list(ts.keys())
    
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
        
        plt.figure(dpi = 300)
        for memnum in range(0, 10):
            ts[currentvariable][memnum,:].plot(alpha = 0.4)
        plt.xlabel("Year")
        plt.ylabel(currentvariable)
        plt.title(currentmodel)
        plt.tight_layout()
        plt.savefig(figdir + currentmodel + '_' + currentvariable + "_rawtimeseries.jpg")

        
        for endmonth in [5,7,9]:
            
            savestr = monstrlist[endmonth-sc+1] + '-' + monstrlist[endmonth]  
                    
            #import aggregated regional time series
            aggts = xr.open_dataset(resultsdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + str(sc) + 'mon_' + monstrlist[endmonth] + 'end_aggregatedtimeseries.nc')               
            
            plt.figure(dpi = 300)
            for memnum in range(0, 10):
                aggts[currentvariable][memnum,:].plot(alpha = 0.4)
            # plt.axvspan(1990, 2009, facecolor = 'green', alpha = 0.25)
            # plt.axvspan(2013, 2032, facecolor = 'orange', alpha = 0.25)
            # plt.axvspan(2031, 2050, facecolor = 'red', alpha = 0.25)
            # plt.axvspan(2045, 2064, facecolor = 'purple', alpha = 0.25)
            plt.xlabel("Year")
            plt.ylabel(currentvariable)
            plt.title(currentmodel)
            plt.tight_layout()
            plt.savefig(figdir + currentmodel + '_' + currentvariable + '_' + aggtype + '_' + savestr + '_aggtimeseries.jpg')
