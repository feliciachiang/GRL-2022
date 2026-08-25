#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:42:39 2022

1) Import EM-Earth North American regional data (#move data from subdirs to maindir)
2) Check  that all data 1950-1969 and 2000-2019 has been downloaded
3) Output monthly 1950-1969 and 2000-2019 netcdf files

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import shutil
import os
import numpy as np
#import fnmatch
#import re
import xarray as xr
#from natsort import natsorted
#import matplotlib.pyplot as plt
import datetime
import pandas as pd

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/EM-Earth/"

# #for each variable
# #var = 'prcp'
# #var = 'tdew'
# #var = 'tmean'
# var = 'trange'

# #move files from subdirectories to main directories
# maindir = datadir + var

# #list subdirectories in main variable folder
# subdirs = os.listdir(maindir)

# for subdir in subdirs: 
#     print(subdir)
    
#     #list files in subdir
#     subdirfiles = os.listdir(maindir + '/' + subdir)
    
#     for subdirfile in subdirfiles:
#         currentfilename = maindir + '/' + subdir + '/' + subdirfile
        
#         shutil.move(currentfilename, maindir)

#import CNA region 

#check that all data has been downloaded by ensemble member-------------------
var = 'prcp'
#var = 'tdew'
#var = 'tmean'
#var = 'trange'

varlist = ['tdew', 'tmean', 'trange']

for var in varlist:
    print(var)

    currentdir = datadir + var
    os.listdir(currentdir)
    os.chdir(currentdir)
    
    #import all data from 1950-1969
    
    # initializing date
    test_date = datetime.datetime.strptime("195001", "%Y%m")
     
    # initializing K
    K = 12*20
     
    date_generated = pd.date_range(test_date, periods=K, freq='M').strftime('%Y%m')
    
    ensmemberlist = ['001', '002', '003', '004', '005']
    
    
    #for each ensemble member
    for ensmember in ensmemberlist:
        #create list
        filelist = ['EM_Earth_probabilistic_daily_' + var + '_NorthAmerica_' + currentdate + '_' + ensmember + '.nc' for currentdate in date_generated]
     
        #try to import- 7305 days for each 20 year period
        ensmemberdata = xr.open_mfdataset(filelist, combine = 'by_coords')
        
        #convert into monthly data and export specific grid cells
        monthly_means = ensmemberdata.resample(time="M").mean()
        
        #export 
        monthly_means[var].to_netcdf(path = datadir + 'monthly/' + 'EM_Earth_probabilistic_monthly_' + var + '_NorthAmerica_195001-196912_' + ensmember + '.nc')
    
    #import all data from 2000-2019
    
    # initializing date
    test_date = datetime.datetime.strptime("200001", "%Y%m")
     
    # initializing K
    K = 12*20
     
    date_generated = pd.date_range(test_date, periods=K, freq='M').strftime('%Y%m')
    
    ensmemberlist = ['001', '002', '003', '004', '005']
    
    
    #for each ensemble member
    for ensmember in ensmemberlist:
        #create list
        filelist = ['EM_Earth_probabilistic_daily_' + var + '_NorthAmerica_' + currentdate + '_' + ensmember + '.nc' for currentdate in date_generated]
     
        #try to import- 7305 days for each 20 year period
        ensmemberdata = xr.open_mfdataset(filelist, combine = 'by_coords')
        
        #convert into monthly data and export specific grid cells
        monthly_means = ensmemberdata.resample(time="M").mean()
        
        #export 
        monthly_means[var].to_netcdf(path = datadir + 'monthly/' + 'EM_Earth_probabilistic_monthly_' + var + '_NorthAmerica_200001-201912_' + ensmember + '.nc')

