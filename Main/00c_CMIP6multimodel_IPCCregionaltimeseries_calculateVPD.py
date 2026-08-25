#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 2022

Calculate VPD from TAS and HURS

Import: 1) Regional time series containing (among others) TAS and HURS variables
        2) Function to calculate VPD (mb)

Output: 1) Regional time series containing SM, PR, TASMAX, and VPD

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
#import fnmatch
#import re
import xarray as xr
#from natsort import natsorted

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/CMIP6/"

#import functions
def calculate_es(tas):
    
    es = 6.11*np.exp((2.5e6/461)*((1/273) - (1/tas)))       
    return es

def calculate_ea(hurs, tas):    
    #calculate saturation vapor pressure
    es = calculate_es(tas)
    
    #calculate actual vapor pressure
    ea = (hurs/100)*es
    
    return ea


def calculate_vpd(hurs, tas):    
    #calculate saturation vapor pressure
    es = calculate_es(tas)
    
    #calculate actual vapor pressure
    ea = (hurs/100)*es
    
    #vapor pressure deficit
    vpd = es - ea
    
    return vpd

#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#region
regionname = 'CNA'

for modelnum in range(0, len(modelnames)):
    currentmodel = modelnames[modelnum]
    print(currentmodel)

    #import regional time series
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries.nc')  
    
    #find actual vapor pressure   
    ts['ea'] = (['member', 'time'], calculate_ea(ts['hurs'], ts['tas']))
    
    #find saturated vapor presure
    ts['es'] = (['member', 'time'], calculate_es(ts['tas']))
    
    #find vpd
    ts['vpd'] = (['member', 'time'], calculate_vpd(ts['hurs'], ts['tas']))   

    #find bowen ratio = sensible heat flux / latent heat flux
    ts['bo'] = (['member', 'time'], ts['hfss']/ts['hfls'])          
    
    #export
    ts.to_netcdf(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')                  

#check raw values for each model and each variable
import matplotlib.pyplot as plt


for modelnum in range(0, len(modelnames)):
    
    # modelnum = -1
    
    # modelnum = modelnum + 1
    
    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    ts = xr.open_dataset(resultsdir + currentmodel + '_' + regionname + 'region_meantimeseries_all.nc')  
    
    for membernum in range(0, 5):

        plt.figure()
        ts['hurs'][membernum,:].plot()
            
        plt.figure()
        ts['mrsos'][membernum,:].plot()
        
        plt.figure()
        ts['pr'][membernum,:].plot()
        
        plt.figure()
        ts['tas'][membernum,:].plot()
        
        plt.figure()
        ts['tasmax'][membernum,:].plot()
        
        plt.figure()
        ts['vpd'][membernum,:].plot()
        
        plt.figure()
        ts['ea'][membernum,:].plot()
        
        plt.figure()
        ts['es'][membernum,:].plot()
    
    
    
    
    
    
