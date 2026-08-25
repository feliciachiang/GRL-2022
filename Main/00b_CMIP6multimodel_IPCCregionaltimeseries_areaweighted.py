#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract average time series data from IPCC AR5 Regions (Used in SREX report amongst others)

Import: 1) IPCC rasters created using IPCC shapefiles
        2) Land area fraction, Grid cell area
        3) Raw CMIP6 datasets

Output: Netcdf files containing average time series data with individual IPCC AR5 regions weighting by grid cell area and land area fraction

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import pandas as pd
import fnmatch
import re
import xarray as xr
from natsort import natsorted
#import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/CMIP6/"

def ESGFinfofromfilelist(filenames, string = "variable"):
    '''Function that returns information from files downloaded from ESGF
    
    string = [variable, frequency, model, forcing, variant, grid, timerange]'''
    
    import numpy as np
    
    #position = 0
    if string == "variable":
        #get the unique variable names from the list of files (filenames)
        names = [name.split("_")[0] for name in filenames]
    elif string == "frequency":
        #get the unique frequency names
        names = [name.split("_")[1] for name in filenames]
    elif string == "model":
        #get the unique model names
        names = [name.split("_")[2] for name in filenames]
    elif string == "forcing":
        #get the unique forcing names
        names = [name.split("_")[3] for name in filenames]
    elif string == "variant":
        #get the unique ensemble member names
        names = [name.split("_")[4] for name in filenames] 
    elif string == "grid":
        #get the unique grid label names
        names = [name.split("_")[5] for name in filenames]         
    elif string == "timerange":    
        #get the time range
        names = [name.split("_")[6] for name in filenames]
    else:
        names = []
        #change to output error: incorrect string
    
    return np.unique(names)


#for each model
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#import resources (such as crop maps)
resourcedir = wdir + "resources/"
#open IPCC regions
regionfilenames = fnmatch.filter(os.listdir(resourcedir + 'referenceRegions/'), '*.nc')
#find region names
regionsearchstr = 'region(.*)_highres.nc';
regionnames = [re.search(regionsearchstr, name).group(1) for name in regionfilenames]

#new longitude and latitude
newlon = np.arange(0.25, 360, 0.5)
newlat = np.arange(-89.75,90, 0.5)

#modelinds = np.array([4,8,9,10])
#modelinds = np.array([2,3,5,6,7,11])

for modelnum in range(0, len(modelnames)):
#for modelnum in modelinds:

    currentmodel = modelnames[modelnum]
    print(currentmodel)

    laffilename = fnmatch.filter(os.listdir(datadir + currentmodel), '*sftlf*')
    areacellafilename = fnmatch.filter(os.listdir(datadir + currentmodel), '*areacella*')
    
    #find land area fraction 
    laf = xr.open_dataset(datadir + currentmodel + "/" + laffilename[0])
    #reindex to match the IPCC region resolution
    laf = laf.reindex(lat = newlat, lon = newlon, method = 'nearest')
    
    #import areacella
    areacella = xr.open_dataset(datadir + currentmodel + "/" + areacellafilename[0])
    areacella = areacella.reindex(lat = newlat, lon = newlon, method = 'nearest')
    
    #find grid cell area normalized by land area fraction
    landareamask = (laf['sftlf']*areacella['areacella']/100)

    #list all files
    modeldir = datadir + currentmodel + '/mon/'
    
    filenames = fnmatch.filter(os.listdir(modeldir), '*.nc')
    filenames = [filename for filename in filenames if '._' not in filename]
    
    #list all unique variables (from GISS model)
    varnames = ESGFinfofromfilelist(filenames, 'variable')
    #list all unique forcings
    forcings = ESGFinfofromfilelist(filenames, 'forcing')
    #list all unique ensemble members
    members = ESGFinfofromfilelist(filenames, 'variant')
    members = natsorted(members)
    
    timelen = 3012
    
    
    #for each IPCC region:
    #for regionnum in range(0, len(regionnames)):
        
    #for CNA region:
    regionnum = 4
    
    print(regionnames[regionnum])
    #open region shapefile
    regionds = xr.open_dataset(resourcedir + 'referenceRegions/' + regionfilenames[regionnum])
    #revise region's latitude and longitude to match the land area fraction and lake fraction reindex lat and lon
    regionds = regionds.sortby('Latitude', ascending = True)
    regionds = regionds.roll(Longitude = 360, roll_coords=True)
    regionds = regionds.assign_coords(Longitude = newlon)
    regionds = regionds.rename({'Latitude': 'lat', 'Longitude': 'lon'})
    
    #plot to check
    #regionds['Region'].where(laf['sftlf']>50).plot()
    
    #find region without ocean and without lake biases
    land = xr.Dataset({'mask': (['lat', 'lon'], regionds['Region'].values)}, 
        coords = {'lat': regionds['lat'].values, 'lon': regionds['lon'].values})
    
    #save region
    #land.to_netcdf(resourcedir + 'AR5regions/' + regionnames[regionnum] + 'region_finalmask.nc')
    
    #create xarray
    _arrsave = xr.Dataset(coords = {'member': members, 'time': pd.date_range("1850-01", periods = timelen, freq = 'M')}) 
    
    #for each ensemble member and variable, output the area-weighted regional average time series
    for varnum in range(0, len(varnames)):
        currentvariable = varnames[varnum]
        print(currentvariable)
                
        #use missingyears function to check if missing years
        varfilenames = [filename for filename in filenames if currentvariable + '_' in filename]

        print(varfilenames)
        
        arr_regionaltimeseries = np.empty((len(members), timelen))
        arr_regionaltimeseries[:] = np.nan

        
        for membernum in range(0, len(members)):
            currentmember = members[membernum]
            print(currentmember)
            
            #currentmember = 'r1i1p1f2'
            
            #use missingyears function to check if missing years
            variantfilenames = fnmatch.filter(varfilenames, '*' + currentmember + '*')
            print(variantfilenames)
            
            filenameswpath = [modeldir + filename for filename in variantfilenames]
            
            #import and combine historical and ssp585 files
            ds = xr.open_mfdataset(natsorted(filenameswpath), combine = 'by_coords', parallel = True, chunks = {'lat': 10, 'lon': 10})
            ds = ds.reindex(lat = newlat, lon = newlon, method = 'nearest')
            
            #calculate area-averaged weights with cos
            #weights = np.cos(np.deg2rad(ds.lat))
            #weights.name = 'weights'
            
            #weight by land area fraction and grid cell area
            weights = landareamask
            weights.name = 'weights'
            
            ds = ds[currentvariable].where(land['mask'] == 1)
    
            #weight by land fraction and grid cell area
            ds_weighted = ds.weighted(weights)
            
            ds_regionalmean = ds_weighted.mean(dim = ('lat', 'lon'), skipna = True)
                        
            #save in array
            arr_regionaltimeseries[membernum, :] = ds_regionalmean
        
        #save in xarray
        _arrsave[currentvariable] = (['member', 'time'], arr_regionaltimeseries)
    
    #export arr_rollingglobalmean
    _arrsave.to_netcdf(resultsdir + currentmodel + '_' + regionnames[regionnum] + 'region_meantimeseries.nc')                  
        
            
            
