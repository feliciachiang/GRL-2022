#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:16:28 2022

1a) Import monthly EM-Earth data 
1b) Import IPCC rasters created using IPCC shapefiles
2) Find area-averaged monthly time series (using latitudes) for Central North America region and export

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import fnmatch
import re
import xarray as xr
from natsort import natsorted
#import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/EM-Earth/"

#import EM-Earth data for lat and lon
ds = xr.open_dataset(datadir + 'monthly/' + 'EM_Earth_probabilistic_monthly_prcp_NorthAmerica_195001-196912_001.nc')
#remove values below -180 longitude
ds_prcp = ds['prcp'][:,:,77:1764]
#sort latitude by ascending values
ds_prcp = ds_prcp.sortby('lat', ascending = True)

#import resources (such as crop maps)
resourcedir = wdir + "resources/"
#open IPCC regions
regionfilenames = fnmatch.filter(os.listdir(resourcedir + 'referenceRegions/'), '*.nc')
regionfilenames = natsorted(regionfilenames)
#find region names
regionsearchstr = 'region(.*)_highres.nc';
regionnames = [re.search(regionsearchstr, name).group(1) for name in regionfilenames]
regionnames = natsorted(regionnames)

#for each IPCC region:
#for regionnum in range(0, len(regionnames)):
    
#for CNA region:
regionnum = 4

print(regionnames[regionnum])
#open region shapefile
regionds = xr.open_dataset(resourcedir + 'referenceRegions/' + regionfilenames[regionnum])
#revise region's latitude and longitude to match the land area fraction and lake fraction reindex lat and lon
regionds = regionds.sortby('Latitude', ascending = True)
regionds = regionds.rename({'Latitude': 'lat', 'Longitude': 'lon'})

#new longitude and latitude
newlon = np.arange(-179.95, 180, 0.1)
newlat = np.arange(-89.95,90, 0.1)

#reindex lat and lon
regionds = regionds.reindex(lat = newlat, lon = newlon, method = 'nearest')

#plot to check
#regionds['Region'].where(laf['sftlf']>50).plot()

#find region without ocean and without lake biases
land = xr.Dataset({'mask': (['lat', 'lon'], regionds['Region'][969:1738,0:1687].values)}, 
    coords = {'lat': ds_prcp['lat'].values, 'lon': ds_prcp['lon'].values})

#save region
#land.to_netcdf(resourcedir + 'AR5regions/' + regionnames[regionnum] + 'region_finalmask.nc')

#define number of members
members = [1,2,3,4,5]
#define variables
varnames = ['prcp', 'tdew', 'tmean', 'trange']
#time period
timeperiodstr = '195001-196912'
#timeperiodstr = '200001-201912'

dsfortime = xr.open_dataset(datadir + 'monthly/' + 'EM_Earth_probabilistic_monthly_prcp_NorthAmerica_' + timeperiodstr + '_001.nc')

#time len
timelen = len(dsfortime['time'].values)

#create xarray for 1950-1969
_arrsave = xr.Dataset(coords = {'member': members, 'time': dsfortime['time'].values}) 

os.chdir(datadir + 'monthly')
varfilenames = os.listdir('.')

#for each ensemble member and variable, output the area-weighted regional average time series
for varnum in range(0, len(varnames)):
    currentvariable = varnames[varnum]
    print(currentvariable)
    
    currentvarfilenames = [filename for filename in varfilenames if currentvariable + '_' in filename]
    currentvarfilenames = [filename for filename in currentvarfilenames if timeperiodstr in filename]
    currentvarfilenames = natsorted(currentvarfilenames)
    
    #print(varfilenames)
    
    arr_regionaltimeseries = np.empty((len(members), timelen))
    arr_regionaltimeseries[:] = np.nan

    
    for membernum in range(0, len(members)):
        currentmember = members[membernum]
        print(currentmember)
                
        #import and combine historical and ssp585 files
        currentds = xr.open_dataset(currentvarfilenames[membernum])
        #cut to exclude past -180 longitude
        currentds = currentds[currentvariable][:,:,77:1764]
        #flip latitude to sort by ascending
        currentds = currentds.sortby('lat', ascending = True)
        
        #calculate area-averaged weights with cos
        weights = np.cos(np.deg2rad(currentds.lat))
        weights.name = 'weights'
                
        currentds = currentds.where(land['mask'] == 1)

        #weight by land fraction and grid cell area
        ds_weighted = currentds.weighted(weights)
        
        ds_regionalmean = ds_weighted.mean(dim = ('lat', 'lon'), skipna = True)
                    
        #save in array
        arr_regionaltimeseries[membernum, :] = ds_regionalmean
    
    #save in xarray
    _arrsave[currentvariable] = (['member', 'time'], arr_regionaltimeseries)

#export arr_rollingglobalmean
_arrsave.to_netcdf(datadir + 'EM_Earth_probabilistic_monthly_CNAregion_' + timeperiodstr + 'meantimeseries.nc')                  


#create functions for tmax
def calculate_tmax(tmean, trange):
    
    tmax = tmean + trange/2   
    return tmax

def calculate_es(tmean):
    #convert tmean from C to K
    #tmean = tmean + 273.15
    
    #calculate saturation vapor presssure
    #es = 6.11*np.exp((2.5e6/461)*((1/273) - (1/tmean)))     
    es = 6.112*np.exp((17.67*tmean)/(tmean + 243.5))
    return es
    
def calculate_ea(tdew):        
    #calculate actual vapor pressure from Bolton (1980) - using celsius
    ea = 6.112*np.exp((17.67*tdew)/(tdew + 243.5))
    return ea

def calculate_vpd(es, ea):    
    #vapor pressure deficit
    vpd = es - ea
    return vpd

#timeperiodstr = '195001-196912'
timeperiodstr = '200001-201912'

#import regional time series
ts = xr.open_dataset(datadir + 'EM_Earth_probabilistic_monthly_CNAregion_' + timeperiodstr + 'meantimeseries.nc')  

#find tmax
ts['tmax'] = (['member', 'time'], calculate_tmax(ts['tmean'], ts['trange']))

#find saturated vapor presure
ts['es'] = (['member', 'time'], calculate_es(ts['tmean']))

#find actual vapor pressure   
ts['ea'] = (['member', 'time'], calculate_ea(ts['tdew']))

#find vpd
ts['vpd'] = (['member', 'time'], calculate_vpd(ts['es'], ts['ea']))   

#export
ts.to_netcdf(datadir + 'EM_Earth_probabilistic_monthly_CNAregion_' + timeperiodstr + 'meantimeseries_allalt.nc')                  

    