# -*- coding: utf-8 -*-
"""
Check that all files (all years) have been downloaded from ESGF
Output xarray netcdfs 

Method:
    1) Calculate global mean temperature for each month of each ensemble time series 
        -Weighted with the cos of the latitude
    2) Calculate rolling 20-year mean and extract annual values
    3) Find the first year of the 20-year period closest to the degree (1, 2, 3, 4)

@author: Felicia Chiang - Contact: felicia.chiang@nasa.gov
"""
#import libraries
import os
import numpy as np
import fnmatch
#import re
import xarray as xr
from natsort import natsorted
#import matplotlib.pyplot as plt

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"
resultsdir = wdir + "HSM project/results/"

datadir = "/Volumes/GISS/CMIP6/"

# #change directory to where functions are stored
# os.chdir(scriptdir + "functions/")
# #import useful functions
# from HSM_CMIP6infofromfilelist import ESGFinfofromfilelist
# from HSM_modelcheck import relevantfilenames, combineyears

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

def relevantfilenames(filenames, variable, member):
    '''Function that returns relevant filenames given the variable, forcing, and member'''
    import fnmatch
    
    #isolate relevant filenames from given variable, forcing, and ensemble member
    #list all files with current variable
    currentfilelist = fnmatch.filter(filenames, variable + '_*.nc')
    #find all files with current forcing
    #currentfilelist = fnmatch.filter(currentfilelist, '*' + forcing + '_*.nc')
    #find all files with current ensemble member
    currentfilelist = fnmatch.filter(currentfilelist, '*' + member + '_*.nc')

    return currentfilelist
    
def combineyears(filenames, path, timesteplen):
    '''Function that takes the files from the filelist and outputs the complete xarray'''
    import xarray as xr
    
    #sort filenames just in case read out of order
    filenames.sort()
    
    #add path to filenames
    filenameswpath = [path + filename for filename in filenames]
    
    #import all using xarray - open_mfdataset function
    ds = xr.open_mfdataset(filenameswpath, combine = 'by_coords', parallel = True, chunks = {'lat': 10, 'lon': 10})
    #same length?
    if (len(ds.indexes['time']) - timesteplen) == 0:
        return ds


# #----------------------EDIT INPUTS HERE--------------------------------------#
# #variable to check
# variable = "pr"

# #----------------------------------------------------------------------------#

#list all models in alpha order
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#for each model
for modelnum in range(0, len(modelnames)):
    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    #check that all have all variables and all years
    currentdir = os.path.join(datadir, currentmodel, 'mon/')

    filenames = fnmatch.filter(os.listdir(currentdir), '*.nc')
    filenames = [filename for filename in filenames if '._' not in filename]
      
    #list all unique variables (from GISS model)
    #varnames = ESGFinfofromfilelist(filenames, 'variable')
    # #list all unique forcings
    # forcings = ESGFinfofromfilelist(filenames, 'forcing')
    #list all unique ensemble members
    members = natsorted(ESGFinfofromfilelist(filenames, 'variant'))
    print(members)
    
    currentvariable = 'tas_'
    
    print(currentvariable)
    
    #list of filenames for current variable
    currentfilelist = [filename for filename in filenames if currentvariable in filename]
    
    window = 20
    
    arr_rollingglobalmean = np.empty((len(members), 251-window+1))
    arr_rollingglobalmean[:] = np.nan
    
    #for each member, check if all years are present in the files downloaded
    for variantnum in range(0, len(members)):
        currentmember = members[variantnum]
        print(currentmember)
        
        #currentmember = 'r1i1p1f2'
        
        #use missingyears function to check if missing years
        variantfilenames = fnmatch.filter(currentfilelist, '*' + currentmember + '*')
        print(variantfilenames)
        
        filenameswpath = [currentdir + filename for filename in variantfilenames]
        
        #import and combine historical and ssp585 files
        ds = xr.open_mfdataset(natsorted(filenameswpath), combine = 'by_coords')
        
        #calculate area-averaged weights
        weights = np.cos(np.deg2rad(ds.lat))
        weights.name = 'weights'
        
        ds = ds[currentvariable[0:len(currentvariable)-1]]
        #ds_land = ds[currentvariable[0:len(currentvariable)-1]].where(laf['sftlf']>50)
        
        ds_weighted = ds.weighted(weights)
        
        #ds[currentvariable[0:len(currentvariable)-1]].where(laf['sftlf']>50)[0,:,:].plot()
        
        ds_globalmean = ds_weighted.mean(dim = ('lat', 'lon'), skipna = True)
        ds_globalmean = ds_globalmean.chunk({'time': None})
        
        #find 20-year rolling mean time series
        ds_first50globalmean = np.mean(ds_globalmean[0:50*12].values)
        ds_rollingglobalmean = ds_globalmean.rolling(time = window*12, center = False).mean()[11::12].values
        
        #find difference from the first 50 years (1850-1899)
        ds_finalanom = ds_rollingglobalmean - ds_first50globalmean
        
        #save in array
        arr_rollingglobalmean[variantnum, :] = ds_finalanom[window-1:len(ds_finalanom)]
    
    
    #export arr_rollingglobalmean
    _arrsave = xr.Dataset({'rollingglobalmeananom': (['member', 'startyear'], arr_rollingglobalmean)},
                          coords = {'member': members, 'startyear': np.arange(1850,2100-window+2)}) 
    
    _arrsave.to_netcdf(resultsdir + currentmodel + '_' + str(window) + 'year_rollingglobalmean_areaweighted.nc')                  
     
    
    #find starting year for 20-year at n degrees
    arr_startyeardeg = np.empty((len(members), 4))
    arr_startyeardeg[:] = np.nan
    
    #export years for each ensemble member (1, 2, 3, 4 degrees)
    for variantnum in range(0, len(members)):
        
        for degree in range(1,5):
            #print(degree)
            #find starting position for each degree for each ensemble member (number of years from 1850)
            arr_startyeardeg[variantnum, degree-1] = (np.abs(arr_rollingglobalmean[variantnum,:]-degree)).argmin()
            
    #export years csv
    features = xr.Dataset({'startyear': (['member', 'degree'], arr_startyeardeg)},
                          coords = {'member': members, 'degree': np.arange(1,5)})  
    
    #export
    features.to_netcdf(resultsdir + currentmodel + '_globalwarmingdegrees.nc') 


for modelnum in range(0, len(modelnames)):
    modelnum = modelnum + 1
    
    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    #import and review
    features = xr.open_dataset(resultsdir + currentmodel + '_globalwarmingdegrees.nc')  
    
    features['startyear'] + 1850         

