# -*- coding: utf-8 -*-
"""
Check that all files (all years) have been downloaded from ESGF
Specific use for: High soil moisture constraining temperatures and VPD 

@author: Felicia Chiang - Contact: felicia.chiang@nasa.gov
"""
#import libraries
import os
import numpy as np
import fnmatch
#import re
#import xarray as xr
import matplotlib.pyplot as plt
from natsort import natsorted

#access data
wdir = "/Users/fchiang/GISS/"
scriptdir = wdir + "HSM 2022 scripts/"

datadir = "/Volumes/GISS/CMIP6/"

# # #change directory to where functions are stored
# import sys
# sys.path.append(scriptdir + "functions/")

# #import useful functions
# #import HSM_CMIP6infofromfilelist
# #HSM_CMIP6infofromfilelist.__file__
# from HSM_CMIP6infofromfilelist import ESGFinfofromfilelist

# #import HSM_modelcheck
# #HSM_modelcheck.__file__
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

#list all models in alpha order
modelnames = sorted([d for d in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, d))])

#only look at the relevant variables
varnames = ['hurs', 'mrsos', 'pr', 'tas', 'tasmax']

#for each model
for modelnum in range(0, len(modelnames)):
    currentmodel = modelnames[modelnum]
    print(currentmodel)
    
    #check that all have all variables and all years
    currentdir = os.path.join(datadir, currentmodel, 'mon/')

    filenames = fnmatch.filter(os.listdir(currentdir), '*.nc')
      
    #list all unique variables (from GISS model)
    #varnames = ESGFinfofromfilelist(filenames, 'variable')
    # #list all unique forcings
    # forcings = ESGFinfofromfilelist(filenames, 'forcing')
    #list all unique ensemble members
    members = natsorted(ESGFinfofromfilelist(filenames, 'variant'))
    print(members)
    
    #varnum = -1


    #for each variable,
    for varnum in range(0, len(varnames)):
        
        #varnum = varnum + 1
        currentvariable = varnames[varnum]
        print(currentvariable)
            
        #for each member, check if all years are present in the files downloaded
        for variantnum in range(0, 5):
            currentmember = members[variantnum]
            print(currentmember)
                        
            #use missingyears function to check if missing years
            currentfilenames = relevantfilenames(filenames, currentvariable, currentmember)
            print(currentfilenames)
            
            timesteplen = len(np.arange(np.datetime64('1850-01-01'), np.datetime64('2101-01-01'), np.timedelta64(1, 'M'), dtype = 'datetime64[M]'))
            
            #import and combine multiple files
            ds = combineyears(currentfilenames, currentdir, timesteplen)
            
            #create simple time series plot of random location
            plt.figure()
            ds[currentvariable][:,50,20].plot()
            
            plt.figure()
            ds[currentvariable][:,30,20].plot()
            
            #create simple map of random time slice
            plt.figure()
            ds[currentvariable][10,:,:].plot()
            
