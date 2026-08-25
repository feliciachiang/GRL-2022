#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function that returns relevant filenames given the variable, forcing, and member
Function that assembles the files from the filelist and outputs the complete xarray

Created: Tues March 16 2021
Last edited: Wed March 09 2022

@author: Felicia Chiang - Contact: felicia.chiang@nasa.gov
"""

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