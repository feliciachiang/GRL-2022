#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function that returns information from files downloaded from the
Earth System Grid Federation (ESGF) system 
Packages used: numpy, re

Created: Tues March 16 2021
Last edited: Tues March 16 2021

@author: Felicia Chiang - Contact: felicia.chiang@nasa.gov
"""

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