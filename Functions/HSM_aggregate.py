#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aggregate n months of data for each year

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

import numpy as np

# ----------------------------------------------------------------------------

def aggregate(xd, sc, endmonth, function):
    #xd is the time series starting in january
    #sc is the number of months to aggregate
    #endmonth is the last month of the aggregate
    
    A1 = np.empty((len(xd)-sc+1, sc))
    A1[:] = np.nan
    for i in range(0, sc):
        A1[:,i] = xd[i:len(xd)-sc+1+i]

    if function == "sum":
        X = np.nansum(A1, 1)
    elif function == "mean":
        X = np.nanmean(A1, 1)
        
    #add nan values to time series for sc > 1    
    nanvector = np.empty(sc-1)
    nanvector[:] = np.nan
    X = np.concatenate((nanvector, X))       
    
    #extract appropriate aggregate from monthly time series
    finalts = X[endmonth::12]
    
    #return annual time series 
    return finalts