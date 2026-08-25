#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to calculate es (saturation vapor pressure) and vpd (vapor pressure deficit)
Based on p350 from Global Physical Climatology By Dennis Hartmann

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

import numpy as np

# ----------------------------------------------------------------------------

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

