#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create map of Central North America region
Output jpg map

@author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import packages
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature

#save figure in
figdir = "/Users/fchiang/GISS/HSM project/figures/"

#create plot
ax = plt.axes(projection=ccrs.PlateCarree())

plt.title('Central North America')
ax.coastlines(resolution = '50m')

ax.set_extent([-140, -50, 20, 70], ccrs.PlateCarree())

#plt.show()



#access data
wdir = "/Volumes/GISS/"

#import resources (such as crop maps)
resourcedir = wdir + "resources/"

reader = shpreader.Reader(resourcedir + 'referenceRegions/referenceRegions.shp')

CNA = [region for region in reader.records() if region.attributes['LAB'] == 'CNA'][0]

#for CNA region:
shape_feature = ShapelyFeature([CNA.geometry], ccrs.PlateCarree(), facecolor = "None", edgecolor = 'blue', lw = 3)
#add shape
ax.add_feature(shape_feature)

#save figure 
plt.savefig(figdir + "CNAreferencemap.jpg", dpi = 600)
