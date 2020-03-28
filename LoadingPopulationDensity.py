# -*- coding: utf-8 -*-
"""
this script loads the population density and maps them onto a grid with 1 
degree resolution

Created on Sat Mar 28 10:59:16 2020

@author: Roland
"""



#%% libraries
import numpy as np
import gdal



#%% get the data
FileName = 'population_density\\gpw_v4_population_density_rev11_2020_1_deg.tif'
ds = gdal.Open(FileName)
geoMatrix = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
arr = band.ReadAsArray()
del band, ds

# get the details of the transform
ulX = geoMatrix[0]
ulY = geoMatrix[3]
xDist = geoMatrix[1]
yDist = geoMatrix[5]
rtnX = geoMatrix[2]
rtnY = geoMatrix[4]
Ny, Nx = arr.shape



#%% aggregate to grid
# resolution 
Res = 1

# make the grid for the world
lon = np.arange(-180, 180, 1)
lat = np.arange(-90,90,1)
Lon, Lat = np.meshgrid(lon, lat)

# fill the grid
