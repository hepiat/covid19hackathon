# -*- coding: utf-8 -*-
"""
this script loads the population density and maps them onto a grid with 1 
degree resolution

Created on Sat Mar 28 10:59:16 2020

@author: Roland, Juhyun
"""



#%% libraries
import numpy as np
import gdal
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw


#%% get the data
FileName = 'gpw_v4_population_density_rev11_2020_1_deg.tif'
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

plt.imshow(np.log10(arr))

#%% aggregate to grid
# resolution 
Res = 1

# make the grid for the world
lon = np.arange(-180, 180, 1)
lat = np.arange(-90,90,1)
Lon, Lat = np.meshgrid(lon, lat)

# fill the grid



#%% new trial
import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.pyplot as plt
#import seaborn as sns

world_path ="./ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"
world = shp.Reader(world_path)

def read_shapefile(sf):
    """
    Read a shapefile into a Pandas dataframe with a 'coords' 
    column holding the geometry information. This uses the pyshp
    package
    """
    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df


df = read_shapefile(world)
df.shape


#%%Draw on the grid
def plot_shape(id, s=None):
    """ PLOTS A SINGLE SHAPE """
    plt.figure()
    ax = plt.axes()
    ax.set_aspect('equal')
    shape_ex = world.shape(id)
    x_lon = np.zeros((len(shape_ex.points),1))
    y_lat = np.zeros((len(shape_ex.points),1))
    for ip in range(len(shape_ex.points)):
        x_lon[ip] = shape_ex.points[ip][0]
        y_lat[ip] = shape_ex.points[ip][1]
    plt.plot(x_lon,y_lat) 
    x0 = np.mean(x_lon)
    y0 = np.mean(y_lat)
    plt.text(x0, y0, s, fontsize=10)
    # use bbox (bounding box) to set plot limits
    plt.xlim(shape_ex.bbox[0],shape_ex.bbox[2])
    return x0, y0



#%%Draw on the grid
def plot_map(world, x_lim = None, y_lim = None, figsize = (11,9)):
    '''
    Plot map with lim coordinates
    '''
    plt.figure(figsize = figsize)
    id=0
    for shape in world.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        plt.plot(x, y, 'k')
        
        if (x_lim == None) & (y_lim == None):
            x0 = np.mean(x)
            y0 = np.mean(y)
            plt.text(x0, y0, id, fontsize=10)
        id = id+1
    
    if (x_lim != None) & (y_lim != None):     
        plt.xlim(x_lim)
        plt.ylim(y_lim)


#%%Draw on the grid
plot_map(world)



#%% get in onto the grid
# define the grid 
lon = np.arange(-180, 180, 1)
lat = np.arange(-90,90,1)
Lat, Lon = np.meshgrid(lon, lat)

# loop through the countries 
MeshAll = np.zeros((Lon.shape[0], Lon.shape[1], df.shape[0]))
for ii in range(df.shape[0]):
    
    # current country
    print(df['ISO_A3'][ii])
    
    # getting the coordinates
    coords = np.array(df['coords'][ii])
    
    # transformation
    coords_trans = np.zeros(coords.shape)
    coords_trans[:,0] = coords[:,1] + 90
    coords_trans[:,1] = coords[:,0] + 180

    # try again without spliting the map
    #img = Image.new('L', (Lat.shape), 0)
    #ImageDraw.Draw(img).polygon(tuple(map(tuple, (coords_trans))), outline=1, fill=1)
    #mm = np.array(img)
    
    # identify the different sections by looking at recurring points
    co = coords_trans
    IDstart = list()
    idstart = 0
    IDend = list()
    while idstart < co.shape[0]:
        co1 = co[idstart,:]
        
        idloop = np.arange(co.shape[0]) * ((co[:,0] == co1[0]) * (co[:,1] == co1[1]))
        idloop = idloop[idloop>0]
        
        # collect the indices
        IDstart.append(idstart)
        IDend.append(idloop[-1])
        
        # reset id start
        idstart = idloop[-1] + 1
            
    # loop over the different polygons in this country
    Mask = np.zeros((Lon.shape[1],Lon.shape[0],len(IDstart)))
    for pp in range(len(IDstart)):    
        # create a new image
        img = Image.new('L', (Lat.shape), 0)
    
        # define the polygon
        ImageDraw.Draw(img).polygon(tuple(map(tuple, coords_trans[IDstart[pp]:(IDend[pp]+1),:])), outline=1, fill=1)
        
        # make the mask
        mask = np.array(img)
        Mask[:,:,pp] = mask
        
    Mask = np.sum(Mask, axis=2) 
    Mask[Mask > 1] = 1
        
    MeshAll[:,:,ii] = np.transpose(Mask)
    
    # show the image for checking
    #plt.figure()
    #plt.imshow(np.transpose(mask), origin='lower')


# save the mask for later usage
np.save('Mask', MeshAll) 
np.save('Countries', df)   



#%% plot a single country for checking
plt.figure()
for pp in range(len(IDstart)):
    
    plt.plot(coords_trans[IDstart[pp]:(IDend[pp]+1),0], coords_trans[IDstart[pp]:(IDend[pp]+1),1])