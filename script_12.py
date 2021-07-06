#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC Training: NWP Data Processing With Python - Script 12: Precipitation 
# Author: Diego Souza 
#-----------------------------------------------------------------------------------------------------------
import pygrib                              # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
import matplotlib.pyplot as plt            # Plotting library
import cartopy, cartopy.crs as ccrs        # Plot maps
import cartopy.io.shapereader as shpreader # Import shapefiles
import numpy as np                         # Scientific computing with Python
import matplotlib                          # Comprehensive library for creating static, animated, and interactive visualizations in Python 
#-----------------------------------------------------------------------------------------------------------  

# Open the GRIB file
grib = pygrib.open("gfs.t00z.pgrb2full.0p50.f024")
 
# Read the instant precip 
grb = grib.select(name='Precipitation rate', typeOfLevel = 'surface')[0]

# Get information from the file    
init  = str(grb.analDate)      # Init date / time
run   = str(grb.hour).zfill(2) # Run
ftime = str(grb.forecastTime)  # Forecast hour
valid = str(grb.validDate)     # Valid date / time 
print('Init: ' + init + ' UTC')
print('Run: ' + run + 'Z')
print('Forecast: +' + ftime)
print('Valid: ' + valid + ' UTC')

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-93.0, -60.00, -25.00, 18.00]

# Read the data for a specific region
precip, lats, lons = grb.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

#-----------------------------------------------------------------------------------------------------------

# Convert from kg m**-2 s**-1 to mm/h
precip = precip * 60 * 60

# To smooth the contours
import scipy.ndimage
precip = scipy.ndimage.zoom(precip, 3)
lats = scipy.ndimage.zoom(lats, 3)
lons = scipy.ndimage.zoom(lons, 3)

#-----------------------------------------------------------------------------------------------------------

# Select the variable
totpr = grib.select(name='Total Precipitation', typeOfLevel = 'surface')[1]

# Read the data for a specific region
totpr = totpr.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)[0]

# To smooth the contours
totpr = scipy.ndimage.zoom(totpr, 3)

#-----------------------------------------------------------------------------------------------------------

# Choose the plot size (width x height, in inches)
fig, axs = plt.subplots(1,2, figsize=(10,5), sharex = False, sharey = False, subplot_kw=dict(projection=ccrs.PlateCarree())) # 1 row x 2 columns

#-----------------------------------------------------------------------------------------------------------

# Define the image extent
axs[0].set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Add a background image
import cartopy.feature as cfeature
land = axs[0].add_feature(cfeature.LAND, facecolor='whitesmoke')
ocean = axs[0].add_feature(cfeature.OCEAN, facecolor='white')

# Add a shapefile
# https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
shapefile = list(shpreader.Reader('BR_UF_2019.shp').geometries())
axs[0].add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray',facecolor='none', linewidth=0.3)

# Add coastlines, borders and gridlines
axs[0].coastlines(resolution='10m', color='black', linewidth=0.8)
axs[0].add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = axs[0].gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Define de contour interval
data_min = 0.1
data_max = 50
interval = 1
levels = np.arange(data_min,data_max,interval)

# Create a custom color palette 
colors = ["#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2", "#0fa00f", "#28be28", "#50f050", "#72f06e", "#b3faaa", "#fff9aa", "#ffe978", "#ffc13c", "#ffa200", "#ff6200", "#ff3300", "#ff1500", "#c00100", "#a50200", "#870000", "#653b32"]
cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_over('#000000')
cmap.set_under('#ffffff')

# Plot the contours
from matplotlib import colors, cm
#norm = colors.LogNorm(0.1, 10, clip='False')

# make the norm:  Note the center is offset so that the land has more
# dynamic range:
#norm = colors.TwoSlopeNorm(vmin=0.1, vcenter=10, vmax=50)
#norm = colors.FuncNorm((_forward, _inverse), vmin=0, vmax=50)
#norm = colors.PowerNorm(gamma=0.5, vmin=0, vmax=50)
#img1 = axs[0].contourf(lons, lats, precip, norm=norm, cmap=cmap, levels=levels, extend='max')    

img1 = axs[0].contourf(lons, lats, precip, cmap=cmap, levels=levels, extend='max') 

# Add a colorbar
#ticks = np.arange(0, 50, 5).tolist()     
ticks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# Add a colorbar
plt.colorbar(img1, label='Instant Precipitation Rate (mm/h)', orientation='horizontal', pad=0.02, fraction=0.05, ticks=ticks, ax=axs[0])

# Add a title
axs[0].set_title('GFS: Instant Precipitation Rate (mm/h)' , fontweight='bold', fontsize=6, loc='left')
axs[0].set_title('Valid: ' + valid, fontsize=6, loc='right')

#-----------------------------------------------------------------------------------------------------------

# Define the image extent
axs[1].set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Add a background image
import cartopy.feature as cfeature
land = axs[1].add_feature(cfeature.LAND, facecolor='whitesmoke')
ocean = axs[1].add_feature(cfeature.OCEAN, facecolor='white')

# Add a shapefile
# https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
shapefile = list(shpreader.Reader('BR_UF_2019.shp').geometries())
axs[1].add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray',facecolor='none', linewidth=0.3)

# Add coastlines, borders and gridlines
axs[1].coastlines(resolution='10m', color='black', linewidth=0.8)
axs[1].add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = axs[1].gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Define de contour interval
data_min = 1
data_max = 100 
interval = 5
levels = np.arange(data_min,data_max,interval)

# Create a custom color palette 
colors = ["#b4f0f0", "#96d2fa", "#78b9fa", "#3c95f5", "#1e6deb", "#1463d2", "#0fa00f", "#28be28", "#50f050", "#72f06e", "#b3faaa", "#fff9aa", "#ffe978", "#ffc13c", "#ffa200", "#ff6200", "#ff3300", "#ff1500", "#c00100", "#a50200", "#870000", "#653b32"]

from matplotlib import cm                                          # Colormap handling utilities
from matplotlib.colors import LinearSegmentedColormap              # Linear interpolation for color maps
my_colors = cm.colors.LinearSegmentedColormap.from_list("",colors) # Create a custom colormap
my_colors = my_colors(np.linspace(0, 1, 256))                      # Create the array
my_colors[0:2,-1] = 0.0 
cmap = LinearSegmentedColormap.from_list(name='my_cmap', colors=my_colors)

cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_over('#000000')
cmap.set_under('#828282')

 # Create the color scale for low values
colors2 = ["#bebebe", "#a5a5a5", "#969696", "#828282"]

my_colors = cm.colors.LinearSegmentedColormap.from_list("",colors2) # Create a custom colormap
my_colors = my_colors(np.linspace(0, 1, 256))                       # Create the array
my_colors[0:80,-1] = 0.0
cmap2 = LinearSegmentedColormap.from_list(name='my_cmap', colors=my_colors)
cmap2.set_over('#828282')

vmin2 = 0.0
vmax2 = 1.0
thick_interval2 = 0.1

# Plot the contours (low values)
data_min2 = 0.0
data_max2 = 2.0
interval2 = 0.2
levels2 = np.arange(data_min2,data_max2,interval2)
img2 = axs[1].contourf(lons, lats, totpr, cmap=cmap2, levels=levels2)

# Plot the contours (high values)
img3 = axs[1].contourf(lons, lats, totpr, cmap=cmap, levels=levels, extend='max')    
img4 = axs[1].contour(lons, lats, totpr, colors='white', linewidths=0.3, levels=levels)
#axs[1].clabel(img4, inline=1, inline_spacing=0, fontsize='8',fmt = '%1.0f', colors= 'black')

ticks = np.arange(0, 100, 5).tolist() 

# Add a colorbar
plt.colorbar(img3, label='Total Precipitation (mm - 24h)', orientation='horizontal', pad=0.02, fraction=0.05,  ticks=ticks, ax=axs[1])

# Add a title
axs[1].set_title('GFS: Total Precipitation (mm - 24h)' , fontweight='bold', fontsize=6, loc='left')
axs[1].set_title('Valid: ' + valid, fontsize=6, loc='right')

#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('image_12.png', bbox_inches='tight', pad_inches=0, dpi=100)

# Show the image
plt.show()