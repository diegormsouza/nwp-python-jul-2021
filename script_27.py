#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: Python and GOES-R Imagery: Script 27 - METAR + NWP + Satellite Plot
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------

# Required modules
from netCDF4 import Dataset                                     # Read / Write NetCDF4 files
from osgeo import gdal                                          # Python bindings for GDAL
import matplotlib.pyplot as plt                                 # Plotting library
import cartopy, cartopy.crs as ccrs                             # Plot maps
import cartopy.io.shapereader as shpreader                      # Import shapefiles
import os                                                       # Miscellaneous operating system interfaces
import numpy as np                                              # Scientific computing with Python
from matplotlib import cm                                       # Colormap handling utilities
from utilities import download_CMI                              # Our function for download
from utilities import reproject                                 # Our function for reproject
from utilities import loadCPT                                   # Import the CPT convert function
import requests                                                 # HTTP library for Python
from datetime import timedelta, date, datetime                  # Basic Dates and time types
from metpy.calc import reduce_point_density                     # Provide tools for unit-aware, meteorological calculations    
from metpy.io import metar                                      # Parse METAR-formatted data
from metpy.plots import current_weather, sky_cover, StationPlot # Contains functionality for making meteorological plots
import pygrib                                                   # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
gdal.PushErrorHandler('CPLQuietErrorHandler')                   # Ignore GDAL warnings

#-----------------------------------------------------------------------------------------------------------

def plot_maxmin_points(lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
    """
    This function will find and plot relative maximum and minimum for a 2D grid. The function
    can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
    values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
    max/min value. The symbol text can be set to a string value and optionally the color of the
    symbol and any plotted value can be set with the parameter color
    lon = plotting longitude values (2D)
    lat = plotting latitude values (2D)
    data = 2D data that you wish to plot the max/min symbol placement
    extrema = Either a value of max for Maximum Values or min for Minimum Values
    nsize = Size of the grid box to filter the max and min values to plot a reasonable number
    symbol = String to be placed at location of max/min value
    color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
    plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
    The max/min symbol will be plotted on the current axes within the bounding frame
    (e.g., clip_on=True)
    """
    from scipy.ndimage.filters import maximum_filter, minimum_filter

    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    for i in range(len(mxy)):
        txt1 = ax.annotate(symbol, xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), color=color, size=20,
                clip_on=True, annotation_clip=True, horizontalalignment='center', verticalalignment='center',
                transform=ccrs.PlateCarree())

        txt2 = ax.annotate('\n' + str(int(data[mxy[i], mxx[i]])), xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), 
                color=color, size=10, clip_on=True, annotation_clip=True, fontweight='bold', horizontalalignment='center', verticalalignment='top',
                transform=ccrs.PlateCarree())  

#----------------------------------------------------------------------------------------------------------- 

# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-93.0, -60.00, -25.00, 18.00] # Min lon, Max lon, Min lat, Max lat

# Datetime to process (today in this example, to match the GFS date)
#date = datetime.today().strftime('%Y%m%d')
#yyyymmddhhmn = date + '0000'
yyyymmddhhmn = '202107020000' # CHANGE THIS DATE TO THE SAME DATE OF YOUR NWP DATA

#-----------------------------------------------------------------------------------------------------------
# Get the Band 13 Data

# Download the file
file_ir = download_CMI(yyyymmddhhmn, 13, input)
#-----------------------------------------------------------------------------------------------------------
# Variable
var = 'CMI'

# Open the file
img = gdal.Open(f'NETCDF:{input}/{file_ir}.nc:' + var)

# Read the header metadata
metadata = img.GetMetadata()
scale = float(metadata.get(var + '#scale_factor'))
offset = float(metadata.get(var + '#add_offset'))
undef = float(metadata.get(var + '#_FillValue'))
dtime = metadata.get('NC_GLOBAL#time_coverage_start')

# Load the data
ds_cmi = img.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize).astype(float)

# Apply the scale, offset and convert to celsius
ds_cmi = (ds_cmi * scale + offset) - 273.15

# Reproject the file
filename_ret = f'{output}/IR_{yyyymmddhhmn}.nc'
reproject(filename_ret, img, ds_cmi, extent, undef)

# Open the reprojected GOES-R image
file = Dataset(filename_ret)

# Get the pixel values
data = file.variables['Band1'][:]

#----------------------------------------------------------------------------------------------------------- 

# Open the GRIB file
grib = pygrib.open("gfs.t00z.pgrb2full.0p50.f000")
 
# Select the variable
prmls = grib.select(name='Pressure reduced to MSL')[0]

# Get information from the file    
init  = str(prmls.analDate)      # Init date / time
run   = str(prmls.hour).zfill(2) # Run
ftime = str(prmls.forecastTime)  # Forecast hour
valid = str(prmls.validDate)     # Valid date / time 
print('Init: ' + init + ' UTC')
print('Run: ' + run + 'Z')
print('Forecast: +' + ftime)
print('Valid: ' + valid + ' UTC')

# Read the data for a specific region
prmls, lats, lons = prmls.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

# Convert to hPa
prmls = prmls / 100

#-----------------------------------------------------------------------------------------------------------

# Input directory
dir = "Samples"; os.makedirs(dir, exist_ok=True)

# Download the METAR File
#date = datetime.today().strftime('%Y%m%d')
date = '20210702' # CHANGE THIS DATE TO THE SAME DATE OF YOUR NWP DATA
url = 'https://thredds-test.unidata.ucar.edu/thredds/fileServer/noaaport/text/metar' 
file_name = 'metar_' + date + '_0000.txt'

# Sends a GET request to the specified url
myfile = requests.get(url + '//' + file_name)

# Download the file
open(dir + '//' + file_name, 'wb').write(myfile.content)

# METAR File
# https://unidata.github.io/MetPy/latest/examples/plots/Station_Plot.html
data_metar = metar.parse_metar_file(dir + '//' + file_name)

# Drop rows with missing winds
data_metar = data_metar.dropna(how='any', subset=['wind_direction', 'wind_speed'])

#-----------------------------------------------------------------------------------------------------------

# Set up the map projection
proj = ccrs.PlateCarree()

# Create the figure and an axes set to the projection.
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Change the DPI of the resulting figure. Higher DPI drastically improves the
# look of the text rendering.
plt.rcParams['savefig.dpi'] = 255

# Use the Cartopy map projection to transform station locations to the map and
# then refine the number of stations plotted by setting a minimum radius
point_locs = proj.transform_points(ccrs.PlateCarree(), data_metar['longitude'].values, data_metar['latitude'].values)
data_metar = data_metar[reduce_point_density(point_locs, 3)]

# Define the color scale based on the channel
colormap = "gray_r" # White to black for IR channels
    
# Plot the image
img1 = ax.imshow(data, origin='upper', vmin=-80, vmax=60, extent=img_extent, cmap=colormap, alpha=1.0)

# Define de contour interval
data_min = 500
data_max = 1050
interval = 2
levels = np.arange(data_min,data_max,interval)

# Plot the contours
img2 = ax.contour(lons, lats, prmls, colors='white', linewidths=0.7, levels=levels)
ax.clabel(img2, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'gray')

# Use definition to plot H/L symbols
plot_maxmin_points(lons, lats, prmls, 'max', 50, symbol='H', color='b',  transform=ccrs.PlateCarree())
plot_maxmin_points(lons, lats, prmls, 'min', 25, symbol='L', color='r', transform=ccrs.PlateCarree())

# Add a shapefile
# https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
shapefile = list(shpreader.Reader('BR_UF_2019.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='gray',facecolor='none', linewidth=0.3)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='black', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

#-----------------------------------------------------------------------------------------------------------
# Station Plot

# Start the station plot by specifying the axes to draw on, as well as the
# lon/lat of the stations (with transform). We also the fontsize to 12 pt.
stationplot = StationPlot(ax, data_metar['longitude'].values, data_metar['latitude'].values,
                          clip_on=True, transform=ccrs.PlateCarree(), fontsize=8)

# Plot the temperature and dew point to the upper and lower left, respectively, of
# the center point. Each one uses a different color.
stationplot.plot_parameter('NW', data_metar['air_temperature'].values, color='red')
stationplot.plot_parameter('SW', data_metar['dew_point_temperature'].values,
                           color='darkgreen')

# A more complex example uses a custom formatter to control how the sea-level pressure
# values are plotted. This uses the standard trailing 3-digits of the pressure value
# in tenths of millibars.
stationplot.plot_parameter('NE', data_metar['air_pressure_at_sea_level'].values,
                           formatter=lambda v: format(10 * v, '.0f')[-3:])

# Plot the cloud cover symbols in the center location. This uses the codes made above and
# uses the `sky_cover` mapper to convert these values to font codes for the
# weather symbol font.
stationplot.plot_symbol('C', data_metar['cloud_coverage'].values, sky_cover)

# Same this time, but plot current weather to the left of center, using the
# `current_weather` mapper to convert symbols to the right glyphs.
stationplot.plot_symbol('W', data_metar['present_weather'].values, current_weather)

# Add wind barbs
#stationplot.plot_barb(data_metar['eastward_wind'].values, data_metar['northward_wind'].values)

# Also plot the actual text of the station id. Instead of cardinal directions,
# plot further out by specifying a location of 2 increments in x and 0 in y.
stationplot.plot_text((2, 0), data_metar['station_id'].values)

#-----------------------------------------------------------------------------------------------------------

# Extract the date from the satellite data
date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))

# Add a title
plt.title('GOES-16 Band 13 ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC' + ' + METAR + PSML (hPa)', fontweight='bold', fontsize=6, loc='left')
plt.title('Reg.: ' + str(extent) , fontsize=6, loc='right')

# Save the image
plt.savefig(f'{output}/image_27.png', bbox_inches='tight', pad_inches=0, dpi=300)

plt.show()