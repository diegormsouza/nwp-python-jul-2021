#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: Python and GOES-R Imagery: Script 26 - METAR + NWP Plot
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------

# Adapted from: https://unidata.github.io/MetPy/latest/examples/plots/Station_Plot.html

import matplotlib.pyplot as plt                                 # Plotting library
import cartopy, cartopy.crs as ccrs                             # Plot maps
import cartopy.io.shapereader as shpreader                      # Import shapefiles
import cartopy.feature as cfeature                              # Common drawing and filtering operations
import os                                                       # Miscellaneous operating system interfaces
import numpy as np                                              # Scientific computing with Python
import requests                                                 # HTTP library for Python
from datetime import timedelta, date, datetime                  # Basic Dates and time types
from metpy.calc import reduce_point_density                     # Provide tools for unit-aware, meteorological calculations    
from metpy.io import metar                                      # Parse METAR-formatted data
from metpy.plots import current_weather, sky_cover, StationPlot # Contains functionality for making meteorological plots
import pygrib                                                   # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
                   
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

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-93.0, -60.00, -25.00, 18.00]

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

# Input and output directories
dir = "Samples"; os.makedirs(dir, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

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
data = metar.parse_metar_file(dir + '//' + file_name)

# Drop rows with missing winds
data = data.dropna(how='any', subset=['wind_direction', 'wind_speed'])

#-----------------------------------------------------------------------------------------------------------

# Choose the plot size (width x height, in inches)
plt.figure(figsize=(8,8))

# Set up the map projection
proj = ccrs.PlateCarree()

# Use the Geostationary projection in cartopy
ax = plt.axes(projection=proj)

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Change the DPI of the resulting figure. Higher DPI drastically improves the
# look of the text rendering.
plt.rcParams['savefig.dpi'] = 255

# Use the Cartopy map projection to transform station locations to the map and
# then refine the number of stations plotted by setting a minimum radius
point_locs = proj.transform_points(ccrs.PlateCarree(), data['longitude'].values, data['latitude'].values)
data = data[reduce_point_density(point_locs, 3)]

# Add some various map elements to the plot to make it recognizable.
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

# Define de contour interval
data_min = 500
data_max = 1050
interval = 2
levels = np.arange(data_min,data_max,interval)

# Plot the contours
img1 = ax.contour(lons, lats, prmls, colors='gray', linewidths=0.7, levels=levels)
ax.clabel(img1, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'gray')

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
stationplot = StationPlot(ax, data['longitude'].values, data['latitude'].values,
                          clip_on=True, transform=ccrs.PlateCarree(), fontsize=8)

# Plot the temperature and dew point to the upper and lower left, respectively, of
# the center point. Each one uses a different color.
stationplot.plot_parameter('NW', data['air_temperature'].values, color='red')
stationplot.plot_parameter('SW', data['dew_point_temperature'].values,
                           color='darkgreen')

# A more complex example uses a custom formatter to control how the sea-level pressure
# values are plotted. This uses the standard trailing 3-digits of the pressure value
# in tenths of millibars.
stationplot.plot_parameter('NE', data['air_pressure_at_sea_level'].values,
                           formatter=lambda v: format(10 * v, '.0f')[-3:])

# Plot the cloud cover symbols in the center location. This uses the codes made above and
# uses the `sky_cover` mapper to convert these values to font codes for the
# weather symbol font.
stationplot.plot_symbol('C', data['cloud_coverage'].values, sky_cover)

# Same this time, but plot current weather to the left of center, using the
# `current_weather` mapper to convert symbols to the right glyphs.
stationplot.plot_symbol('W', data['present_weather'].values, current_weather)

# Add wind barbs
#stationplot.plot_barb(data['eastward_wind'].values, data['northward_wind'].values)

# Also plot the actual text of the station id. Instead of cardinal directions,
# plot further out by specifying a location of 2 increments in x and 0 in y.
stationplot.plot_text((2, 0), data['station_id'].values)

# Add a title
plt.title('METAR + GFS PSML (hPa) | ' + date + ' 00:00 UTC', fontsize=8, loc='center')
#-----------------------------------------------------------------------------------------------------------

# Save the image
plt.savefig(f'{output}/image_26.png', bbox_inches='tight', pad_inches=0, dpi=300)

plt.show()