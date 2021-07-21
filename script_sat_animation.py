#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: Python and GOES-R Imagery: Script 12 - Cropping the Full Disk + Animation
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset                    # Read / Write NetCDF4 files
import matplotlib.pyplot as plt                # Plotting library
from datetime import datetime                  # Basic Dates and time types
from datetime import timedelta, date, datetime # Basic Dates and time types
import cartopy, cartopy.crs as ccrs            # Plot maps
import cartopy.io.shapereader as shpreader     # Import shapefiles
import os                                      # Miscellaneous operating system interfaces
from utilities import download_CMI             # Our own utilities
from utilities import geo2grid, latlon2xy, convertExtent2GOESProjection      # Our own utilities
#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "/content/Samples"; os.makedirs(input, exist_ok=True)
output = "/content/Animation"; os.makedirs(output, exist_ok=True)

# Desired extent
extent = [-64.0, -36.0, -40.0, -15.0] # Min lon, Max lon, Min lat, Max lat

# Datetime to process
yyyymmddhhmn = '202102181800'
band = '13'

# Initial time and date
yyyy = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
mm = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%m')
dd = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%d')
hh = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
mn = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

date_ini = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)))
date_end = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)) + timedelta(hours=2))

while (date_ini <= date_end):

    # Date structure
    yyyymmddhhmn = datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M')
    print(yyyymmddhhmn)

    # Download the file
    file_name = download_CMI(yyyymmddhhmn, band, input)

    #-----------------------------------------------------------------------------------------------------------
    # Open the GOES-R image
    file = Dataset(f'{input}/{file_name}.nc')
                      
    # Convert lat/lon to grid-coordinates
    lly, llx = geo2grid(extent[1], extent[0], file)
    ury, urx = geo2grid(extent[3], extent[2], file)
            
    # Get the pixel values
    data = file.variables['CMI'][ury:lly, llx:urx] - 273.15      
    #-----------------------------------------------------------------------------------------------------------
    # Compute data-extent in GOES projection-coordinates
    img_extent = convertExtent2GOESProjection(extent)
    #-----------------------------------------------------------------------------------------------------------
    # Choose the plot size (width x height, in inches)
    plt.figure(figsize=(10,10))

    # Use the Geostationary projection in cartopy
    ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))

    # Define the color scale based on the channel
    colormap = "gray_r" # White to black for IR channels
            
    # Plot the image
    img = ax.imshow(data, origin='upper', vmin=-80, vmax=60, extent=img_extent, cmap=colormap)

    # Add a shapefile
    # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
    shapefile = list(shpreader.Reader('BR_UF_2019.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='white',facecolor='none', linewidth=0.3)

    # Add coastlines, borders and gridlines
    ax.coastlines(resolution='10m', color='white', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.5)
    ax.gridlines(color='white', alpha=0.5, linestyle='--', linewidth=0.5)

    # Add a colorbar
    plt.colorbar(img, label='Brightness Temperatures (°C)', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

    # Extract date
    date = (datetime.strptime(file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))

    # Add a title
    plt.title('GOES-16 Band 13 ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
    plt.title('Reg.: ' + str(extent) , fontsize=10, loc='right')
    #-----------------------------------------------------------------------------------------------------------
    # Save the image
    plt.savefig(f'{output}/G16_B{band}_{yyyymmddhhmn}.png', bbox_inches='tight', pad_inches=0, dpi=100)

    # Show the image
    plt.show()

    # Increment the date_ini
    date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(minutes=10))
