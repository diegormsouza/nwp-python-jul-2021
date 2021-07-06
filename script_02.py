#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC Training: NWP Data Processing With Python - Script 2: Basic Plot
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------
import pygrib                        # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
import matplotlib.pyplot as plt      # Plotting library
#----------------------------------------------------------------------------------------------------------- 

# Open the GRIB file
grib = pygrib.open("gfs.t00z.pgrb2full.0p50.f000")
 
# Select the variable
grb = grib.select(name='2 metre temperature')[0]

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-93.0, -60.00, -25.00, 18.00]

# Read the data for the selected extent
tmtmp, lats, lons = grb.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(8,8))
 
# Plot the image
plt.imshow(tmtmp, origin='lower', cmap='jet')
#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('image_2.png')

# Show the image
plt.show()