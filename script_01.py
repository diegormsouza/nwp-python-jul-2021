#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: NWP Data Processing With Python - Script 1: Knowing the Available Variables
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------
# Required modules
import pygrib                        # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
#-----------------------------------------------------------------------------------------------------------

# Image path
path = ('gfs.t00z.pgrb2full.0p50.f000')

# Open the GRIB file
grib = pygrib.open(path)

# Print all variables in the terminal and save them in a text file 
f = open("variables.txt", "w") # Create and open the file
for variables in grib:
    # Put the variables in the txt file
    print(variables, file=f)
    # Print the variables in the terminal
    print(variables)
f.close() # Close the file