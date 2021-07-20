#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: Download SST, SST-Anomaly or SST-7-Day-Trend Data from a FTP
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------
# Required modules
from datetime import datetime, timedelta # Basic Dates and time types
import os                                # Miscellaneous operating system interfaces
import time as t                         # Time access and conversion                                          
from ftplib import FTP                   # FTP protocol client
#-----------------------------------------------------------------------------------------------------------

print('---------------------------------------')
print('SST NOAA FTP Download - Script started.')
print('---------------------------------------')

# Start the time counter
start_time = t.time()  

# Data description:

# SST: NETCDF / ~12 MB / 5 km / Since 1985
# Example FTP link
# ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/sst/2021/
# Example file name
# coraltemp_v3.1_20210102.nc

# SST Anomaly: NETCDF / ~12 MB / 5 km / Since 1985
# Example FTP link
# ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/ssta/2021/
# Example file name
# ct5km_ssta_v3.1_20210101.nc

# 7-Day SST Trend: NETCDF / ~12 MB / 5 km / Since 2019
# Example FTP link
# ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/sst-trend-7d/2021/
# Example file name
# ct5km_sst-trend-7d_v3.1_20210101.nc

#-----------------------------------------------------------------------------------------------------------

# Download directory
dir = "Samples"; os.makedirs(dir, exist_ok=True)

# Select model 
product = 'SST-A' # options: 'SST', 'SST-A', 'SST-7'

# Desired year (four digit)
#year = datetime.today().strftime('%Y') # This will get the current year (four digit)
# Or, select the year you want
year = '2021' 

# Desired month (two digit)
#month = datetime.today().strftime('%m') # This will get the current month (two digit)
# Or, select the month you want
month = '05'

# Desired day (two digit)
#day = datetime.today().strftime('%d') # This will get the current day (two digit)
# Or, select the day you want
day = '15' 

#-----------------------------------------------------------------------------------------------------------

# FTP Address
ftp = FTP('ftp.star.nesdis.noaa.gov') 

# FTP Credentials 
ftp.login('', '') 

# Access the FTP folder, based on the desired model
if (product == 'SST'):
    # FTP Path
    path = ('pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/sst/' + year + '/')
    # File Name
    file_name = 'coraltemp_v3.1_' + year + month + day + '.nc' 
elif (product == 'SST-A'):
    # FTP Path
    path = ('pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/ssta/' + year + '/')
    # File Name
    file_name = 'ct5km_ssta_v3.1_' + year + month + day + '.nc' 
elif (product == 'SST-7'):
    # FTP Path
    path = ('pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/sst-trend-7d/' + year + '/')
    # File Name
    file_name = 'ct5km_sst-trend-7d_v3.1_' + year + month + day + '.nc' 

# Enter the FTP Path
ftp.cwd(path)

#-----------------------------------------------------------------------------------------------------------
      
print('\n---------------------')
print('Downloading FTP File:') 
print('---------------------')
print('Product: ' + product)
print('Date: ' + year + '-' + month + '-' + day)
print('File Name: ' + file_name)

# Download the file
ftp.retrbinary("RETR " + file_name, open(dir + '//' + file_name, 'wb').write)  

# Quit the FPT connection
ftp.quit()

#-----------------------------------------------------------------------------------------------------------

# End the time counter
print('\nTotal Processing Time:', round((t.time() - start_time),2), 'seconds.') 
