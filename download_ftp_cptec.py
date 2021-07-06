#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: NWP Data Processing With Python - NWP Download with Python (CPTEC FTP - WRF, ETA, BAM)
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------
# Required modules
from datetime import datetime, timedelta # Basic Dates and time types
import os                                # Miscellaneous operating system interfaces
import time as t                         # Time access and conversion                                          
from ftplib import FTP                   # FTP protocol client
#-----------------------------------------------------------------------------------------------------------

print('------------------------------------')
print('CPTEC FTP Download - Script started.')
print('------------------------------------')

# Start the time counter
start_time = t.time()  

# Data description:

# WRF: GRIB2 / ~250 MB / 5 km  / Lats: -55,75 ~ +14,3 / Lons: 276,25 ~ +350,0 / 1h interval + 03 days / Since Feb 08, 2019
# Example FTP link
# http://ftp.cptec.inpe.br/modelos/tempo/WRF/ams_05km/brutos/2021/05/15/00/
# Example file name
# WRF_cpt_05KM_2021051500_2021051500.grib2

# ETA: GRIB1 / ~12 MB  / 40 km / Lats: -50,2 ~ +12,20 / Lons: -83,00001 ~ -25,79 / 1 hour interval + 10 days / Since Jul 08, 2020
# Example FTP link
# http://ftp.cptec.inpe.br/modelos/tempo/Eta/ams_40km/brutos/2021/05/15/00/
# Example file name
# eta_40km_2021051600+2021051600.grb

# BAM: GRIB1 / ~1.5 GB / 20 km / Lats: -90,00001 ~ +90,00001 / Lons: 0.00 ~ +359,82 / 6 h interval + 10 days / Since Jul 01, 2018
# Example FTP link
# http://ftp.cptec.inpe.br/modelos/tempo/BAM/TQ0666L064/brutos/2021/05/15/00/
# Example file name
# GPOSNMC20210516002021052000P.fct.TQ0666L064.grb

#-----------------------------------------------------------------------------------------------------------

# Download directory
dir = "Samples"; os.makedirs(dir, exist_ok=True)

# Select model 
model = 'WRF' # options: 'WRF', 'ETA', 'BAM'

# Desired year (four digit)
year = datetime.today().strftime('%Y') # This will get the current year (four digit)
# Or, select the year you want
#year = '2021' 

# Desired month (two digit)
month = datetime.today().strftime('%m') # This will get the current month (two digit)
# Or, select the month you want
#month = '05'

# Desired day (two digit)
day = datetime.today().strftime('%d') # This will get the current day (two digit)
# Or, select the day you want
#day = '15' 

# Desired run
run = '00'

# Desired forecast hours
hour_ini = 0  # Init time  
hour_end = 0  # End time
hour_int = 1  # Interval

#-----------------------------------------------------------------------------------------------------------

# Create the run date for later use
date_run = datetime.strptime(year + month + day + run, '%Y%m%d%H')
date_run_string = date_run.strftime('%Y%m%d%H')

# FTP Address
ftp = FTP('ftp.cptec.inpe.br') 

# FTP Credentials 
ftp.login('', '') 

# Access the FTP folder, based on the desired model
if (model == 'WRF'):
    # FTP Path
    path = ('modelos/tempo/WRF/ams_05km/brutos/' + year + '/' + month + '/' + day + '/' + run + '/')
elif (model == 'ETA'):
    # FTP Path
    path = ('modelos/tempo/Eta/ams_40km/brutos/' + year + '/' + month + '/' + day + '/' + run + '/')
elif (model == 'BAM'):
    # FTP Path
    path = ('modelos/tempo/BAM/TQ0666L064/brutos/' + year + '/' + month + '/' + day + '/' + run + '/')

# Enter the FTP Path
ftp.cwd(path)

#-----------------------------------------------------------------------------------------------------------

# Download loop
for hour in range(hour_ini, hour_end + 1, hour_int):
    
    # Get the download file name
    date_forecast = (date_run + timedelta(hours=hour)).strftime('%Y%m%d%H')

    # File name creation, based on the desired model
    if (model == 'WRF'):
        # File Name
        file_name = 'WRF_cpt_05KM_' + date_run_string + '_' + date_forecast + '.grib2' 
    elif (model == 'ETA'):
        # File Name
        file_name = 'eta_40km_' + date_run_string + '+' + date_forecast + '.grb'
    elif (model == 'BAM'):
        # File Name
        file_name = 'GPOSNMC' + date_run_string + date_forecast + 'P.fct.TQ0666L064.grb'
    
    print('\n---------------------')
    print('Downloading FTP File:') 
    print('---------------------')
    print('Model: ' + model)
    print('Date: ' + date_run_string)
    print('Run: ' + run)
    print('Forecast Hour: ' + date_forecast)
    print('File Name: ' + file_name)

    # Download the file
    ftp.retrbinary("RETR " + file_name, open(dir + '//' + file_name, 'wb').write)  

# Quit the FPT connection
ftp.quit()

#-----------------------------------------------------------------------------------------------------------

# End the time counter
print('\nTotal Processing Time:', round((t.time() - start_time),2), 'seconds.') 