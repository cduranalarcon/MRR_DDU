## Import Libraries
import sys, os
sys.path.append("lib/") # adding lib path
from MRR_functions import check_RawFiles as crf  # Check Raw files for corrupted characters
from MRR_functions import raw2snow # Process raw data into Doppler moments using MK2012
from MRR_functions import QuickL_report #MRR QuickLooks for daily reports
import time
import matplotlib.pyplot as plt #Graphic package
from MRR_functions import mkfolders # create output folders
from netCDF4 import Dataset # Read and write ncCDF files
import numpy as np
from MRR_functions import CopyDataset_simpler

#Define folter station name
Short_name_station = "DDU" 
name_station = '_'.join(Short_name_station.split(' '))

# Data path
path = "Data/"+name_station+"/"

#Graphic parameters
font = {'family'    :   'serif',
        'weight'    :   'normal',
        'size'      :   20}
plt.rc('font', **font)        

#Decription of the output files
Descr = "MRR data at " + name_station + ", first MRR processed with MK12 method v.0.103."

#Creating output folders for daily reports
if os.path.exists(os.path.dirname('Daily-schedules/')) == False: 
    os.mkdir(os.path.dirname('Daily-schedules/')) 

#Creating general output folders
mkfolders()

#Compute previous day
today = time.time()
yesterday = time.strftime("%Y%m%d", time.gmtime(today-3600*24))

#Temporal Resolution
TRES = 60 #Seconds

# Z-S relationship, Z = A*S^B
## Experimental parameters for DDU according to Grazioli et al. (2017)
A = 76
B = 0.91 

#Processing
year=yesterday[:4] 
month=yesterday[4:6] 
day=yesterday[6:8]
file_in = path+"RawSpectra/"+str(year)+str(month).zfill(2)+"/"+str(month).zfill(2)+str(day).zfill(2)+".raw"
file_out = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES0.nc"
file_out2 = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES_report.nc"

temp_file = path+"temp/rawfile.raw"

if (os.path.isfile(file_in) == True): #Processes the data if file exists

    n_errors = crf(file_in, temp_file) #Check for particular characters in the rawfiles that stops the post-processing

    if n_errors > 0: 
        time.sleep(5) #Gives times to the temporal file, ti be created
        print "Processing data of the day "+yesterday+" UTC"
        raw2snow(temp_file,file_out, TRES = TRES, Descr = Descr) # Convert Raw into Doppler Moments using MK2012
        time.sleep(5)  

        ds = Dataset(file_out,'a')
        Ze = ds.variables["Ze"][:]
        Ze = np.ma.masked_where(Ze < -14, Ze) #Min Ze threshold for precipitation at 1minute
        S = ds.createVariable('S', 'f', ('time', 'range',),fill_value=-9999.)
        S[:] = ((10**(Ze/10.))/(1.*A))**(1./B) 
        S.description = "Snowfall rate derived from S-Ze relationship in Grazioli et al. (2017, TC)"
        S.units = "mm/h"
        ds.close()                                                

    else:
        time.sleep(5) #Gives times to the temporal file, ti be created
        print "Processing data of the day "+yesterday+" UTC"
        raw2snow(file_in,file_out, TRES = TRES, Descr = Descr)

        ds = Dataset(file_out,'a')
        Ze = ds.variables["Ze"][:]
        Ze = np.ma.masked_where(Ze < -14, Ze) #Min Ze threshold for precipitation at 1minute
        S = ds.createVariable('S', 'f', ('time', 'range',),fill_value=-9999.)
        S[:] = ((10**(Ze/10.))/(1.*A))**(1./B) 
        S.description = "Snowfall rate derived from S-Ze relationship in Grazioli et al. (2017, TC)"
        S.units = "mm/h"
        ds.close()                        
        time.sleep(5)                              
else:
    print "Raw file not found "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)

ds = Dataset(file_out,'r')
CopyDataset_simpler(file_out2,ds)
ds.close()

QuickL_report(year=yesterday[:4], month=yesterday[4:6], day=yesterday[6:8], name_station = 'DDU',TRES = TRES, Ze_ranges = [-20, 20],W_ranges = [-3, 3], SW_ranges = [0, 1], S_ranges = [0, 0.5], cmap = 'jet',format = 'png',dpi=200)        

os.remove(file_out)