#MRR processing

#By: Claudio Duran-Alarcon - PhD student, IGE

#Directors: Brice Boudevillain, Alexis Berne

#Project: ANR-APRES3

#Importing libraries

import sys, os
sys.path.append("lib/IMProToo-master/IMProToo/") # adding IMProToo path (see Maahn and Kollias (2012))
sys.path.append("lib/") # adding lib path
#import core as IMProToo
#import numpy as np
from MRR_functions import check_RawFiles as crf 
from MRR_functions import MRRQlooks
from MRR_functions import raw2snow
#from MRR_functions import Time_integ
#from MRR_functions import save_Time_integ
import time#, warnings, pylab, copy
from datetime import datetime
from netCDF4 import Dataset
from shutil import copyfile
#warnings.filterwarnings("ignore") #Ignore IMProToo warmings

#######################
#MODIFIED  BY the USER#
################################################################
Short_name_station = "DDU" #Define folter station name
name_station = '_'.join(Short_name_station.split(' '))
path = "Data/"+name_station+"/"
years =  [2017] # e.g.  [2014];  [2014,2015];  range(2014,2016)   
months = [2] # e.g.  [12];  [4,5,6];  range(1,13)   
days  = range(9,9+1)# e.g.  [1];  [1,2,4];  range(1,5)   

#Description
Descr = "MRR data at " + name_station + ", MRR processed with MK12 method v.0.103."
################################################################

# Process data
if os.path.exists(path+"MK_processed/") == False : os.mkdir(path+"MK_processed/")
if os.path.exists(path+"temp/") == False : os.mkdir(path+"temp/")

for year in years:
    for month in months:
        for day in days: 
            file_in = path+"RawSpectra/"+str(year)+str(month).zfill(2)+"/"+str(month).zfill(2)+str(day).zfill(2)+".raw"
            file_out = path+"MK_processed/"+str(year)+str(month).zfill(2)+"/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+".nc"
            temp_file = path+"temp/rawfile.raw"
            
            if (os.path.isfile(file_in) == True): #Processes the data only once
                
                if (os.path.isfile(file_out) == False):
                    
                    if os.path.exists(os.path.dirname(file_out)) == False: os.mkdir(os.path.dirname(file_out))

                    n_errors = crf(file_in, temp_file) #Check for particular characters in the rawfiles that stops the post-processing

                    if n_errors > 0: 
                        time.sleep(10) #Gives times to the temporal file, ti be created
                        raw2snow(temp_file,file_out, TRES = 60, Descr = Descr)
                    else:
                        raw2snow(file_in,file_out, TRES = 60, Descr = Descr)
                else:
                    print "NetCDF file ready "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)
            else:
                print "Raw file not found "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)
                
                
#Quicklooks

if os.path.exists(path+"Plots/") == False : os.mkdir(path+"Plots/")

for year in years:
    for month in months:
        for day in days: 
            file_in = path+"MK_processed/"+str(year)+str(month).zfill(2)+"/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+".nc"
            fig_out = path+"Plots/"+str(year)+str(month).zfill(2)+"/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)
            
            if (os.path.isfile(file_in) == True): #Processes the data only if file exists
                
                if (os.path.isfile(fig_out+".png") == False): #Processes the data only once
                    
                    if os.path.exists(os.path.dirname(fig_out)) == False: os.mkdir(os.path.dirname(fig_out))
                    			
                    MRRQlooks(file_in, fig_out,year, month, day, Ze_ranges = [-15, 20], W_ranges = [0, 4], SW_ranges = [0, 2], format='png',dpi=300,name_station = Short_name_station)
                    #pylab.show()
                    
                else:
                    print "Figure ready "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)
            else:
                print "NetCDF file not found "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)                    


# Extract subset of data
ds = Dataset(file_out)
times= ds.variables['time']
ranges = ds.variables['range']
height = ds.variables['height']
Ze = ds.variables['Ze']
W = ds.variables['W']
SW = ds.variables['spectralWidth']
SNR = ds.variables['SNR']

NC2transfer = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_sub.nc"

sub_ds = Dataset(NC2transfer, 'w', format="NETCDF3_CLASSIC")
sub_ds.description = ds.description
sub_ds.history = ds.history
sub_ds.author = ds.author
sub_ds.source = ds.source

## dimensions
sub_ds.createDimension('time', None)
sub_ds.createDimension('range', 31)

## variables
times_sub = sub_ds.createVariable('time', 'int', ('time',),fill_value=-9999.)
ranges_sub = sub_ds.createVariable('range', 'int', ('range',),fill_value=-9999.)
height_sub = sub_ds.createVariable('height', 'f', ('time', 'range',),fill_value=-9999.)
Ze_sub = sub_ds.createVariable('Ze', 'f', ('time', 'range',),fill_value=-9999.)
W_sub = sub_ds.createVariable('W', 'f', ('time', 'range',),fill_value=-9999.)
SW_sub = sub_ds.createVariable('spectralWidth', 'f', ('time', 'range',),fill_value=-9999.)
SNR_sub = sub_ds.createVariable('SNR', 'f', ('time', 'range',),fill_value=-9999.)

## data
ranges_sub[:] = ranges[:]
times_sub[:] = times[:]
height_sub[:] = height[:]
Ze_sub[:] = Ze[:]
W_sub[:] = W[:]
SW_sub[:] = SW[:]
SNR_sub[:] = SNR[:]

ranges_sub.description = ranges.description
times_sub.description = times.description
height_sub.description = height.description
Ze_sub.description = Ze.description
W_sub.description = W.description
SW_sub.description = SW.description
SNR_sub.description = SNR.description

sub_ds.close()
ds.close()
copyfile(fig_out+".png", "Daily-schedules/"+os.path.basename(fig_out+".png"))
##Compress file
#os.chdir("Daily-schedules")
#os.system('C:/Users/duran/Documents/PhD/GIT/MRR/lib/WinRAR/RAR.exe  a -m5 -ep -df MRR_data.rar')

'''
#Send email
os.chdir("MRR")
                
sendEmail = '/'.join(str(os.path.abspath('lib/Sendemail/sendEmail.exe')).split('\\'))

#attached1 = '/'.join(str(os.path.abspath('Daily-schedules/MRR_data.rar')).split('\\'))
attached1 = '/'.join(str(os.path.abspath(NC2transfer)).split('\\'))
attached2 = '/'.join(str(os.path.abspath(fig_out+".png")).split('\\'))

FROM = "claudio.duran@univ-grenoble-alpes.fr"
TO = "claudioduran@ug.uchile.cl"
SERVER = "mail.lthe.fr:25"
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y")
MMS = "Daily MRR report from DDU, Antarctica, for the pasts UTC day "+dt_string

os.system(sendEmail+' -a '+attached1+' '+attached2+' -f '+FROM+' -t '+TO+' -u "DAILY MRR REPORT" -s '+SERVER+' -m '+MMS)
os.remove("attached1")
os.remove("attached2")
'''