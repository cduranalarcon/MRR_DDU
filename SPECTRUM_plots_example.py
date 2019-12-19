# -*- coding: utf-8 -*-
"""
Created on Wed May 02 14:48:18 2018

@author: C. Duran-Alarcon
"""

from netCDF4 import Dataset
from calendar import timegm
import pylab, time
import numpy as np


#input (modified by the user)
###################
name = "2018-01-04 02:40:00" #format "%Y-%m-%d %H%M%S" #UTC time
name2 = name[0:4]+"-"+name[5:7]+"-"+name[8:10]+"_"+name[11:13]+name[14:16]+name[17:19] 
filename = "c:/temp/20180104/IMProToo_MRR_0104.nc" #directory of the netcdf file
e_factor = 0.001 #exaggeration of the curves 
path_out = "C:/temp/20180104/"+name2+"_MK12.png"
xmin = -6
xmax = 12
ymin = 0 
ymax = 3.2 
##################

#Font format
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
        
pylab.rc('font', **font)

##Date time
utc_time = time.strptime(name, "%Y-%m-%d %H:%M:%S")
epoch_time = timegm(utc_time)    

#Open dataset
dataset = Dataset(filename)

#reading variables
##W = dataset.variables["W"][:]
##spectralWidth = dataset.variables["spectralWidth"][:]
##SNR = dataset.variables["SNR"][:]
##quality = dataset.variables["quality"][:]
H = dataset.variables["height"][:]
Time = dataset.variables["time"][:]

#t = Time[1000]

tw=3
pixint=np.where(min(abs(Time-epoch_time))==abs(Time-epoch_time))[0][0]
eta = dataset.variables["eta"][pixint-tw+1:pixint+1,:,:]
etaMask = dataset.variables["etaMask"][pixint-tw+1:pixint+1,:,:]
Ze2 = np.nanmean(dataset.variables["Ze"][pixint-tw+1:pixint+1,:],axis=0)



eta=np.nanmean(eta,axis=0)
etaMask=np.nanmean(etaMask,axis=0)

# conversion into dZe/dv
lamb=299792458./24.15e9 #wavelenght
K2 = 0.92 
eta=10**18*(lamb**4*eta/(np.pi**5*0.92))
dv=0.1893669
vf=dv*(np.linspace(0,191,192)-64)

#masking noise
etaM=np.ma.masked_where(etaMask==1, eta)                
etaM2=np.ma.masked_where((1-etaMask)==1, eta)                

f, axarr = pylab.subplots(1,2, sharey='row',figsize=(12,8))

Ze=[]
Ze.append(-9999)
Ze.append(-9999)
Ze.append(-9999)

#Plotting
pylab.axes(axarr[0])

for ih in range(3,31):
    pylab.plot(vf,(etaM[ih,:]-np.nanmin(etaM2[ih,64:128]))*e_factor+0.1*ih,color="black")
    pylab.plot(vf,(eta[ih,:]-np.nanmin(etaM2[ih,64:128]))*e_factor+0.1*ih,":",color="black")
    #print np.nanmin(etaM[ih,:]-np.nanmean(etaM2[ih,:])),np.nanmax(etaM[ih,:]-np.nanmean(etaM2[ih,:]))
    if (np.nansum(etaM[ih,:]) != np.nan):
        Ze.append(10*np.log10(np.nansum(etaM[ih,:])))
    else:
        Ze.append(-9999)

Ze = np.array(Ze)
Ze = np.ma.masked_where(Ze==-9999, Ze)                
 
W = np.zeros(shape=np.shape(Ze))

nv = 0
for v in vf:
    W = W + etaM[:,nv].filled(0)*v
    nv = nv + 1

W = W/10.**(Ze/10.)

pylab.axis([xmin,xmax,ymin,ymax])
pylab.xlabel("Doppler velocity [m s"+r'$^{-1}$'+"]")
pylab.ylabel("Height [m]")
pylab.title("MK12 spectra "+ name)
pylab.axes(axarr[1])
p1,=pylab.plot(W,H[0]/1000., color='red',label="W")    
pylab.xlabel("W [m s"+r'$^{-1}$'+"]")
pylab.axis([-2,10,ymin,ymax])
ax=axarr[1].twiny()
p2,=ax.plot(Ze,H[0]/1000., color='blue',label="Ze")  
ax.set_xlabel("Ze [dBZe]")  
pylab.axis([10,40,ymin,ymax])
pylab.legend(handles=[p1, p2],frameon=False)
pylab.savefig(path_out,format="png",bbox_inches = 'tight', dpi=600)
pylab.show()

for i in range(31):
    print i, Ze2.filled(-9999)[i]