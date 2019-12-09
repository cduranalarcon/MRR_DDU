#######################################################################################################
#######################################################################################################
## set of MRR function 
##
##check_RawFiles: reads raw files line by line and remove extrange characteres (frequent problem in Windows)
## 
##raw2snow: Process raw files using Maahn and Kollias method
##
##MRRQlooks: Make figures of Doppler momments
##
##Time_integ: Temporal integration of processed MRR data
##
##By: Claudio Duran-ALarcon IGE
##
#######################################################################################################
#######################################################################################################

# check for the special character '~' that stops the IMProToo routine
def check_RawFiles(filename, file_out):
    import numpy as np
    import os 
    
    if os.path.isfile(file_out) == True: os.remove(file_out)

    dataset = open(filename,"r")
    dataset0 = open(file_out,"w")
    
    n_errors = 0
    
    for i in np.linspace(1,650000,650000):
        line = dataset.readline()
        if np.size(line.split('~'))==1:
            dataset0.write(line)
        if np.size(line.split('~'))>1:
            n_errors = n_errors + 1

    dataset.close()
    dataset0.close()
    
    return n_errors

#Convert Raw MRR data into Doppler moments using IMProToo    
#TRES = Temporal RESolution (60 s by default)
def raw2snow(file_in,file_out, TRES = 60, Descr = "MRR data", author = "C. Duran-Alarcon, IGE - UGA",ncForm="NETCDF3_CLASSIC"): #convert raw files in MRR Doppler moments
    import pylab, os, time, warnings, sys
    path_IMProToo="lib/IMProToo-master/IMProToo/" #See more in https://github.com/maahn/IMProToo
    sys.path.append(path_IMProToo)
    import core as IMProToo # Minimum width of a peak set as 4 bins(instead of 3) to remove clutter at the expense of sensitivity      
    warnings.filterwarnings("ignore") #Ignore IMProToo warmings
    rawData = IMProToo.mrrRawData(file_in) # Read RawData
    processedSpec = IMProToo.MrrZe(rawData) #create the IMProToo object and load rawData
    processedSpec.averageSpectra(TRES)# integrates the data in 60 seconds by default.
    processedSpec.co["ncCreator"] = author
    processedSpec.co["ncDescription"] = Descr
    processedSpec.co["dealiaseSpectrum"] = True 
    processedSpec.rawToSnow() # Converts RawData into Radar moments
    processedSpec.writeNetCDF(file_out,ncForm=ncForm) # Saves the processed data in a nc file    

#Quicklooks of the MRR moments    
def MRRQlooks(file_in, file_out,year, month, day, Ze_ranges = [-15, 20], W_ranges = [-3, 3], SW_ranges = [0, 2],format='png',dpi=300,name_station = '',max_h=3000, height_plot = False):
    import numpy as np
    from netCDF4 import Dataset
    import matplotlib.dates as mdate
    import pylab

    font = {'family'    :   'serif',
            'weight'    :   'normal',
            'size'      :   20}

    pylab.rc('font', **font)        

    ds = Dataset(file_in)

    Ze = ds.variables["Ze"][:]
    W = ds.variables["W"][:]
    SW = ds.variables["spectralWidth"][:]
    Time = ds.variables["time"][:]
    H = ds.variables["height"][:]

    h = np.linspace(np.nanmin(H),np.nanmax(H),31)

    mat = [Ze, W, SW]

    secs = mdate.epoch2num(Time)

    fig = pylab.figure(figsize=(18,13))
    
    n = 3
    for i in range(1,4):
        ax = fig.add_subplot(n, 1, i)
        if i == 1: vmMn = [Ze_ranges[0], Ze_ranges[1],r'$Z_e$']  
        if i == 2: vmMn = [W_ranges[0], W_ranges[1],r'$W$']    
        if i == 3: vmMn = [SW_ranges[0], SW_ranges[1],r'$\sigma$']    

        pylab.pcolor(secs,h,np.transpose(mat[i-1]), vmin = vmMn[0],vmax = vmMn[1])

        if i == 1:
            pylab.colorbar(label = r'$Z_e$'+" [dB"+r'$Z_e$'+"]")    
            pylab.title(name_station+" - MRR, "+str(year).zfill(4)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2))
        if i == 2:
            pylab.colorbar(label = r'$W$'+" [m s"+r'$^{-1}$'+"]")    
        if i == 3: 
            pylab.xlabel('Time [UTC]') 
            pylab.colorbar(label = r'$\sigma$'+" [m s"+r'$^{-1}$'+"]")
        else: pylab.xlabel('')
        pylab.ylabel('Height [m a.g.l.]')

        pylab.text(secs[int(np.size(secs)/round(24*(secs[-1]-secs[0])))],h[3*np.size(h)/4],vmMn[2],size='x-large')
        date_formatter = mdate.DateFormatter('%H')
        ax.xaxis.set_major_formatter(date_formatter)
        if round(24*(secs[-1]-secs[0]))< 24:
            pylab.axis([round(secs[0],1),round(secs[-1],1),0,max_h])
        else:
            pylab.axis([round(secs[0]),round(secs[-1]),0,max_h])
            pylab.xticks(np.linspace(round(secs[0]),round(secs[-1]),13))

    pylab.savefig(file_out+".png",bbox_inches='tight',format='png',dpi=300)
     
    if height_plot == True:
        fig, ax = pylab.subplots(figsize=(18,5))
    
        pylab.pcolor(secs,range(31),np.transpose(H), vmin = np.nanmin(H),vmax = np.nanmax(H))
    
        pylab.title(name_station+" - MRR, "+str(year).zfill(4)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2))
        pylab.xlabel('Time [UTC]') 
        pylab.colorbar(label = 'Range [m a.g.l.]')
        pylab.ylabel('Range [m]')
    
        pylab.text(secs[int(np.size(secs)/round(24*(secs[-1]-secs[0])))],range(31)[3*np.size(h)/4],'Range',size='x-large')
        date_formatter = mdate.DateFormatter('%H')
        ax.xaxis.set_major_formatter(date_formatter)
        if round(24*(secs[-1]-secs[0]))< 24:
            pylab.axis([round(secs[0],1),round(secs[-1],1),0,30])
        else:
            pylab.axis([round(secs[0]),round(secs[-1]),0,30])
            pylab.xticks(np.linspace(round(secs[0]),round(secs[-1]),13))
    
        pylab.savefig(file_out+'_h.png',bbox_inches='tight',format=format,dpi=dpi)
    ds.close()
#Intersection between two lists    
def intersect(a, b):
    return list(set(a) & set(b))

#Temporal integration from 1min to x-h. Default 1-h
def Time_integ(tres = 1.,file_in = None, date_ini_str = None,date_end_str = None, scale = 1,offset = 0):
    import numpy as np
    import datetime, calendar, pylab,sys, time
    from netCDF4 import Dataset
    #tres: Output temporal resolution
    if file_in != None:
        ds = Dataset(file_in) #Oppen Dataset
        
        #Read variables
        Time = ds.variables["time"][:]
        Ze = scale*ds.variables["Ze"][:]+offset
        W = ds.variables["W"][:]
        SW = ds.variables["spectralWidth"][:]
        
        #Masking variables
        Ze=np.ma.masked_where(Ze == -9999., Ze)                
        W=np.ma.masked_where(W == -9999., W)                
        SW=np.ma.masked_where(SW == -9999., SW)   
        
        #Define Automatically the start and end dates of the datasets
        if date_ini_str == None: date_ini_str = str(time.gmtime(Time[0])[0])+str(time.gmtime(Time[0])[1]).zfill(2)+str(time.gmtime(Time[0])[2]).zfill(2)+" 00:00"
        if date_end_str == None: date_end_str = str(time.gmtime(Time[-1]+3600*24)[0])+str(time.gmtime(Time[-1]+3600*24)[1]).zfill(2)+str(time.gmtime(Time[-1]+3600*24)[2]).zfill(2)+" 00:00"

        date_ini_int=[int(date_ini_str[0:4]),int(date_ini_str[4:6]),int(date_ini_str[6:8]),int(date_ini_str[9:11]),int(date_ini_str[12:14])]
        date_end_int=[int(date_end_str[0:4]),int(date_end_str[4:6]),int(date_end_str[6:8]),int(date_end_str[9:11]),int(date_end_str[12:14])]

        d_ini=datetime.datetime(date_ini_int[0],date_ini_int[1],date_ini_int[2],date_ini_int[3],date_ini_int[4])
        d_end=datetime.datetime(date_end_int[0],date_end_int[1],date_end_int[2],date_end_int[3],date_end_int[4])

        date_ini = calendar.timegm(d_ini.timetuple())
        date_end = calendar.timegm(d_end.timetuple())

        nsam=int((date_end-date_ini)/(3600.*tres))

        new_times=np.linspace(date_ini+3600.*tres,date_end,int(nsam)) 

        Ze_xh=np.zeros(shape=[int(nsam),28])-9999.
        Ze_xh_full=np.zeros(shape=[int(nsam),28])-9999.
        W_xh=np.zeros(shape=[int(nsam),28])-9999.
        W_xh_full=np.zeros(shape=[int(nsam),28])-9999.
        SW_xh=np.zeros(shape=[int(nsam),28])-9999.
        SW_xh_full=np.zeros(shape=[int(nsam),28])-9999.
        mat_h=np.zeros(shape=[int(nsam),28])
        virga =np.zeros(shape=[int(nsam)])+2
        tyears =np.zeros(shape=[int(nsam)])
        tmonths =np.zeros(shape=[int(nsam)])

        Ze_xh=np.ma.masked_where(Ze_xh == -9999., Ze_xh)                
        Ze_xh_full=np.ma.masked_where(Ze_xh_full == -9999., Ze_xh_full)                
        W_xh=np.ma.masked_where(W_xh == -9999., W_xh)                
        W_xh_full=np.ma.masked_where(W_xh_full == -9999., W_xh_full)                
        SW_xh=np.ma.masked_where(SW_xh == -9999., SW_xh)
        SW_xh_full=np.ma.masked_where(SW_xh_full == -9999., SW_xh_full)                

        n=0

        for t in new_times:    
            pix = np.where((Time<=t) & (Time>t-3600.*tres))    

            if np.size(pix) > 0: #Before 1

                profZe=10*np.log10(np.ma.sum(10**(Ze[pix[0],:]/10.),axis=0)/(60.*tres))
                profW=np.nanmean(W[pix[0],:],axis=0)#/(60.*tres)
                profSW=np.nanmean(SW[pix[0],:],axis=0)#/(60.*tres)

                tyears[n] = time.gmtime(t)[0]
                tmonths[n] = time.gmtime(t)[1]
                Ze_xh[n,:] = profZe.filled(-9999.)#-profZe[1]
                W_xh[n,:] = profW.filled(-9999.)
                SW_xh[n,:] = profSW.filled(-9999.)
                if Ze_xh[n,0] == -9999.: 
                    virga[n] = 1#Virga    
                else:
                    virga[n] = 0 #Not virga    
                        
                if (profZe.mask[0] == False) & (np.nansum(profZe.mask[0:28]) == 0):
                    Ze_xh_full[n,:] = profZe.filled(-9999.)
                    W_xh_full[n,:] = profW.filled(-9999.)
                    SW_xh_full[n,:] = profSW.filled(-9999.)                   
            n=n+1

        Ze_xh = np.ma.masked_where(Ze_xh == -9999., Ze_xh)                
        Ze_xh_full = np.ma.masked_where(Ze_xh_full == -9999., Ze_xh_full)                
        W_xh = np.ma.masked_where(Ze_xh.mask, W_xh)                
        W_xh_full = np.ma.masked_where(Ze_xh_full.mask, W_xh_full)                
        SW_xh = np.ma.masked_where(Ze_xh.mask, SW_xh)                
        SW_xh_full = np.ma.masked_where(Ze_xh_full.mask, SW_xh_full)
        virga = np.array(virga, dtype = int)                
        virga = np.ma.masked_where(virga == 2, virga)                

        MAT = [new_times,tyears,tmonths, Ze_xh,W_xh,SW_xh,Ze_xh_full,W_xh_full,SW_xh_full,virga]

        return MAT
    else:
        print "###>Warning: Filename not defined!###"
        
#Saving time agregated profiles as nc file 
def save_Time_integ(file_out,tres = 1, mat = None, station=""):
    import numpy as np
    from netCDF4 import Dataset
    if mat != None:
        root_grp = Dataset(file_out, 'w', format="NETCDF3_CLASSIC")
        root_grp.description = station+' MRR moments, '+str(int(tres))+'-h temporal resolution'

        # dimensions
        root_grp.createDimension('time', None)
        root_grp.createDimension('range', 28)

        # variables
        time = root_grp.createVariable('time', 'int', ('time',),fill_value=-9999.)
        years = root_grp.createVariable('years', 'int', ('time',),fill_value=-9999.)
        months = root_grp.createVariable('months', 'int', ('time',),fill_value=-9999.)
        Range = root_grp.createVariable('range', 'int', ('range',),fill_value=-9999.)
        Ze = root_grp.createVariable('Ze', 'f', ('time', 'range',),fill_value=-9999.)
        W = root_grp.createVariable('W', 'f', ('time', 'range',),fill_value=-9999.)
        SW = root_grp.createVariable('spectralWidth', 'f', ('time', 'range',),fill_value=-9999.)
        #Ze_full = root_grp.createVariable('Ze_full', 'f', ('time', 'range',),fill_value=-9999.)
        #W_full = root_grp.createVariable('W_full', 'f', ('time', 'range',),fill_value=-9999.)
        #SW_full = root_grp.createVariable('spectralWidth_full', 'f', ('time', 'range',),fill_value=-9999.)
        virga_mask = root_grp.createVariable('virga_mask', 'int', ('time',),fill_value=-9999)
        S = root_grp.createVariable('S', 'f', ('time', 'range',),fill_value=-9999.)
        
        # data
        Range[:] = range(28) 
        time[:] = mat[0]
        years[:] = mat[1]
        months[:] = mat[2]
        Ze[:] = mat[3].filled(-9999.)
        W[:] = mat[4].filled(-9999.)
        SW[:] = mat[5].filled(-9999.)
        #Ze_full[:] = mat[6].filled(-9999.)
        #W_full[:] = mat[7].filled(-9999.)
        #SW_full[:] = mat[8].filled(-9999.)
        virga_mask[:] = mat[9].filled(-9999)
        S[:] = ((10**(mat[3].filled(-9999.)/10.))/76.)**(1/0.91)

        # Variable Attributes
        Range.description = "Range bins"
        time.description = "Measurement time. Following Meteks and IMProToo convention, the dataset at e.g. 11:00 contains the means between 10:00:00 and 11:59:00 (if delta t = 1h)!"
        years.description = "Number of the year"
        months.description = "Number of the month"
        Ze.description = "Mean equivalent reflectivity factor. In absence of precipitation profiles the value of 0 mm^6/m^3 was considered for the average within the time interval."    
        W.description = "Mean Doppler Velocity."
        SW.description = "Mean Doppler Velocity."
        #Ze_full.description = "Equivalent reflectivity factor, full extended vertical profiles. In absence of precipitation profiles the value of 0 mm^6/m^3 was considered in the average within the time interval."    
        #W_full.description = "Mean Doppler Velocity, full extended vertical profiles"
        #SW_full.description = "Mean spectral width, full extended vertical profiles"
        virga_mask.description = "virga/not virga mask, 0: Surface precip (at lowest level), 1: virga"
        S.description = "Snowfall rate derived from S-Ze relationship in Grazioli et al. (2017, TC)"

        Range.units = "#"
        time.units = "seconds since 1970-01-01"
        years.units = "Number of the year"
        months.units = "Number of the month"
        Ze.units = "dBz"
        W.units = "m/s"
        SW.units = "m/s"
        #Ze_full.units = "dBz"
        #W_full.units = "m/s"
        #SW_full.units = "m/s"
        virga_mask.units = "bool"
        S.units = "mm/h"

        root_grp.close() 
    else:
        print "Data matrix not found!"
   

def densplot(x,y,title="a)",bins = None, vmax = None, Range = None):
    import numpy as np
    import pylab
    
    def fmt(x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        #return r'${} \times 10^{{{}}}$'.format(a, b)
        #print a == "0.0"
        print a,b
        if a != "0.0": 
            return r'$'+str(a)+'\xc2\xb7'.decode('utf8')+' 10^{{'+str(b)+'}}$'
        else:
            return r'$'+str(0)+'$' 

    h=np.histogram2d(x,y,bins=bins,range=Range)
    histo_temp=np.ma.masked_where(h[0] <= 0, h[0])
    extent = np.min(h[1]), np.max(h[1]), np.min(h[2]), np.max(h[2])

    histo_temp=100*histo_temp/(1.*np.nansum(histo_temp))

    if vmax==None: vmax=np.max(histo_temp)

    im=pylab.imshow(np.transpose(np.flip(((histo_temp)),1)),interpolation="nearest",
                    aspect='auto',extent=extent,vmin=0,vmax=vmax)
    
    #im=pylab.pcolor(np.linspace(0,2,bins[0]),np.linspace(0,2.9,bins[1]),
    #                np.transpose(histo_temp),vmin=0,vmax=vmax)
    
    if title != '': pylab.title(title)
    print extent
    return im   