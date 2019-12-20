#######################################################################################################
#######################################################################################################
## set of MRR function 
##
##check_RawFiles: reads raw files line by line and remove extrange characteres (frequent problem in Windows)
## 
##raw2snow: Process raw files using Maahn and Kollias method (2012)
##
##MRRQlooks: Make figures of Doppler momments
##
##Time_integ: Temporal integration of processed MRR data
##
##densplot: Compute and plot the joint distribution (e.g. probability of variable vs. height)
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
    import os, time, warnings, sys
    import matplotlib.pyplot as plt
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

#Quicklooks of the MRR moments 1st part.
#Checks if the input data already exist.
#Produce a Quicklook of from a nc file, if the plots does not already exist.  
#Create output directories
def QuickL(year, month, day, path, name_station = 'DDU',TRES = 60, Ze_ranges = [-20, 20],W_ranges = [-3, 3], SW_ranges = [0, 1], S_ranges = [0, 0.5], cmap = 'jet',format = 'png',dpi=300):
    #Import libraries
    import sys, os
    sys.path.append("lib/") # adding lib path
    from MRR_functions import MRRQlooks # Automatic Quicklooks
    import matplotlib.pyplot as plt
    
    file_in = path+"MK_processed/"+str(year)+str(month).zfill(2)+"/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES.nc"
    fig_out = path+"Plots/DopplerMoments_Precip/"+str(year)+str(month).zfill(2)+"/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES"

    if (os.path.isfile(file_in) == True): #PRocess the data only if the input exist

        if (os.path.isfile(fig_out+'.'+format) == False): #Processes the data only once

            if os.path.exists(os.path.dirname(fig_out)) == False: os.mkdir(os.path.dirname(fig_out)) #Creates the output directory
            ##See Quicklooks of the MRR moments 2nd part.    
            MRRQlooks(file_in, fig_out,year, month, day, 
                      Ze_ranges = Ze_ranges, 
                      W_ranges = W_ranges, SW_ranges = SW_ranges,
                      S_ranges = S_ranges,
                      format = format,
                      dpi=dpi, 
                      name_station = name_station,
                      cmap = cmap) 
            plt.show()

        else:
            print "Figure ready "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)+" "+str(int(TRES))+" seconds temporal resolution"
    else:
        print "NetCDF file not found "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2) +" "+str(int(TRES))+" seconds temporal resolution"   

#Quicklooks of the MRR moments 1st part only report.
#Checks if the input data already exist.
#Produce a Quicklook of from a nc file, if the plots does not already exist.  
#Create output directories
def QuickL_report(year, month, day, name_station = 'DDU',TRES = 60, Ze_ranges = [-20, 20],W_ranges = [-3, 3], SW_ranges = [0, 1], S_ranges = [0, 0.5], cmap = 'jet',format = 'png',dpi=300):
    #Import libraries
    import sys, os
    sys.path.append("lib/") # adding lib path
    from MRR_functions import MRRQlooks # Automatic Quicklooks
    import matplotlib.pyplot as plt
    
    file_in = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES_report.nc"
    fig_out = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_"+str(int(TRES))+"TRES"

    if (os.path.isfile(file_in) == True): #PRocess the data only if the input exist

        if (os.path.isfile(fig_out+'.'+format) == False): #Processes the data only once

            if os.path.exists(os.path.dirname(fig_out)) == False: os.mkdir(os.path.dirname(fig_out)) #Creates the output directory
            ##See Quicklooks of the MRR moments 2nd part.    
            MRRQlooks(file_in, fig_out,year, month, day, 
                      Ze_ranges = Ze_ranges, 
                      W_ranges = W_ranges, SW_ranges = SW_ranges,
                      S_ranges = S_ranges,
                      format = format,
                      dpi=dpi, 
                      name_station = name_station,
                      cmap = cmap) 
            plt.show()

        else:
            print "Figure ready "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2)+" "+str(int(TRES))+" seconds temporal resolution"
    else:
        print "NetCDF file not found "+str(year)+"/"+str(month).zfill(2)+"/"+str(day).zfill(2) +" "+str(int(TRES))+" seconds temporal resolution"   
        
    
#Quicklooks of the MRR moments 2nd part.
def MRRQlooks(file_in, file_out,year, month, day, Ze_ranges = [-20, 20], W_ranges = [-3, 3], SW_ranges = [0, 1],S_ranges = [0, 0.5],format='png',dpi=300,name_station = '',max_h=3000, height_plot = False, cmap = 'jet'):
    import numpy as np
    from netCDF4 import Dataset
    import matplotlib.dates as mdate
    import matplotlib.pyplot as plt

    font = {'family'    :   'serif',
            'weight'    :   'normal',
            'size'      :   20}

    plt.rc('font', **font)        

    ds = Dataset(file_in)
    
    Ze = ds.variables["Ze"][:]
    W = ds.variables["W"][:]
    SW = ds.variables["spectralWidth"][:]
    S = ds.variables["S"][:]
    Time = ds.variables["time"][:]
    H = ds.variables["height"][:]

    h = np.linspace(np.nanmin(H),np.nanmax(H),31)

    mat = [Ze, W, SW, S]

    secs = mdate.epoch2num(Time)

    fig = plt.figure(figsize=(18,18))
    
    n = 4
    for i in range(1,5):
        ax = fig.add_subplot(n, 1, i)
        if i == 1: vmMn = [Ze_ranges[0], Ze_ranges[1],r'$Z_e$']  
        if i == 2: vmMn = [W_ranges[0], W_ranges[1],r'$W$']    
        if i == 3: vmMn = [SW_ranges[0], SW_ranges[1],r'$\sigma$']    
        if i == 4: vmMn = [S_ranges[0], S_ranges[1],'Snowfall rate']    
            
        plt.pcolor(secs,h,np.transpose(mat[i-1]), vmin = vmMn[0],vmax = vmMn[1],cmap = cmap)

        if i == 1:
            plt.colorbar(label = r'$Z_e$'+" [dB"+r'$Z_e$'+"]")    
            plt.title(name_station+" - MRR, "+str(year).zfill(4)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2))
        if i == 2:
            plt.colorbar(label = r'$W$'+" [m s"+r'$^{-1}$'+"]")    
        if i == 3: 
            plt.colorbar(label = r'$\sigma$'+" [m s"+r'$^{-1}$'+"]")
        if i == 4:
            plt.xlabel('Time [UTC]') 
            plt.colorbar(label = 'Snowfall rate'+" [mm h"+r'$^{-1}$'+"]")    
            
        else: plt.xlabel('')
        plt.ylabel('Height [m a.g.l.]')

        plt.text(secs[int(np.size(secs)/round(24*(secs[-1]-secs[0])))],h[3*np.size(h)/4],vmMn[2],size='x-large', color = "black")
        date_formatter = mdate.DateFormatter('%H')
        ax.xaxis.set_major_formatter(date_formatter)
        if round(24*(secs[-1]-secs[0]))< 24:
            plt.axis([round(secs[0],1),round(secs[-1],1),0,max_h])
        else:
            plt.axis([round(secs[0]),round(secs[-1]),0,max_h])
            plt.xticks(np.linspace(round(secs[0]),round(secs[-1]),13))

    plt.savefig(file_out+".png",bbox_inches='tight',format='png',dpi=300)
     
    if height_plot == True:
        fig, ax = plt.subplots(figsize=(18,5))
    
        plt.pcolor(secs,range(31),np.transpose(H), vmin = np.nanmin(H),vmax = np.nanmax(H),cmap = cmap)
    
        plt.title(name_station+" - MRR, "+str(year).zfill(4)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2))
        plt.xlabel('Time [UTC]') 
        plt.colorbar(label = 'Range [m a.g.l.]')
        plt.ylabel('Range [m]')
    
        plt.text(secs[int(np.size(secs)/round(24*(secs[-1]-secs[0])))],range(31)[3*np.size(h)/4],'Range',size='x-large')
        date_formatter = mdate.DateFormatter('%H')
        ax.xaxis.set_major_formatter(date_formatter)
        if round(24*(secs[-1]-secs[0]))< 24:
            plt.axis([round(secs[0],1),round(secs[-1],1),0,30])
        else:
            plt.axis([round(secs[0]),round(secs[-1]),0,30])
            plt.xticks(np.linspace(round(secs[0]),round(secs[-1]),13))
    
        plt.savefig(file_out+'_h.png',bbox_inches='tight',format=format,dpi=dpi)
    ds.close()
#Intersection between two lists    
def intersect(a, b):
    return list(set(a) & set(b))

#Temporal integration from 1min to x-h. Default 1-h
def Time_integ(tres = 1.,file_in = None, date_ini_str = None,date_end_str = None, scale = 1,offset = 0):
    import numpy as np
    import datetime, calendar,sys, time
    import matplotlib.pyplot as plt
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
    from netCDF4 import Dataset # Read and write ncCDF files
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

#####
def densplot(var1,var2,title="a)",bins = None,vmin = 0, vmax = None, Range = None, cbar_label='%',xlab = '',ylab = '',cmap = 'jet'):
    import numpy as np
    import matplotlib.pyplot as plt
    from netCDF4 import Dataset # Read and write ncCDF files
    
    #Reshape matriw to vectors
    x = np.reshape(var1, np.size(var1))
    y = np.reshape(var2, np.size(var2))
    
    #Index of no masked data
    pix = np.where(x.mask == False)

    #Only actual data
    x = x[pix]
    y = y[pix]
    
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
    
    #2D histogram
    h=np.histogram2d(x,y,bins=bins,range=Range)
    histo_temp=np.ma.masked_where(h[0] <= 0, h[0])
    
    #Convert counts in perceptage
    histo_temp=100*histo_temp/(1.*np.nansum(histo_temp))

    if vmax==None: vmax=np.max(histo_temp)

    # Plotting 2D histogram    
    im=plt.pcolor(np.linspace(h[1][0],h[1][-1],bins[0]),
                    np.linspace(h[2][0],h[2][-1],bins[1]), 
                    np.transpose(histo_temp),
                    vmin=vmin,vmax=vmax, cmap = cmap)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.colorbar(label = cbar_label)
    plt.title(title)
    
    #Mean vertical profile
    pmean = np.nanmean(var1,axis = 0)
    
    #Percentile profiles    
    p20 = []
    p50 = []
    p80 = []

    for i in range(np.size(var1[1,:])):
    
        pix2 = np.where(var1[:,i].mask == False)
        if np.size(pix2)>0:
            p20.append(np.percentile(var1[:,i][pix2],20))
            p50.append(np.percentile(var1[:,i][pix2],50))
            p80.append(np.percentile(var1[:,i][pix2],80))
                
        else:
            p20.append(-9999.)
            p50.append(-9999.)
            p80.append(-9999.)
    
    # Masking and computing Percentiles
    p20 = np.ma.masked_where(np.array(p20) == -9999,np.array(p20))
    p50 = np.ma.masked_where(np.array(p50) == -9999,np.array(p50))
    p80 = np.ma.masked_where(np.array(p80) == -9999,np.array(p80))
    
    #Plotting vertical profiles
    plt.plot(pmean,var2[0], color = "black", linewidth = 2)
    plt.plot(p20,var2[0], color = "grey", linewidth = 2)
    plt.plot(p50,var2[0],'--', color = "black", linewidth = 2)
    plt.plot(p80,var2[0], color = "grey", linewidth = 2)
    return im   

# Saves Multiple netCDF files
def MFDataset_save(fileout,ds_merged,format="NETCDF3_CLASSIC"):
    from netCDF4 import Dataset # Read and write ncCDF files
    #Create new database
    ds_merged2 = Dataset(fileout, 'w', format=format) 
    ds_merged2.description = ds_merged.description+' Time merged'
    ds_merged2.history = ds_merged.history
    ds_merged2.author = ds_merged.author
    ds_merged2.source = ds_merged.source + "CDA processing"
    ds_merged2.properties = ds_merged.properties

    # dimensions
    ds_merged2.createDimension('time', None)
    ds_merged2.createDimension('range', 31)
    ds_merged2.createDimension('velocity', 192)

    # variables
    time = ds_merged2.createVariable('time', 'int', ('time',),fill_value=-9999.)
    range = ds_merged2.createVariable('range', 'int', ('range',),fill_value=-9999.)
    quality = ds_merged2.createVariable('quality', 'int', ('time', 'range',),fill_value=-9999.)
    etaMask = ds_merged2.createVariable('etaMask', 'int', ('time', 'range','velocity'),fill_value=-9999.)
    velocity = ds_merged2.createVariable('velocity', 'f', ('velocity'),fill_value=-9999.) 
    height = ds_merged2.createVariable('height', 'f', ('time','range'),fill_value=-9999.) 
    eta = ds_merged2.createVariable('eta', 'f', ('time','range','velocity'),fill_value=-9999.)  
    TF = ds_merged2.createVariable('TF', 'f', ('time','range'), fill_value=-9999.) 
    Ze = ds_merged2.createVariable('Ze', 'f', ('time','range'), fill_value=-9999.) 
    spectralWidth = ds_merged2.createVariable('spectralWidth', 'f', ('time','range'), fill_value=-9999.) 
    skewness = ds_merged2.createVariable('skewness', 'f', ('time','range'), fill_value=-9999.) 
    kurtosis = ds_merged2.createVariable('kurtosis', 'f', ('time','range'), fill_value=-9999.) 
    peakVelLeftBorder = ds_merged2.createVariable('peakVelLeftBorder', 'f', ('time', 'range'), fill_value=-9999.) 
    peakVelRightBorder = ds_merged2.createVariable('peakVelRightBorder', 'f', ('time','range'), fill_value=-9999.) 
    leftSlope = ds_merged2.createVariable('leftSlope', 'f', ('time','range'), fill_value=-9999.) 
    rightSlope = ds_merged2.createVariable('rightSlope', 'f', ('time','range'), fill_value=-9999.) 
    W = ds_merged2.createVariable('W', 'f', ('time','range'), fill_value=-9999.) 
    etaNoiseAve = ds_merged2.createVariable('etaNoiseAve', 'f', ('time','range'), fill_value=-9999.) 
    etaNoiseStd = ds_merged2.createVariable('etaNoiseStd', 'f', ('time','range'), fill_value=-9999.) 
    SNR = ds_merged2.createVariable('SNR', 'f', ('time','range'), fill_value=-9999.) 
    S = ds_merged2.createVariable('S', 'f', ('time','range'), fill_value=-9999.) 

    #Input Variables

    time[:] = ds_merged.variables['time'][:]
    range[:] = ds_merged.variables['range'][:]
    quality[:] = ds_merged.variables['quality'][:]
    etaMask[:] = ds_merged.variables['etaMask'][:]
    velocity[:] = ds_merged.variables['velocity'][:]
    height[:] = ds_merged.variables['height'][:]
    eta[:] = ds_merged.variables['eta'][:]
    TF[:] = ds_merged.variables['TF'][:]
    Ze[:] = ds_merged.variables['Ze'][:]
    spectralWidth[:] = ds_merged.variables['spectralWidth'][:]
    skewness[:] = ds_merged.variables['skewness'][:]
    kurtosis[:] = ds_merged.variables['kurtosis'][:]
    peakVelLeftBorder[:] = ds_merged.variables['peakVelLeftBorder'][:]
    peakVelRightBorder[:] = ds_merged.variables['peakVelRightBorder'][:]
    leftSlope[:] = ds_merged.variables['leftSlope'][:]
    rightSlope[:] = ds_merged.variables['rightSlope'][:]
    W[:] = ds_merged.variables['W'][:]
    etaNoiseAve[:] = ds_merged.variables['etaNoiseAve'][:]
    etaNoiseStd[:] = ds_merged.variables['etaNoiseStd'][:]
    SNR[:] = ds_merged.variables['SNR'][:]
    S[:] = ds_merged.variables['S'][:]

    # Variable Attributes

    time.description = ds_merged.variables['time'].description
    range.description = ds_merged.variables['range'].description
    quality.description = ds_merged.variables['quality'].description
    etaMask.description = ds_merged.variables['etaMask'].description
    velocity.description = ds_merged.variables['velocity'].description
    height.description = ds_merged.variables['height'].description
    eta.description = ds_merged.variables['eta'].description
    TF.description = ds_merged.variables['TF'].description
    Ze.description = ds_merged.variables['Ze'].description
    spectralWidth.description = ds_merged.variables['spectralWidth'].description
    skewness.description = ds_merged.variables['skewness'].description
    kurtosis.description = ds_merged.variables['kurtosis'].description
    peakVelLeftBorder.description = ds_merged.variables['peakVelLeftBorder'].description
    peakVelRightBorder.description = ds_merged.variables['peakVelRightBorder'].description
    leftSlope.description = ds_merged.variables['leftSlope'].description
    rightSlope.description = ds_merged.variables['rightSlope'].description
    W.description = ds_merged.variables['W'].description
    etaNoiseAve.description = ds_merged.variables['etaNoiseAve'].description
    etaNoiseStd.description = ds_merged.variables['etaNoiseStd'].description
    SNR.description = ds_merged.variables['SNR'].description
    S.description = ds_merged.variables['S'].description

    time.units = ds_merged.variables['time'].units
    range.units = ds_merged.variables['range'].units
    quality.units = ds_merged.variables['quality'].units
    etaMask.units = ds_merged.variables['etaMask'].units
    velocity.units = ds_merged.variables['velocity'].units
    height.units = ds_merged.variables['height'].units
    eta.units = ds_merged.variables['eta'].units
    TF.units = ds_merged.variables['TF'].units
    Ze.units = ds_merged.variables['Ze'].units
    spectralWidth.units = ds_merged.variables['spectralWidth'].units
    skewness.units = ds_merged.variables['skewness'].units
    kurtosis.units = ds_merged.variables['kurtosis'].units
    peakVelLeftBorder.units = ds_merged.variables['peakVelLeftBorder'].units
    peakVelRightBorder.units = ds_merged.variables['peakVelRightBorder'].units
    leftSlope.units = ds_merged.variables['leftSlope'].units
    rightSlope.units = ds_merged.variables['rightSlope'].units
    W.units = ds_merged.variables['W'].units
    etaNoiseAve.units = ds_merged.variables['etaNoiseAve'].units
    etaNoiseStd.units = ds_merged.variables['etaNoiseStd'].units
    SNR.units = ds_merged.variables['SNR'].units
    S.units = ds_merged.variables['S'].units
    time.timezone = ds_merged.variables['time'].timezone    

    ds_merged2.close()

# Saves a ncfile with less variables to send by email
def CopyDataset_simpler(fileout,ds_merged,format="NETCDF3_CLASSIC"):
    from netCDF4 import Dataset # Read and write ncCDF files
    #Create new database
    ds_merged2 = Dataset(fileout, 'w', format=format) 
    ds_merged2.description = ds_merged.description+' Time merged'
    ds_merged2.history = ds_merged.history
    ds_merged2.author = ds_merged.author
    ds_merged2.source = ds_merged.source + "CDA processing"
    ds_merged2.properties = ds_merged.properties

    # dimensions
    ds_merged2.createDimension('time', None)
    ds_merged2.createDimension('range', 31)
    #ds_merged2.createDimension('velocity', 192)

    # variables
    time = ds_merged2.createVariable('time', 'int', ('time',),fill_value=-9999.)
    range = ds_merged2.createVariable('range', 'int', ('range',),fill_value=-9999.)
    quality = ds_merged2.createVariable('quality', 'int', ('time', 'range',),fill_value=-9999.)
    #etaMask = ds_merged2.createVariable('etaMask', 'int', ('time', 'range','velocity'),fill_value=-9999.)
    #velocity = ds_merged2.createVariable('velocity', 'f', ('velocity'),fill_value=-9999.) 
    height = ds_merged2.createVariable('height', 'f', ('time','range'),fill_value=-9999.) 
    #eta = ds_merged2.createVariable('eta', 'f', ('time','range','velocity'),fill_value=-9999.)  
    #TF = ds_merged2.createVariable('TF', 'f', ('time','range'), fill_value=-9999.) 
    Ze = ds_merged2.createVariable('Ze', 'f', ('time','range'), fill_value=-9999.) 
    spectralWidth = ds_merged2.createVariable('spectralWidth', 'f', ('time','range'), fill_value=-9999.) 
    #skewness = ds_merged2.createVariable('skewness', 'f', ('time','range'), fill_value=-9999.) 
    #kurtosis = ds_merged2.createVariable('kurtosis', 'f', ('time','range'), fill_value=-9999.) 
    #peakVelLeftBorder = ds_merged2.createVariable('peakVelLeftBorder', 'f', ('time', 'range'), fill_value=-9999.) 
    #peakVelRightBorder = ds_merged2.createVariable('peakVelRightBorder', 'f', ('time','range'), fill_value=-9999.) 
    #leftSlope = ds_merged2.createVariable('leftSlope', 'f', ('time','range'), fill_value=-9999.) 
    #rightSlope = ds_merged2.createVariable('rightSlope', 'f', ('time','range'), fill_value=-9999.) 
    W = ds_merged2.createVariable('W', 'f', ('time','range'), fill_value=-9999.) 
    #etaNoiseAve = ds_merged2.createVariable('etaNoiseAve', 'f', ('time','range'), fill_value=-9999.) 
    #etaNoiseStd = ds_merged2.createVariable('etaNoiseStd', 'f', ('time','range'), fill_value=-9999.) 
    SNR = ds_merged2.createVariable('SNR', 'f', ('time','range'), fill_value=-9999.) 
    S = ds_merged2.createVariable('S', 'f', ('time','range'), fill_value=-9999.) 

    #Input Variables

    time[:] = ds_merged.variables['time'][:]
    range[:] = ds_merged.variables['range'][:]
    quality[:] = ds_merged.variables['quality'][:]
    #etaMask[:] = ds_merged.variables['etaMask'][:]
    #velocity[:] = ds_merged.variables['velocity'][:]
    height[:] = ds_merged.variables['height'][:]
    #eta[:] = ds_merged.variables['eta'][:]
    #TF[:] = ds_merged.variables['TF'][:]
    Ze[:] = ds_merged.variables['Ze'][:]
    spectralWidth[:] = ds_merged.variables['spectralWidth'][:]
    #skewness[:] = ds_merged.variables['skewness'][:]
    #kurtosis[:] = ds_merged.variables['kurtosis'][:]
    #peakVelLeftBorder[:] = ds_merged.variables['peakVelLeftBorder'][:]
    #peakVelRightBorder[:] = ds_merged.variables['peakVelRightBorder'][:]
    #leftSlope[:] = ds_merged.variables['leftSlope'][:]
    #rightSlope[:] = ds_merged.variables['rightSlope'][:]
    W[:] = ds_merged.variables['W'][:]
    #etaNoiseAve[:] = ds_merged.variables['etaNoiseAve'][:]
    #etaNoiseStd[:] = ds_merged.variables['etaNoiseStd'][:]
    SNR[:] = ds_merged.variables['SNR'][:]
    S[:] = ds_merged.variables['S'][:]

    # Variable Attributes

    time.description = ds_merged.variables['time'].description
    range.description = ds_merged.variables['range'].description
    quality.description = ds_merged.variables['quality'].description
    #etaMask.description = ds_merged.variables['etaMask'].description
    #velocity.description = ds_merged.variables['velocity'].description
    height.description = ds_merged.variables['height'].description
    #eta.description = ds_merged.variables['eta'].description
    #TF.description = ds_merged.variables['TF'].description
    Ze.description = ds_merged.variables['Ze'].description
    spectralWidth.description = ds_merged.variables['spectralWidth'].description
    #skewness.description = ds_merged.variables['skewness'].description
    #kurtosis.description = ds_merged.variables['kurtosis'].description
    #peakVelLeftBorder.description = ds_merged.variables['peakVelLeftBorder'].description
    #peakVelRightBorder.description = ds_merged.variables['peakVelRightBorder'].description
    #leftSlope.description = ds_merged.variables['leftSlope'].description
    #rightSlope.description = ds_merged.variables['rightSlope'].description
    W.description = ds_merged.variables['W'].description
    #etaNoiseAve.description = ds_merged.variables['etaNoiseAve'].description
    #etaNoiseStd.description = ds_merged.variables['etaNoiseStd'].description
    SNR.description = ds_merged.variables['SNR'].description
    S.description = ds_merged.variables['S'].description

    time.units = ds_merged.variables['time'].units
    range.units = ds_merged.variables['range'].units
    quality.units = ds_merged.variables['quality'].units
    #etaMask.units = ds_merged.variables['etaMask'].units
    #velocity.units = ds_merged.variables['velocity'].units
    height.units = ds_merged.variables['height'].units
    #eta.units = ds_merged.variables['eta'].units
    #TF.units = ds_merged.variables['TF'].units
    Ze.units = ds_merged.variables['Ze'].units
    spectralWidth.units = ds_merged.variables['spectralWidth'].units
    #skewness.units = ds_merged.variables['skewness'].units
    #kurtosis.units = ds_merged.variables['kurtosis'].units
    #peakVelLeftBorder.units = ds_merged.variables['peakVelLeftBorder'].units
    #peakVelRightBorder.units = ds_merged.variables['peakVelRightBorder'].units
    #leftSlope.units = ds_merged.variables['leftSlope'].units
    #rightSlope.units = ds_merged.variables['rightSlope'].units
    W.units = ds_merged.variables['W'].units
    #etaNoiseAve.units = ds_merged.variables['etaNoiseAve'].units
    #etaNoiseStd.units = ds_merged.variables['etaNoiseStd'].units
    SNR.units = ds_merged.variables['SNR'].units
    S.units = ds_merged.variables['S'].units
    time.timezone = ds_merged.variables['time'].timezone    

    ds_merged2.close()
    
    
#Create output directories    
def mkfolders():
    import os

    if os.path.exists(os.path.dirname('Data/DDU/MK_processed/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/MK_processed/')) 
    if os.path.exists(os.path.dirname('Data/DDU/Plots/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/Plots/')) 
    if os.path.exists(os.path.dirname('Data/DDU/temp/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/temp/')) 
    if os.path.exists(os.path.dirname('Data/DDU/Plots/DopplerMoments_Precip/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/Plots/DopplerMoments_Precip/')) 
    if os.path.exists(os.path.dirname('Data/DDU/Plots/Vertical_profiles/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/Plots/Vertical_profiles/')) 
    if os.path.exists(os.path.dirname('Data/DDU/MK_processed/TimeMerged/')) == False: 
        os.mkdir(os.path.dirname('Data/DDU/MK_processed/TimeMerged/')) 
    
#Make a list of day between two dates     
def find_dates_btw(date_ini,date_end):
    import time
    from calendar import timegm
    import numpy as np
    
    utc_ini = time.strptime(date_ini, "%Y-%m-%d")
    epoch_ini = timegm(utc_ini)
    utc_end = time.strptime(date_end, "%Y-%m-%d")
    epoch_end = timegm(utc_end)

    ndays = int((epoch_end-epoch_ini)//(24*3600.) + 1 )

    epochs = np.linspace(epoch_ini,epoch_end,ndays)
    timestamps = [time.strftime("%Y%m%d", time.gmtime(epochs[i])) for i in range(ndays)]
    
    return [epochs,timestamps]


#input (modified by the user)
###################
#name = "2018-01-04 02:40:00" #format "%Y-%m-%d %H%M%S" #UTC time
#name2 = name[0:4]+"-"+name[5:7]+"-"+name[8:10]+"_"+name[11:13]+name[14:16]+name[17:19] 
#filename = "c:/temp/20180104/IMProToo_MRR_0104.nc" #directory of the netcdf file
#e_factor = 0.001 #exaggeration of the curves 
#path_out = "C:/temp/20180104/"+name2+"_MK12.png"
#xmin = -6
#xmax = 12
#ymin = 0 
#ymax = 3.2 
#mMZe = [-20,20]
#mMW = [-2,2]
##################

def plot_spectra(date_central="2017-02-08 02:40:00",tw=3,e_factor=0.01,xmin=-6,xmax=12,ymin=0,ymax=3.2,mMZe=[-20,20],mMW=[-2,2],path='Data/DDU/MK_processed/',TRES=60):
    #tw = temporal window
    #date_central = "%Y-%m-%d %H:%M:%S"
    from netCDF4 import Dataset
    from calendar import timegm
    import pylab, time
    import numpy as np
    
    #Font format
    font = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 15}

    pylab.rc('font', **font)

    ##Date time
    utc_time = time.strptime(date_central, "%Y-%m-%d %H:%M:%S")
    epoch_time = timegm(utc_time)    

    dday = time.strftime("%Y%m%d", time.gmtime(epoch_time))
    
    filename = path + dday[:6]+'/'+'DDU_'+dday+'_'+str(int(TRES))+'TRES.nc'
    
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
    pylab.title("MK12 spectra "+ date_central)
    pylab.axes(axarr[1])
    p1,=pylab.plot(W,H[0]/1000., color='red',label="W")    
    pylab.xlabel("W [m s"+r'$^{-1}$'+"]")
    pylab.axis([mMW[0],mMW[1],ymin,ymax])
    ax=axarr[1].twiny()
    p2,=ax.plot(Ze,H[0]/1000., color='blue',label="Ze")  
    ax.set_xlabel("Ze [dBZe]")  
    pylab.axis([mMZe[0],mMZe[1],ymin,ymax])
    pylab.legend(handles=[p1, p2],frameon=False)
    #pylab.savefig(path_out,format="png",bbox_inches = 'tight', dpi=600)
    pylab.show()

    #for i in range(31):
    #    print i, Ze2.filled(-9999)[i]