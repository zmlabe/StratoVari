"""
Plot maps of PAMIP data for each month from November to April to calculate
the variability of responses per variability 
Notes
-----
    Author : Zachary Labe
    Date   : 8 July 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import cmocean
import itertools
import palettable.cubehelix as cm

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Maps of STD Responses- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','U300','SLP','Z500','Z30','T1000','T2M','THICK','P']

######################
def readDataPeriods(varnames,sliceq):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'Future','surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,'Past','surface')
    
    ### Select ensemble mean period
    if sliceq == 'Mean':
        varfuture = varfuture[:,:,:,:]
        varpast = varpast[:,:,:,:]
    elif sliceq == 'A':
        varfuture = varfuture[:100,:,:,:]
        varpast = varpast[:100,:,:,:]
    elif sliceq == 'B':
        varfuture = varfuture[100:200,:,:,:]
        varpast = varpast[100:200,:,:,:]
    elif sliceq == 'C':
        varfuture = varfuture[200:,:,:,:]
        varpast = varpast[200:,:,:,:]
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Rearrange months (N,D,J,F,M,A)
    varfuturem = np.append(varfuture[:,-2:,:,:],varfuture[:,:4,:,:],
                           axis=1)
    varpastm = np.append(varpast[:,-2:,:,:],varpast[:,:4,:,:],axis=1)
    
    ### Standard deviation of responses
    stdm = np.nanstd(varfuturem-varpastm,axis=0)
    print(varfuturem.shape)
    
    return stdm,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
for v in range(len(varnames)):
    stdm,lat,lon = readDataPeriods(varnames[v],'Mean')
    stda,lat,lon = readDataPeriods(varnames[v],'A')
    stdb,lat,lon = readDataPeriods(varnames[v],'B')
    stdc,lat,lon = readDataPeriods(varnames[v],'C')

    var = list(itertools.chain(*[stdm,stda,stdb,stdc]))
    
    ### Plot Variables
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T2M' or varnames[v] == 'T1000':
        limit = np.arange(0,8.1,1)
        barlim = np.arange(0,9,4)
    elif varnames[v] == 'SLP':
        limit = np.arange(0,15.1,1)
        barlim = np.arange(0,16,5)
    elif varnames[v] == 'Z500':
        limit = np.arange(0,150.1,10)
        barlim = np.arange(0,151,75) 
    elif varnames[v] == 'Z30':
        limit = np.arange(0,500.1,50)
        barlim = np.arange(0,501,250) 
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        limit = np.arange(0,30.1,2)
        barlim = np.arange(0,31,10)
    elif varnames[v] == 'SWE':
        limit = np.arange(0,10.1,0.5)
        barlim = np.arange(0,11,5)
    elif varnames[v] == 'P':
        limit = np.arange(0,5.1,0.25)
        barlim = np.arange(0,5,5) 
    elif varnames[v] == 'THICK':
        limit = np.arange(0,100.1,5)
        barlim = np.arange(0,101,50)
    elif varnames[v] == 'EGR':
        limit = np.arange(0,0.201,0.02)
        barlim = np.arange(0,0.3,0.2)
        
    lonq,latq = np.meshgrid(lon,lat)
    
    fig = plt.figure()
    for i in range(len(var)):
        ax1 = plt.subplot(4,6,i+1)
    
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        varn, lons_cyclic = addcyclic(var[i], lon)
        varn, lons_cyclic = shiftgrid(180., varn, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
                  
        circle = m.drawmapboundary(fill_color='white',
                                   color='dimgrey',linewidth=0.7)
        circle.set_clip_on(False)
        
        cs = m.contourf(x,y,varn,limit,extend='max')

        m.drawcoastlines(color='dimgrey',linewidth=0.6)
        
        if any([i==0,i==6,i==12,i==18]):
            ax1.tick_params(labelleft='on')       
        else:
            ax1.tick_params(labelleft='off') 
        if i < 18:
            ax1.tick_params(labelbottom='off') 
        if any([i==0,i==6,i==12]):
            ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                            width=2,color='dimgrey')
            ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')  
        else:
            if i < 24 and i != 18:
                ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')
                if i < 18:
                    ax1.tick_params(axis='y',direction='out',which='major',
                                    pad=3,width=0,color='dimgrey')
                    ax1.tick_params(axis='x',direction='out',which='major',
                                    pad=3,width=0,color='dimgrey')    
   
        cmap = cm.cubehelix2_16.mpl_colormap                              
        cs.set_cmap(cmap) 
    
        labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                        fontsize=13,color='dimgrey',rotation=0,
                        ha='center',va='center')
        if i==0:     
            plt.annotate(r'\textbf{Mean}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==6:     
            plt.annotate(r'\textbf{A}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==12:     
            plt.annotate(r'\textbf{B}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==18:     
            plt.annotate(r'\textbf{C}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
    

    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.02])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'T2M' or varnames[v] == 'T1000':
        cbar.set_label(r'\textbf{std. dev. [$^\circ$C]}',
                                 fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{std. dev. [m]}',
                                 fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z30':
        cbar.set_label(r'\textbf{std. dev. [m]}',
                                 fontsize=11,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{std. dev. [hPa]}',
                                 fontsize=11,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        cbar.set_label(r'\textbf{std. dev. [m/s]}',
                                 fontsize=11,color='dimgray')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{std. dev. [mm]}',
                                 fontsize=11,color='dimgray')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{std. dev. [mm/day]}',
                                 fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{std. dev. [m]}',
                                 fontsize=11,color='dimgray') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{std. dev. [1/day]}',
                                 fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=6)

    plt.subplots_adjust(hspace=0.0,bottom=0.14,top=0.93,wspace=0.0)
    
    plt.savefig(directoryfigure + '%s_Maps_STDResponse_100yr.png' % varnames[v],
                dpi=300)
    print('Completed: Script done!')


        