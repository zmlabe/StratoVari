"""
Plot maps of PAMIP data for each month from November to April using
the ensemble mean (300) for signal-to-noise ratio

Notes
-----
    Author : Zachary Labe
    Date   : 1 July 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import cmocean
import itertools

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
print('\n' '----Plotting Monthly SNR- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Read data
def readDataPeriods(varnames):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'SIT_Fu','surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,'SIT_Cu','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Rearrange months (N,D,J,F,M,A)
    varfuturem = np.append(varfuture[:,-2:,:,:],varfuture[:,:4,:,:],
                           axis=1)
    varpastm = np.append(varpast[:,-2:,:,:],varpast[:,:4,:,:],axis=1)
    
    ### Calculate ensemble mean
    future = np.nanmean(varfuturem,axis=0)
    climo = np.nanmean(varpastm,axis=0)
    
    ### Calculate anomalies
    anomall = varfuturem - varpastm
    anomm = future - climo
    
    ### Calculate standard deviation
    anomstd = np.nanstd(anomall,axis=0)
    
    ### Calculate signal to noise (mean/std)
    sig = np.abs(anomm)/anomstd
    
    return sig,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
t2,lat,lon = readDataPeriods('T2M')
slp,lat,lon = readDataPeriods('SLP')
z500,lat,lon = readDataPeriods('Z500')
u300,lat,lon = readDataPeriods('U200')
u10,lat,lon = readDataPeriods('U10')

### Append variables for plotting
var = list(itertools.chain(*[t2,slp,z500,u300,u10]))

### Plot Variables
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(0,1.1,0.1)
barlim = np.round(np.arange(0,1.1,0.2),2)
    
lonq,latq = np.meshgrid(lon,lat)

fig = plt.figure()
for i in range(len(var)):
    ax1 = plt.subplot(5,6,i+1)

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
    
    if any([i==0,i==6,i==12,i==18,i==18]):
        ax1.tick_params(labelleft='on')       
    else:
        ax1.tick_params(labelleft='off') 
    if i < 24:
        ax1.tick_params(labelbottom='off') 
    if any([i==0,i==6,i==12,i==18]):
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey')  
    else:
        if i < 30 and i != 24:
            ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey')
            if i < 24:
                ax1.tick_params(axis='y',direction='out',which='major',
                                pad=3,width=0,color='dimgrey')
                ax1.tick_params(axis='x',direction='out',which='major',
                                pad=3,width=0,color='dimgrey')    
   
    cmap = cmocean.cm.rain
    cs.set_cmap(cmap)  

    labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
    if i < 6:
        ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=13,color='dimgrey',rotation=0,
                    ha='center',va='center')
    if i==0:     
        plt.annotate(r'\textbf{T2M}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
    elif i==6:     
        plt.annotate(r'\textbf{SLP}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
    elif i==12:     
        plt.annotate(r'\textbf{Z500}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
    elif i==18:     
        plt.annotate(r'\textbf{U200}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
    elif i==24:     
        plt.annotate(r'\textbf{U10}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  


cbar_ax = fig.add_axes([0.312,0.09,0.4,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{Signal-to-Noise Ratio}',fontsize=9,
                         color='dimgray',labelpad=0)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.dividers.set_color('dimgrey')
cbar.dividers.set_linewidth(1.2)
cbar.outline.set_linewidth(1.2)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(hspace=0.0,bottom=0.14,top=0.93,wspace=0.0)

plt.savefig(directoryfigure + 'SNR_Variables_Pi.png',dpi=300)
print('Completed: Script done!')


    