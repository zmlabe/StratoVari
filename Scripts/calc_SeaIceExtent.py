"""
Script calculates sea ice extent in every experiment!

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2019
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import palettable.cubehelix as cm
import cmocean 
import itertools

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Forcing/'
directoryoutput = '/home/zlabe/Documents/Research/StratoVari/Data/SIE/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Arctic SIE- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments and functions
months = np.arange(1,12+1,1)

def readSICData(experiment):
    ### Call function for 4d variable data
    lat,lon,lev,sicq = MO.readExperiAll('SIC',experiment,'surface')
    
    ### Calculate ensemble mean
    sicm = np.nanmean(sicq[:,:,:,:],axis=0)
    
    ### Remove missing data
    sicm[np.where(sicm <= -1e10)] = np.nan
    
    ### Mask data 
    directorymask = '/seley/zlabe/simu/masks/'
    filenamemask = 'domain.camocn.1.9x2.5_gx1v6_090403.nc'
    datam = Dataset(directorymask + filenamemask)
    mask = datam.variables['frac'][:]
    datam.close()
    
    ### Set missing data
    sic = sicm * mask
    sic[sic<0] = 0
    sic[sic>100] = 100
    
    ### Slice for Arctic data (Northern Hemisphere)
    latq = np.where(lat >= 0)[0]
    lat = lat[latq]
    sic = sic[:,latq,:]

    return sic,lat,lon

def calcExtent(sic,lat2):
    """
    Calculate sea ice extent from sea ice concentration grids
    """
    ### Extent is a binary 0 or 1 for 15% SIC threshold
    thresh=15
    sic[np.where(sic<thresh)]=np.nan
    sic[np.where(sic>=thresh)]=1
    
    ext = np.zeros((sic.shape[0]))
    valyr = np.zeros((sic.shape))
    for ti in range(ext.shape[0]):
        for i in range(lat.shape[0]):
            for j in range(lon.shape[0]):
                if sic[ti,i,j] == 1.0:
                   ### Area 1.9x2.5 grid cell [58466.1 = (278.30) * (210.083)]
                   valyr[ti,i,j] = 58466.1 * np.cos(np.radians(lat2[i,j]))
        ext[ti] = np.nansum(valyr[ti,:,:])/1e6
        
    return ext
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
sicf,lat,lon = readSICData('Future')
sicpd,lat,lon = readSICData('Current')
sicpi,lat,lon = readSICData('Past')

### Meshgrid lat and lon
lon2,lat2 = np.meshgrid(lon,lat)

### Calculate extent
extf = calcExtent(sicf,lat2)
extpd = calcExtent(sicpd,lat2)
extpi = calcExtent(sicpi,lat2)

print('Completed: Data processed!')
###############################################################################
###############################################################################
###############################################################################
### Create subplots of sea ice anomalies 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
fig = plt.figure()
ax = plt.subplot(111)

plt.plot(extf,linewidth=4,marker='o',markersize=6,
         color=cmocean.cm.balance(0.8),label=r'\textbf{Future [2$\bf{^{\circ}}$C]}',clip_on=False)
plt.plot(extpd,linewidth=4,marker='o',markersize=6,
         color=cmocean.cm.balance(0.2),label=r'\textbf{Present-Day}',clip_on=False)
plt.plot(extpi,linewidth=4,marker='o',markersize=6,
         color='darkslateblue',label=r'\textbf{Pre-Industrial}',clip_on=False)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.yticks(np.arange(0,30,2),list(map(str,np.arange(0,30,2))))
plt.ylim([0,20])
xlabels = [r'JAN',r'FEB',r'MAR',r'APR',r'MAY',r'JUN',
           r'JUL',r'AUG',r'SEP',r'OCT',r'NOV',r'DEC',]
plt.xticks(np.arange(0,30,1),xlabels)
plt.xlim([0,11])

ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.45)

plt.ylabel(r'\textbf{Sea Ice Extent [$\bf{\times}$10$^{6}$\ \textbf{km}$^2$]}',
                     color='k',fontsize=10,labelpad=9)
plt.legend(shadow=False,fontsize=8,loc='upper center',
           bbox_to_anchor=(0.5, 0.1),fancybox=True,frameon=False,ncol=3)

plt.savefig(directoryfigure + 'SIE_Climos.png',dpi=900)

###############################################################################
###############################################################################
###############################################################################
### Save data output
np.savetxt(directoryoutput + 'Arctic_SIE_PAMIP-1.6.txt',extf)
np.savetxt(directoryoutput + 'Arctic_SIE_PAMIP-1.1.txt',extpd)
np.savetxt(directoryoutput + 'Arctic_SIE_PAMIP-1.5.txt',extpi)