"""
Calculate time series of polar vortex for the entire year

Notes
-----
    Author : Zachary Labe
    Date   : 25 June 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
import cmocean
import scipy.stats as sts

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
print('\n' '----Plotting Polar Vortex Time Series- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10']
simuh = 'Past' # Enter simulation time (Current,Past)
letters = [r'Mean',r'A',r'B',r'C']
###############################################################################
if simuh == 'Current':
    simuq = 'Cu'
elif simuh == 'Past':
    simuq = 'Pi'
else:
    print(ValueError('Wrong simulation selected!'))
###############################################################################
###############################################################################
###############################################################################
### Call function for 4d variable data
lat,lon,lev,varfuture = MO.readExperiAll(varnames[0],'Future','surface')
lat,lon,lev,varpast = MO.readExperiAll(varnames[0],simuh,'surface')

### Create 2d array of latitude and longitude
lon2,lat2 = np.meshgrid(lon,lat)

### Remove missing data 
varfuture[np.where(varfuture < -1e10)] = np.nan
varpast[np.where(varpast < -1e10)] = np.nan

### Slice ensembles if needed
varfuture = varfuture[:,:,:,:]
varpast = varpast[:,:,:,:]

### Calculate polar vortex strength using 60N
latq = np.where((lat >= 59.5) & (lat <= 60.5))[0]
latu = lat[latq].squeeze()
varfutureu = varfuture[:,:,latq,:].squeeze()
varpastu = varpast[:,:,latq,:].squeeze()

### Calculate zonal mean
varfuturez = np.nanmean(varfutureu[:,:,:],axis=2)
varpastz = np.nanmean(varpastu[:,:,:],axis=2)

### Calculate ensemble mean
futurem = np.nanmean(varfuturez,axis=0)
pastm = np.nanmean(varpastz,axis=0)

### Calculate anomaly
anom = varfuturez - varpastz
anomm = futurem - pastm

### Calculate 2 sigma
stdanom = np.nanstd(anom,axis=0)

### Rearrange time series
anomt = np.append(anom[:,8:],anom[:,:8],axis=1)
anommt = np.append(anomm[8:],anomm[:8],axis=0)
stdanomt = np.append(stdanom[8:],stdanom[:8],axis=0)

###############################################################################
###############################################################################
###############################################################################    
### Plot time series
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})

### Adjust axes in time series plots 
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

plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.plot(anommt+stdanomt,linestyle='-',color='dimgrey',
         alpha=0.4,linewidth=0,clip_on=False)
plt.plot(anommt-stdanomt,linestyle='-',color='dimgrey',
         alpha=0.4,linewidth=0,clip_on=False)
plt.plot(anommt,linestyle='-',color=cmocean.cm.balance(0.01),
         linewidth=3,clip_on=False)
plt.fill_between(np.arange(0,12,1),anommt+stdanomt*2.,anommt-stdanomt*2.,
                 color=cmocean.cm.ice(0.7),alpha=0.7,clip_on=False,
                 zorder=1)
plt.fill_between(np.arange(0,12,1),anommt+stdanomt*1.,anommt-stdanomt*1.,
                 color=cmocean.cm.ice(0.4),alpha=0.7,clip_on=False,
                 zorder=2)

plt.yticks(np.arange(-60,61,5),list(map(str,np.arange(-60,61,5))),
           fontsize=9)
xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr',r'May',
           r'Jun',r'Jul',r'Aug']
plt.xticks(np.arange(0,12,1),xlabels,fontsize=9)
plt.ylabel(r'\textbf{U10 [m/s]}',color='dimgrey',fontsize=12)  
plt.ylim([-45,45])
plt.xlim([0,11])

plt.savefig(directoryfigure + 'PolarVortex_TimeSeries_%s.png' % (simuh),
            dpi=300)