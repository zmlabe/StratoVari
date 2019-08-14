"""
Script calculates the average net surface heat flux over the Arctic following
methods from Peings et al. 2014

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2019
"""

### Import modules
import numpy as np
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
import matplotlib.pyplot as plt
import cmocean

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Forcing/'
directoryoutput = '/home/zlabe/Documents/Research/StratoVari/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating Average RNET - %s----' % titletime)

### Alott time series
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

### Read in SIC data
def readSICData(experiment):
    ### Call function for 4d variable data
    lat,lon,lev,sicq = MO.readExperiAll('SIC',experiment,'surface')
    
    ### Calculate ensemble mean
    sicm = np.nanmean(sicq[:,:,:,:],axis=0)
    
    ### Remove missing data
    sicm[np.where(sicm <= -1e10)] = np.nan
    
    ### Slice for Arctic data (Northern Hemisphere)
    latq = np.where(lat >= 40)[0]
    lat = lat[latq]
    sic = sicm[:,latq,:]

    sic[np.where(sic <= 10)] = 0.    
    sic[np.where(sic > 10)] = 1.
    
    return sic,lat,lon

def readHeatFluxData(experiment):
    ### Call function for 4d variable data
    lat,lon,lev,flxq = MO.readExperiAll('RNET',experiment,'surface')
    
    ### Calculate ensemble mean
    flxm = np.nanmean(flxq[:,:,:,:],axis=0) * -1 # upwards is positive
    
    ### Remove missing data
    flxm[np.where(flxm <= -1e10)] = np.nan
    
    ### Slice for Arctic data (Northern Hemisphere)
    latq = np.where(lat >= 40)[0]
    lat = lat[latq]
    flx = flxm[:,latq,:]
    
    return flx,lat,lon

###########################################################################
###########################################################################
###########################################################################
### Read in data
sicpd,lat,lon = readSICData('Current')
sicpi,lat,lon = readSICData('Past')

flxf,lat,lon = readHeatFluxData('Future')
flxpd,lat,lon = readHeatFluxData('Current')
flxpi,lat,lon = readHeatFluxData('Past')

### Create meshgrid for lat/lon
lon2,lat2 = np.meshgrid(lon,lat)

### Mask where historical sea ice is present to calculate heat flux
heatanompd = (flxf - flxpd) * sicpd
heatanompi = (flxf - flxpi) * sicpi

### Only calculate where heat flux is changing over sea ice areas
heatanompd[np.where(heatanompd==0.)] = np.nan
heatanompi[np.where(heatanompi==0.)] = np.nan

### Calculate weighted average
avepd = UT.calc_weightedAve(heatanompd,lat2)
avepi = UT.calc_weightedAve(heatanompi,lat2)

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

plt.axhline(y=0,linestyle='-',color='dimgrey',linewidth=2)
plt.plot(avepd,linewidth=4,marker='o',markersize=6,
         color=cmocean.cm.balance(0.2),label=r'\textbf{$\bf{\Delta}$Pd}',clip_on=False)
plt.plot(avepi,linewidth=4,marker='o',markersize=6,
         color=cmocean.cm.balance(0.8),label=r'\textbf{$\bf{\Delta}$Pi}',clip_on=False)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.yticks(np.arange(-40,90,10),list(map(str,np.arange(-40,90,10))))
plt.ylim([-10,40])
xlabels = [r'JAN',r'FEB',r'MAR',r'APR',r'MAY',r'JUN',
           r'JUL',r'AUG',r'SEP',r'OCT',r'NOV',r'DEC',]
plt.xticks(np.arange(0,30,1),xlabels)
plt.xlim([0,11])

ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.45)

plt.ylabel(r'\textbf{Net Surface Heat Flux [W/m$\bf{^{2}}$]}',
                     color='k',fontsize=10,labelpad=9)
plt.legend(shadow=False,fontsize=8,loc='upper center',
           bbox_to_anchor=(0.5, 1.01),fancybox=True,frameon=False,ncol=2)

plt.savefig(directoryfigure + 'RNET_Anomalies.png',dpi=900)

###############################################################################
###############################################################################
###############################################################################
### Save data output
np.savetxt(directoryoutput + 'Arctic_RNETA_PI.txt',avepi)
np.savetxt(directoryoutput + 'Arctic_RNETA_PD.txt',avepd)

