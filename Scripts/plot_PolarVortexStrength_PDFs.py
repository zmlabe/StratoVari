"""
Calculate PDFs for polar vortex response

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
print('\n' '----Plotting PDF Polar Vortex Subsamples- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10']
period = 'JFM' # Enter temporal period (DJF,JFM,JFMA,ND)
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

### List of experiments
runs = [varfuture,varpast]

### Separate per monthly periods
if period == 'DJF':
    varmo = np.empty((len(runs),varpast.shape[0]-1,varpast.shape[2],
                      varpast.shape[3]))
    for i in range(len(runs)):
        varmo[i,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'surface',17)  
elif period == 'JFM':
    varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                      varpast.shape[3]))
    for i in range(len(runs)):
        varmo[i,:,:,:] = np.nanmean(runs[i][:,:3,:,:],axis=1)
elif period == 'JFMA':
    varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                      varpast.shape[3]))
    for i in range(len(runs)):
        varmo[i,:,:,:] = np.nanmean(runs[i][:,:4,:,:],axis=1)
elif period == 'ND':
    varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                      varpast.shape[3]))
    for i in range(len(runs)):
        varmo[i,:,:,:] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
else:
    ValueError('Wrong period selected! (DJF,JFM,JFMA,ND)')
    
### Remove missing data 
varmo[np.where(varmo < -1e10)] = np.nan

###############################################################################
###############################################################################
###############################################################################
### Slice data for 60N
latq = np.where((lat >= 59.5) & (lat <= 60.5))[0]
latu = lat[latq].squeeze()
varmou = varmo[:,:,latq,:].squeeze()

### Calculate zonal mean
varmoz = np.nanmean(varmou[:,:,:],axis=2)

### Calculate anomalies
anom = varmoz[0,:] - varmoz[1,:]

### Remove nans
mask = ~np.isnan(anom)
anom = anom[mask]

### Fit a distribution
num_bins = np.arange(-50,50,1)
mA,sA = sts.norm.fit(anom[:100])
mB,sB = sts.norm.fit(anom[100:200])
mC,sC = sts.norm.fit(anom[200:])
mm,sm = sts.norm.fit(anom[:])

A = sts.norm.pdf(num_bins,mA,sA)
B = sts.norm.pdf(num_bins,mB,sB)
C = sts.norm.pdf(num_bins,mC,sC)
meann = sts.norm.pdf(num_bins,mm,sm)

plt.figure()
plt.plot(num_bins,A,color='darkblue',linewidth=2.0,label=r'A')
plt.plot(num_bins,B,color='darkgreen',linewidth=2.0,label=r'B')
plt.plot(num_bins,C,color='darkorange',linewidth=2.0,label=r'C')
plt.plot(num_bins,meann,color='k',linewidth=2.0,label=r'Mean',
         linestyle='--',dashes=(1,0.3))

l = plt.legend(shadow=False,fontsize=7,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(0.72,1),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for text in l.get_texts():
    text.set_color('k')
    
### Statistical tests on distribution
tA,pA = sts.ks_2samp(A,meann)
tB,pB = sts.ks_2samp(B,meann)
tC,pC = sts.ks_2samp(C,meann)

print('\n\nP-value between A and mean --> %s!' % np.round(pA,4))
print('P-value between B and mean --> %s!' % np.round(pB,4))
print('P-value between C and mean --> %s!' % np.round(pC,4))
    
plt.savefig(directoryfigure + 'PDFs_PolarVortex_%s_%s.png' % \
            (period,simuh),dpi=300)

