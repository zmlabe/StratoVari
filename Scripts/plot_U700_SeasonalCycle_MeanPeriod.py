"""
Plot U700 seasonal cycle from PAMIP data over 300 ensemble members

Notes
-----
    Author : Zachary Labe
    Date   : 3 July 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
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
print('\n' '----Plotting Seasonal Cycle of U700- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)
letters = [r'Mean',r'A',r'B',r'C']

###############################################################################
###############################################################################
###############################################################################
### Read data
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
    
    ### Calculate anomalies
    anom = varfuturem - varpastm
    
    ### Calculate ensemble mean
    anomm = np.nanmean(anom,axis=0)
    pastm = np.nanmean(varpastm,axis=0)
    
    ### Calculate zonal mean
    anommz = np.nanmean(anomm,axis=2)
    climoz = np.nanmean(pastm,axis=2)
    
    ### Calculate significance for each month
    stat = np.empty((varpastm.shape[1],len(lat)))
    pvalue = np.empty((varpastm.shape[1],len(lat)))
    for i in range(varpastm.shape[1]):
        stat[i],pvalue[i] = UT.calc_indttest(
                                    np.nanmean(varfuturem[:,i,:],axis=2),
                                    np.nanmean(varpastm[:,i,:],axis=2))
    
    return anommz,climoz,lat,lon,pvalue

### Call functions
anomm,climom,lat,lon,pvalm = readDataPeriods('U700','Mean')
anoma,climoa,lat,lon,pvala = readDataPeriods('U700','A')
anomb,climob,lat,lon,pvalb = readDataPeriods('U700','B')
anomc,climoc,lat,lon,pvalc = readDataPeriods('U700','C')

zanom = [anomm,anoma,anomb,anomc]
zclimo = [climom,climoa,climob,climoc]    
zpval = [pvalm,pvala,pvalb,pvalc]

###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
for i in range(len(zanom)):
    anom = zanom[i]
    climo = zclimo[i]
    pval = zpval[i]
    
    ax1 = plt.subplot(2,2,i+1)
    
    if i == 0:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelbottom=False)    
        ax1.yaxis.set_ticks_position('left')
    elif i == 1:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelbottom=False)    
    elif i == 2:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
    elif i == 3:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
    
    ### Set limits
    limit = np.arange(-1,1.01,0.05)
    barlim = np.arange(-1,1.1,0.5)
    time = np.arange(6)
    
    cs = plt.contourf(time,lat,anom.transpose(),limit,extend='both')
    cs1 = plt.contour(time,lat,climo.transpose(),np.arange(0,81,2),
                      linewidths=1,colors='k') 
    cs2 = plt.contourf(time,lat,pval.transpose(),
                       colors='None',hatches=['/////'])     
    
    cs.set_cmap(cmocean.cm.balance)
    
    xlabels = [r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
    plt.xticks(np.arange(0,6,1),xlabels,fontsize=6)
    plt.yticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
    plt.xlim([0,5])
    plt.ylim([0,90])
    
    if i == 0 or i == 2:
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                         color='dimgray')    
    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
            xytext=(0.02,0.88),xycoords='axes fraction',
            color='k',fontsize=11)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()
plt.subplots_adjust(bottom=0.2,hspace=0.05,wspace=0.08,right=0.95)

plt.savefig(directoryfigure + 'Daily_Jet_Past_100yr.png',dpi=300)
print('Completed: Script done!')