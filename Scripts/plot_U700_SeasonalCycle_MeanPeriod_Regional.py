"""
Plot U700 seasonal cycle from PAMIP data over 300 ensemble members for 
three different regions

Notes
-----
    Author : Zachary Labe
    Date   : 5 July 2019
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
print('\n' '----Plotting Seasonal Cycle of U700 Regional- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)
letters = [r'Mean',r'A',r'B',r'C']
labelregion = [r'Atlantic',r'Eurasia',r'Pacific']

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
    
    ### Calculate over region
    lonq1 = np.where((lon >= 0) & (lon <= 120))[0]
    lonq2 = np.where((lon > 120) & (lon <= 240))[0]
    lonq3 = np.where((lon > 240) & (lon <= 360))[0]
    
    varfuturem1 = varfuturem[:,:,:,lonq1]
    varfuturem2 = varfuturem[:,:,:,lonq2]
    varfuturem3 = varfuturem[:,:,:,lonq3]
    varpastm1 = varpastm[:,:,:,lonq1]
    varpastm2 = varpastm[:,:,:,lonq2]
    varpastm3 = varpastm[:,:,:,lonq3]
    
    ### Calculate anomalies
    anom1 = varfuturem1 - varpastm1
    anom2 = varfuturem2 - varpastm2
    anom3 = varfuturem3 - varpastm3
    
    ### Calculate ensemble mean
    anomm1 = np.nanmean(anom1,axis=0)
    pastm1 = np.nanmean(varpastm1,axis=0)
    anomm2 = np.nanmean(anom2,axis=0)
    pastm2 = np.nanmean(varpastm2,axis=0)
    anomm3 = np.nanmean(anom3,axis=0)
    pastm3 = np.nanmean(varpastm3,axis=0)
    
    ### Calculate zonal mean
    anommz1 = np.nanmean(anomm1,axis=2)
    climoz1 = np.nanmean(pastm1,axis=2)
    anommz2 = np.nanmean(anomm2,axis=2)
    climoz2 = np.nanmean(pastm2,axis=2)
    anommz3 = np.nanmean(anomm3,axis=2)
    climoz3 = np.nanmean(pastm3,axis=2)
    
    ### Calculate significance for each month
    stat1 = np.empty((varpastm1.shape[1],len(lat)))
    pvalue1 = np.empty((varpastm1.shape[1],len(lat)))
    for i in range(varpastm1.shape[1]):
        stat1[i],pvalue1[i] = UT.calc_indttest(
                                    np.nanmean(varfuturem1[:,i,:],axis=2),
                                    np.nanmean(varpastm1[:,i,:],axis=2))
    stat2 = np.empty((varpastm2.shape[1],len(lat)))
    pvalue2 = np.empty((varpastm2.shape[1],len(lat)))
    for i in range(varpastm2.shape[1]):
        stat2[i],pvalue2[i] = UT.calc_indttest(
                                    np.nanmean(varfuturem2[:,i,:],axis=2),
                                    np.nanmean(varpastm2[:,i,:],axis=2))
    stat3 = np.empty((varpastm3.shape[1],len(lat)))
    pvalue3 = np.empty((varpastm3.shape[1],len(lat)))
    for i in range(varpastm3.shape[1]):
        stat3[i],pvalue3[i] = UT.calc_indttest(
                                    np.nanmean(varfuturem3[:,i,:],axis=2),
                                    np.nanmean(varpastm3[:,i,:],axis=2))
    
    return lat,lon,anommz1,anommz2,anommz3,climoz1,climoz2,climoz3,pvalue1,pvalue2,pvalue3

### Call functions
latm,lonm,anommz1m,anommz2m,anommz3m,climoz1m,climoz2m,climoz3m,pvalue1m,pvalue2m,pvalue3m = readDataPeriods('U700','Mean')
lata,lona,anommz1a,anommz2a,anommz3a,climoz1a,climoz2a,climoz3a,pvalue1a,pvalue2a,pvalue3a = readDataPeriods('U700','A')
latb,lonb,anommz1b,anommz2b,anommz3b,climoz1b,climoz2b,climoz3b,pvalue1b,pvalue2b,pvalue3b = readDataPeriods('U700','B')
latc,lonc,anommz1c,anommz2c,anommz3c,climoz1c,climoz2c,climoz3c,pvalue1c,pvalue2c,pvalue3c = readDataPeriods('U700','C')

zanom = [anommz1m,anommz2m,anommz3m,
         anommz1a,anommz2a,anommz3a,
         anommz1b,anommz2b,anommz3b,
         anommz1c,anommz2c,anommz3c]
zclimo = [climoz1m,climoz2m,climoz3m,
          climoz1a,climoz2a,climoz3a,
          climoz1b,climoz2b,climoz3b,
          climoz1c,climoz2c,climoz3c]    
zpval = [pvalue1m,pvalue2m,pvalue3m,
         pvalue1a,pvalue2a,pvalue3a,
         pvalue1b,pvalue2b,pvalue3b,
         pvalue1c,pvalue2c,pvalue3c]
latqq = [latm,latm,latm,lata,lata,lata,latb,latb,latb,latc,latc,latc]

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
    lat = latqq[i]
    
    ax1 = plt.subplot(4,3,i+1)
    
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
    elif any([i==1,i==2,i==4,i==5,i==7,i==8]):
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
    elif i == 3 or i == 6:
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
    elif i == 9:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
    elif i == 10 or i == 11:
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
    cs1 = plt.contour(time,lat,climo.transpose(),np.arange(0,81,3),
                      linewidths=0.7,colors='dimgrey') 
    cs2 = plt.contourf(time,lat,pval.transpose(),
                       colors='None',hatches=['///////'])     
    
    cs.set_cmap(cmocean.cm.balance)
    
    xlabels = [r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
    plt.xticks(np.arange(0,6,1),xlabels,fontsize=5)
    plt.yticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=5)
    plt.xlim([0,5])
    plt.ylim([15,90])
       
    if i==0:     
        plt.annotate(r'\textbf{Mean}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                         color='dimgray')  
    elif i==3:     
        plt.annotate(r'\textbf{A}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                         color='dimgray')  
    elif i==6:     
        plt.annotate(r'\textbf{B}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center') 
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                         color='dimgray')  
    elif i==9:     
        plt.annotate(r'\textbf{C}',
            xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')  
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                         color='dimgray')  
    if i < 3:
        ax1.annotate(r'\textbf{%s}' % labelregion[i],
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=13,color='dimgrey',rotation=0,
                    ha='center',va='center')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()
plt.subplots_adjust(bottom=0.18,hspace=0.17,wspace=0.11,right=0.95)

plt.savefig(directoryfigure + 'Daily_Jet_Past_Regional_100yr.png',dpi=300)
print('Completed: Script done!')