"""
Plot differences in heat fluxes for SIT experiments minus SIC experiments for
our AGU presentation

Notes
-----
    Author : Zachary Labe
    Date   : 5 December 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import palettable.cubehelix as cm
import cmocean 
import itertools
from netCDF4 import Dataset

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Differences in Heat Fluxes (SIT-SIC)- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
def readDataPeriods(varnames):
    ### Call function for 4d variable data
    if varnames == 'RNET':
        datatf = Dataset('/seley/zlabe/simu/PAMIP-1.10-300yr/monthly/RNET_1701-2000.nc')
        rnet_sitf = datatf[varnames][:-12]
        lat = datatf['latitude'][:]
        lon = datatf['longitude'][:]
        rnet_sitf = np.reshape(rnet_sitf,(300,12,96,144))
        datatf.close()
        
        datath = Dataset('/seley/zlabe/simu/PAMIP-1.9-300yr/monthly/RNET_1701-2000.nc')
        rnet_sith = datath[varnames][:-12]
        rnet_sith = np.reshape(rnet_sith,(300,12,96,144))
        datath.close()
        
        datacf = Dataset('/seley/zlabe/simu/PAMIP_Fu/monthly/RNET_1701-2000.nc')
        rnet_sicf = datacf[varnames][:]
        rnet_sicf = np.reshape(rnet_sicf,(300,12,96,144))
        datacf.close()
        
        datach = Dataset('/seley/zlabe/simu/PAMIP_Cu/monthly/RNET_1701-2000.nc')
        rnet_sich = datach[varnames][:]
        rnet_sich = np.reshape(rnet_sich,(300,12,96,144))
        datach.close()
        
        ### Rearrange months (N,D,J,F,M,A)
        rnetdiff = np.nanmean((rnet_sitf - rnet_sith) - (rnet_sicf - rnet_sich),axis=0)
        diff = np.append(rnetdiff[-2:,:,:],rnetdiff[:4,:,:],axis=0)
        diff = np.nanmean(diff[1:5,:,:],axis=0)
        
        #######################################################################
        #######################################################################
        #######################################################################  
    else:    
        lat,lon,lev,varfuture = MO.readExperiAll(varnames,'SIT_Fu','surface')
        lat,lon,lev,varpast = MO.readExperiAll(varnames,'SIT_Cu','surface')
        
        ### Remove missing data
        varfuture[np.where(varfuture <= -1e10)] = np.nan
        varpast[np.where(varpast <= -1e10)] = np.nan
        
        ### Rearrange months (N,D,J,F,M,A)
        varfuturem = np.append(varfuture[:,-2:,:,:],varfuture[:,:4,:,:],
                               axis=1)
        varpastm = np.append(varpast[:,-2:,:,:],varpast[:,:4,:,:],axis=1)
        varfuturem = np.nanmean(varfuturem[:,1:5,:,:],axis=1)
        varpastm = np.nanmean(varpastm[:,1:5,:,:],axis=1)
        #######################################################################
        #######################################################################
        ####################################################################### 
        lat,lon,lev,varf_sic = MO.readExperiAll(varnames,'Future','surface')
        lat,lon,lev,varh_sic = MO.readExperiAll(varnames,'Current','surface')
        
        ### Remove missing data
        varf_sic[np.where(varf_sic <= -1e10)] = np.nan
        varh_sic[np.where(varh_sic <= -1e10)] = np.nan
        
        ### Rearrange months (N,D,J,F,M,A)
        varf_sicm = np.append(varf_sic[:,-2:,:,:],varf_sic[:,:4,:,:],
                               axis=1)
        varh_sicm = np.append(varh_sic [:,-2:,:,:],varh_sic [:,:4,:,:],axis=1)
        varf_sicm = np.nanmean(varf_sicm[:,1:5,:,:],axis=1)
        varh_sicm = np.nanmean(varh_sicm[:,1:5,:,:],axis=1)
        
        diff = np.nanmean((varfuturem - varpastm),axis=0) - \
                   np.nanmean((varf_sicm - varh_sicm),axis=0)

    return diff,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
#rnet,lat,lon = readDataPeriods('RNET')
#shflx,lat,lon = readDataPeriods('SHFLX')
#lhflx,lat,lon = readDataPeriods('LHFLX')
#long,lat,lon = readDataPeriods('FLNS')
#
#varq = np.array([rnet,shflx,lhflx,long])
#
#### Create variable names 
#varnamesn = list(map(str,np.repeat(['RNET'],1))) + \
#            list(map(str,np.repeat(['SHFLX'],1))) + \
#            list(map(str,np.repeat(['LHFLX'],1))) + \
#            list(map(str,np.repeat(['LWN'],1)))
#letters = list(map(chr, range(97, 123)))

###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJFM
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure(figsize=(8,4))

for v in range(len(varq)):
    ax = plt.subplot(1,4,v+1)
    
    ### Upward fluxes are positive
    if varnamesn[v] == 'LHFLX':
        varqq = varq[v] * 1.
    elif varnamesn[v] == 'SHFLX':
        varqq = varq[v] * 1.
    elif varnamesn[v] == 'RNET':
        varqq = varq[v] * -1.
    elif varnamesn[v] == 'LWN':
        varqq = varq[v] * 1.
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'LHFLX':
        limit = np.arange(-8,8.1,0.25)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'SHFLX':
        limit = np.arange(-8,8.1,0.25)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'RNET':
        limit = np.arange(-8,8.1,0.25)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'THFLX':
        limit = np.arange(-8,8.1,0.25)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'LWN':
        limit = np.arange(-8,8.1,0.25)
        barlim = np.arange(-8,9,4) 
    
    m = Basemap(projection='npstere',boundinglat=51,lon_0=0,resolution='l',
        round =True,area_thresh=10000)
    var, lons_cyclic = addcyclic(varqq, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    circle = m.drawmapboundary(fill_color='white',
                                color='dimgrey',linewidth=0.7)
    circle.set_clip_on(False)
    
    var, lons_cyclic = addcyclic(varqq, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    cs = m.contourf(x,y,var,limit,extend='both')

    m.drawcoastlines(color='darkgray',linewidth=0.3)
    m.fillcontinents(color='dimgrey')
     
    if varnamesn[v] == 'LHFLX':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'SHFLX':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'RNET':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'LWN':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    ax.annotate(r'\textbf{$\Delta$%s}' % varnamesn[v],xy=(0,0),xytext=(0.5,1.07),
                 textcoords='axes fraction',color='k',
                 fontsize=18,rotation=0,ha='center',va='center')       
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.85,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=8)
            
    ###########################################################################
    cbar_ax = fig.add_axes([0.30,0.14,0.4,0.03])                  
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac=0.07,drawedges=False)    
    cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=10,color='k')
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(labelsize=8,pad=5,labelcolor='k') 
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs,ha='center')
    cbar.ax.tick_params(axis='x', size=.001)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    
plt.tight_layout()
       
plt.savefig(directoryfigure + 'NDJFM_SeaIce-HeatFlux_SIT-SIC_AGU2.png',dpi=900)
print('Completed: Script done!')