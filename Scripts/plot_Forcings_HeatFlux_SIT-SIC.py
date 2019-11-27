"""
Plot differences in heat fluxes for SIT experiments minus SIC experiments

Notes
-----
    Author : Zachary Labe
    Date   : 27 November 2019
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
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Forcing/'

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
        
        diff = np.nanmean((varfuturem - varpastm),axis=0) - \
                   np.nanmean((varf_sicm - varh_sicm),axis=0)

    return diff,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
rnet,lat,lon = readDataPeriods('RNET')
shflx,lat,lon = readDataPeriods('SHFLX')
lhflx,lat,lon = readDataPeriods('LHFLX')
long,lat,lon = readDataPeriods('FLNS')

varq = list(itertools.chain(*[rnet,shflx,lhflx,long]))

### Create variable names 
varnamesn = list(map(str,np.repeat(['RNET'],6))) + \
            list(map(str,np.repeat(['SHFLX'],6))) + \
            list(map(str,np.repeat(['LHFLX'],6))) + \
            list(map(str,np.repeat(['LWN'],6)))
labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
letters = list(map(chr, range(97, 123)))

###########################################################################
###########################################################################
###########################################################################
### Plot variable data for Nov-Apr
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for v in range(len(varq)):
    ax = plt.subplot(4,6,v+1)
    
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
        limit = np.arange(-8,8.1,0.5)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'SHFLX':
        limit = np.arange(-8,8.1,0.5)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'RNET':
        limit = np.arange(-8,8.1,0.5)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'THFLX':
        limit = np.arange(-8,8.1,0.5)
        barlim = np.arange(-8,9,4) 
    elif varnamesn[v] == 'LWN':
        limit = np.arange(-8,8.1,0.5)
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
    elif varnamesn[v] == 'THFLX':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'LWN':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    if any([v == 0,v == 6,v == 12,v == 18]):
        ax.annotate(r'\textbf{$\Delta$%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=13,rotation=90,ha='center',va='center')
    if v < 6:
        ax.annotate(r'\textbf{%s}' % labelmonths[v],
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=15,color='dimgrey',rotation=0,
                    ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=6)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 5:
        cbar_ax = fig.add_axes([0.92,0.71,0.01,0.15])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')     
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 11:
        cbar_ax = fig.add_axes([0.92,0.52,0.01,0.15])              
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')     
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 17:
        cbar_ax = fig.add_axes([0.92,0.33,0.01,0.15])                 
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')     
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 23:
        cbar_ax = fig.add_axes([0.92,0.14,0.01,0.15])              
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')     
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
fig.subplots_adjust(wspace=0,hspace=0)
       
plt.savefig(directoryfigure + 'SeaIce-HeatFlux_SIT-SIC.png',dpi=900)
print('Completed: Script done!')