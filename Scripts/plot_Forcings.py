"""
Plot forcing fields of sea ice and heat fluxes for November through April
using PAMIP experiments

Notes
-----
    Author : Zachary Labe
    Date   : 28 June 2019
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
print('\n' '----Plotting Monthly Maps for Forcing Fields- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
experiment = 'Current'

######################
def readDataPeriods(varnames,experiment):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'Future','surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,experiment,'surface')
    
    ### Select ensemble mean period 
    varfuture = varfuture[:,:,:,:]
    varpast = varpast[:,:,:,:]
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Rearrange months (N,D,J,F,M,A)
    varfuturem = np.append(varfuture[:,-2:,:,:],varfuture[:,:4,:,:],
                           axis=1)
    varpastm = np.append(varpast[:,-2:,:,:],varpast[:,:4,:,:],axis=1)

    return varfuturem,varpastm,lat,lon,lev
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
shf,shh,lat,lon,lev = readDataPeriods('SHFLX',experiment)
lhf,lhh,lat,lon,lev = readDataPeriods('LHFLX',experiment)
longf,longh,lat,lon,lev = readDataPeriods('FLNS',experiment)

### Calculate anomalies 
anomsh = np.nanmean(shf - shh,axis=0)
anomlh = np.nanmean(lhf - lhh,axis=0)
anomlong = np.nanmean(longf - longh,axis=0)
varq = list(itertools.chain(*[anomsh,anomlh,anomlong]))

### Create variable names 
varnamesn = list(map(str,np.repeat(['SHFLX'],6))) + \
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

fig = plt.figure(figsize=(11,5))

for v in range(len(varnamesn)):
    ax = plt.subplot(3,6,v+1)
    
    ### Upward fluxes are positive
    if varnamesn[v] == 'LHFLX':
        varqq = varq[v] * 1.
    elif varnamesn[v] == 'SHFLX':
        varqq = varq[v] * 1.
    elif varnamesn[v] == 'THFLX':
        varqq = varq[v] * -1.
    elif varnamesn[v] == 'RNET':
        varqq = varq[v] * -1.
    elif varnamesn[v] == 'LWN':
        varqq = varq[v] * 1.
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'LHFlX':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnamesn[v] == 'SHFLX':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnamesn[v] == 'RNET':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnamesn[v] == 'THFLX':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnamesn[v] == 'LWN':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    
    m = Basemap(projection='npstere',boundinglat=51,lon_0=0,resolution='l',
                round =True,area_thresh=10000)
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
    
    if varnamesn[v] == 'SIC':
        cmap = cm.jim_special_16.mpl_colormap     
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'LHFLX':
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
    if any([v == 0,v == 6,v == 12]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=16,rotation=90,ha='center',va='center')
    if v < 6:
        ax.annotate(r'\textbf{%s}' % labelmonths[v],
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=20,color='dimgrey',rotation=0,
                    ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=9)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 5:
        cbar_ax = fig.add_axes([0.92,0.65,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=9,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 11:
        cbar_ax = fig.add_axes([0.92,0.395,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')    
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k') 
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=9,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 17:
        cbar_ax = fig.add_axes([0.92,0.14,0.015,0.2])                   
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'LHFlX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'SHFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')   
        elif varnamesn[v] == 'THFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k')     
        elif varnamesn[v] == 'LWN':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=11,color='k') 
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=9,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
fig.subplots_adjust(wspace=0,hspace=0)
       
plt.savefig(directoryfigure + 'SeaIce-HeatFlux_%s.png' % experiment,dpi=900)
print('Completed: Script done!')