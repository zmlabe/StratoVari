"""
Plot difference in heat flux exchanges between the sea ice concentration
and sea ice thickness experiments

Notes
-----
    Author : Zachary Labe
    Date   : 22 November 2019
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
print('\n' '----Plotting Monthly Maps for RNET exchanges- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Functions
def readDataPeriods(varnames):
    ### Call function for 4d variable data
    if varnames == 'SIT':
        dataf = Dataset('/seley/ypeings/simu/PAMIP-1.10-300yr/monthly/SIT_PAMIP-1.10.nc')
        lat = dataf.variables['TLAT'][:]
        lon = dataf.variables['TLON'][:]
        lev = 'surface'
        varfuture = dataf.variables['hi'][:]
        dataf.close()
        
        datah = Dataset('/seley/ypeings/simu/PAMIP-1.9-300yr/monthly/SIT_PAMIP-1.9.nc')
        varpast = datah.variables['hi'][:]
        datah.close()
        
        ### Select ensemble mean period 
        varfuture = varfuture[:,:,:]
        varpast = varpast[:,:,:]
        
        ### Remove missing data
        varfuture[np.where(varfuture <= -1e10)] = np.nan
        varfuture[np.where(varfuture <= 0)] = np.nan
        varpast[np.where(varpast <= -1e10)] = np.nan   
        varpast[np.where(varpast <= 0)] = np.nan

        ### Starts in July
        varfuturem = varfuture[5:11,:,:]
        varpastm = varpast[5:11,:,:]
        
    else:    
        lat,lon,lev,varfuture = MO.readExperiAll(varnames,'Future','surface')
        lat,lon,lev,varpast = MO.readExperiAll(varnames,'Current','surface')
        
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

def readHeatFlux(vari):
    datatf = Dataset('/seley/zlabe/simu/PAMIP-1.10-300yr/monthly/RNET_1701-2000.nc')
    rnet_sitf = datatf[vari][:-12] * -1
    rnet_sitf = np.reshape(rnet_sitf,(300,12,96,144))
    datatf.close()
    
    datath = Dataset('/seley/zlabe/simu/PAMIP-1.9-300yr/monthly/RNET_1701-2000.nc')
    rnet_sith = datath[vari][:-12] * -1
    rnet_sith = np.reshape(rnet_sith,(300,12,96,144))
    datath.close()
    
    datacf = Dataset('/seley/zlabe/simu/PAMIP_Fu/monthly/RNET_1701-2000.nc')
    rnet_sicf = datacf[vari][:] * -1
    rnet_sicf = np.reshape(rnet_sicf,(300,12,96,144))
    datacf.close()
    
    datach = Dataset('/seley/zlabe/simu/PAMIP_Cu/monthly/RNET_1701-2000.nc')
    rnet_sich = datach[vari][:] * -1
    rnet_sich = np.reshape(rnet_sich,(300,12,96,144))
    datach.close()
    
    ### Rearrange months (N,D,J,F,M,A)
    rnetdiff = np.nanmean((rnet_sitf - rnet_sith) - (rnet_sicf - rnet_sich),axis=0)
    diff = np.append(rnetdiff[-2:,:,:],rnetdiff[:4,:,:],
                           axis=0)
    
    print('Completed: Read heat flux data!')
    return diff
    
############################################################################
############################################################################
############################################################################
#### Read in data
sitf,sith,latice,lonice,levice = readDataPeriods('SIT')
sicf,sich,lat,lon,lev = readDataPeriods('SIC')
rnet = readHeatFlux('RNET')

### Calculate anomalies 
anomsic = np.nanmean(sicf - sich,axis=0)
anomheat = rnet
anomsit = sitf - sith
varq = list(itertools.chain(*[anomsic,anomsit,anomheat]))

### Create variable names 
varnamesn = list(map(str,np.repeat(['$\Delta$SIC'],6))) + \
            list(map(str,np.repeat(['$\Delta$SIT'],6))) + \
            list(map(str,np.repeat(['$\Delta$NETFLX'],6)))
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
    ax = plt.subplot(3,6,v+1)
    
    varqq = varq[v]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'SITF':
        limit = np.arange(0,5.1,0.5)
        barlim = np.arange(0,6,1) 
    elif varnamesn[v] == 'SITH':
        limit = np.arange(0,5.1,0.5)
        barlim = np.arange(0,6,1) 
    elif varnamesn[v] == '$\Delta$SIT':
        limit = np.arange(-2,2.1,0.25)
        barlim = np.arange(-2,3,1)
    elif varnamesn[v] == '$\Delta$SIC':
        limit = np.arange(-50,51,5)
        barlim = np.arange(-50,51,25)
    elif varnamesn[v] == '$\Delta$NETFLX':
        limit = np.arange(-10,10.1,0.25)
        barlim = np.arange(-10,11,5) 
    
    if v > 5 and v < 12:
        m = Basemap(projection='npstere',boundinglat=51,lon_0=0,resolution='l',
                    round =True,area_thresh=10000)
        circle = m.drawmapboundary(fill_color='white',
                                    color='dimgrey',linewidth=0.7)
        circle.set_clip_on(False)    
        cs = m.contourf(lonice,latice,varqq,limit,extend='both',latlon=True)
    else:
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
    
    if varnamesn[v] == '$\Delta$SIC':
        cmap = cmocean.cm.balance_r  
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'SITF':
        cmap = 'cubehelix' 
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'SITH':
        cmap = 'cubehelix' 
        cs.set_cmap(cmap)  
    elif varnamesn[v] == '$\Delta$SIT':
        cmap = cmocean.cm.balance_r   
        cs.set_cmap(cmap)  
    elif varnamesn[v] == '$\Delta$NETFLX':
        cmap = cmocean.cm.balance   
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    if any([v == 0,v == 6,v == 12]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
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
        cbar_ax = fig.add_axes([0.92,0.68,0.01,0.15])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == '$\Delta$SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=7.5,color='k')
        elif varnamesn[v] == '$\Delta$NETFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')  
        else:
            cbar.set_label(r'\textbf{m}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 11:
        cbar_ax = fig.add_axes([0.92,0.43,0.01,0.15])              
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='max',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == '$\Delta$SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=7.5,color='k')
        elif varnamesn[v] == '$\Delta$NETFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')  
        else:
            cbar.set_label(r'\textbf{m}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 17:
        cbar_ax = fig.add_axes([0.92,0.17,0.01,0.15])                 
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == '$\Delta$SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=7.5,color='k')
        elif varnamesn[v] == '$\Delta$NETFLX':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=7.5,color='k')  
        else:
            cbar.set_label(r'\textbf{m}',fontsize=7.5,color='k')
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=5,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
fig.subplots_adjust(wspace=0,hspace=0)
       
plt.savefig(directoryfigure + 'HeatFluxExchanges_Diff_SITSIC.png',dpi=900)
print('Completed: Script done!')