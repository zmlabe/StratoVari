"""
Plot maps of PAMIP data for each month from November to April using
the ensemble mean (300) for number of ensemble members to get statistical 
significance

Notes
-----
    Author : Zachary Labe
    Date   : 27 June 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_MinEns_Maps as MENS
import palettable.cubehelix as cm

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
print('\n' '----Plotting Monthly Maps for MENS- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','U30','U50','U300','SLP','Z500','Z200','Z30','T2M','THICK']

######################
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

    return varfuturem,varpastm,lat,lon,lev
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
for v in range(len(varnames)):
    future,climo,lat,lon,lev = readDataPeriods(varnames[v],'Mean')
    nmin = MENS.computeMinEns(future,climo,0.001)
    
    ### Plot Variables
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    limit = np.arange(0,300.1,25)
    barlim = np.arange(0,301,50)
        
    lonq,latq = np.meshgrid(lon,lat)
    
    fig = plt.figure(figsize=(11,5))
    for i in range(len(nmin)):
        ax1 = plt.subplot(1,6,i+1)
    
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        varn, lons_cyclic = addcyclic(nmin[i], lon)
        varn, lons_cyclic = shiftgrid(180., varn, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
                  
        circle = m.drawmapboundary(fill_color='white',
                                   color='dimgrey',linewidth=0.7)
        circle.set_clip_on(False)

        cs = m.contourf(x,y,varn,limit,extend='max')

        m.drawcoastlines(color='dimgrey',linewidth=0.6)
   
        cmap = cm.cubehelix1_16.mpl_colormap            
        cs.set_cmap(cmap)   
    
        labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                        fontsize=20,color='dimgrey',rotation=0,
                        ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.3,0.4,0.04])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=True)
    
    cbar.set_label(r'\textbf{Minimum Ensemble Members}',
                   fontsize=13,color='dimgray')  
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.dividers.set_color('dimgrey')
    cbar.dividers.set_linewidth(1.2)
    cbar.outline.set_linewidth(1.2)
    cbar.ax.tick_params(labelsize=10)

    plt.subplots_adjust(hspace=0.0,bottom=0.14,top=0.93,wspace=0.0)
    
    plt.savefig(directoryfigure + '%s_MENS_Maps.png' % varnames[v],dpi=300)
    print('Completed: Script done!')


        