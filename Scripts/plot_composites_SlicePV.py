"""
Plot maps of PAMIP data for DJF data comparing different simulations. 
Statistical test uses the FDR method with alpha_FDR=0.05. Composites are 
sorted by the strength of the polar vortex response. 

Notes
-----
    Author : Zachary Labe
    Date   : 26 September 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
import calc_SlicePolarVortexYears as PV
import cmocean
import string

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Map Comparison- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','U200','Z500','SLP']
varnamesn = np.repeat(varnames,3)
experi = [r'\textbf{$\bf{\Delta}$Mean}',r'\textbf{$\bf{\Delta+}$PV}',
          r'\textbf{$\bf{\Delta-}$PV}']
letters = list(string.ascii_lowercase)
readallinfo = True
period = 'JFM'

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/SortedPVYears/Comparison/%s/' % period

######################
def readDataPeriods(varnames,simulations,period,stat):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,simulations[0],'surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,simulations[1],'surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data [ensembles,months,lat,lon]
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Slice Polar Vortex Years
    if stat == 'all':
        varfuture = varfuture
        varpast = varpast
        print('\n---Number of NEW ensembles is %s!---\n' % (varfuture.shape[0]))
    else:        
        PVq = PV.polarVortexStats(simulations[0],simulations[1],'U10','surface',
                                  stat,period)
        varfuture = varfuture[PVq,:,:,:]
        varpast = varpast[PVq,:,:,:]
        print('\n---Number of NEW ensembles is %s!---\n' % (varfuture.shape[0]))

    ### Calculate over DJF
    if period == 'OND':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,-3:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,-3:,:,:],axis=1)
    elif period == 'D':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,-1:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,-1:,:,:],axis=1)
    elif period == 'DJF':
        print('Calculating over %s months!' % period)
        runs = [varfuture,varpast]
        var_mo = np.empty((2,varpast.shape[0]-1,varpast.shape[2],varpast.shape[3]))
        for i in range(len(runs)):
            var_mo[i,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,'surface',1) 
        varfuturem = var_mo[0]
        varpastm = var_mo[1]
    elif period == 'JFM':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:3,:,:],axis=1)
    elif period == 'JF':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:2,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:2,:,:],axis=1)
    elif period == 'FMA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:4,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:4,:,:],axis=1)
    elif period == 'FM':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:3,:,:],axis=1)
    elif period == 'J':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:1,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:1,:,:],axis=1)
    elif period == 'F':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:2,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:2,:,:],axis=1)
    elif period == 'M':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,2:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,2:3,:,:],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))
    
    ### Calculate anomalies
    anompi = varfuturem - varpastm
    
    ### Calculate ensemble mean
    anompim = np.nanmean(anompi,axis=0)
    zdiffruns = anompim
    
    ### Calculate climatologies
    zclimo = np.nanmean(varpastm,axis=0)
    
    ### Calculate significance for each month (pick method)
    pruns = UT.calc_FDR_ttest(varfuturem[:,:,:],varpastm[:,:,:],0.05) #FDR
    
    return zdiffruns,zclimo,pruns,lat,lon,lev

###########################################################################
###########################################################################
###########################################################################
### Read in data
if readallinfo == True:
    vari = np.empty((len(varnames),3,96,144)) # [variables,simulations,lat,lon]
    clim = np.empty((len(varnames),3,96,144)) # [variables,simulations,lat,lon]
    pval = np.empty((len(varnames),3,96,144)) # [variables,simulations,lat,lon]
    for v in range(len(varnames)):
        diffp,climop,pp,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['SIT_Fu','SIT_Cu'],
                                                         period,'all')
        diffcu,climocu,pcu,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['SIT_Fu','SIT_Cu'],
                                                         period,'>50')
        diffsit,climosit,psit,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['SIT_Fu','SIT_Cu'],
                                                         period,'<50')
        
        vari[v,:,:,:] = np.asarray([diffp,diffcu,diffsit])
        clim[v,:,:,:] = np.asarray([climop,climocu,climosit])
        pval[v,:,:,:] = np.asarray([pp,pcu,psit])
    
    ### Reshape for subplot [subplots,lat,lon]
    varq = np.reshape(vari,
                     (vari.shape[0]*vari.shape[1],vari.shape[2],vari.shape[3]))
    cliq = np.reshape(clim,
                     (clim.shape[0]*clim.shape[1],clim.shape[2],clim.shape[3]))
    pvaq = np.reshape(pval,
                     (pval.shape[0]*pval.shape[1],pval.shape[2],pval.shape[3]))
    varnamesq = np.tile(varnames,3)
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for JFM
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rcParams['hatch.color'] = 'k'

fig = plt.figure()

for v in range(12):
    ax = plt.subplot(4,3,v+1)
    
    ### Retrieve variables and pvalues
    var = varq[v]
    pvar = pvaq[v]
    climo = cliq[v]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'U200':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    if varnamesn[v] == 'U10':
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,6)
    elif varnamesn[v] == 'Z500':
        limit = np.arange(-50,50.1,5)
        barlim = np.arange(-50,51,25) 
    elif varnamesn[v] == 'SLP':
        limit = np.arange(-4,4.1,0.25)
        barlim = np.arange(-4,5,4)
    elif varnamesn[v] == 'Z30':
        limit = np.arange(-100,100.1,5)
        barlim = np.arange(-100,101,50) 
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
 
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    
    cs = m.contourf(x,y,var,limit,extend='both')
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['.....'])    

    m.drawcoastlines(color='dimgrey',linewidth=0.6)
    
    if varnamesn[v] == 'U10':
        cmap = cmocean.cm.balance          
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'U200':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'Z500':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'SLP':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    if any([v == 0,v == 3,v == 6,v == 9]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=14,rotation=90,ha='center',va='center')
    if any([v == 0,v == 1,v == 2]):
        ax.annotate(r'\textbf{%s}' % experi[v],xy=(0,0),xytext=(0.5,1.12),
                     textcoords='axes fraction',color='dimgrey',
                     fontsize=13,rotation=0,ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=7)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 2:
        cbar_ax = fig.add_axes([0.77,0.72,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U200' or varnamesn[v] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')  
        elif varnamesn[v] == 'SLP':
            cbar.set_label(r'\textbf{hPa}',fontsize=9,color='k',labelpad=-0.4)       
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 5:
        cbar_ax = fig.add_axes([0.77,0.51,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U200' or varnamesn[v] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.4)  
        elif varnamesn[v] == 'SLP':
            cbar.set_label(r'\textbf{hPa}',fontsize=9,color='k',labelpad=-0.4)     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 8:
        cbar_ax = fig.add_axes([0.77,0.29,0.013,0.182])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U200' or varnamesn[v] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')  
        elif varnamesn[v] == 'SLP':
            cbar.set_label(r'\textbf{hPa}',fontsize=9,color='k',labelpad=-0.4)     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 11:
        cbar_ax = fig.add_axes([0.77,0.07,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U200' or varnamesn[v] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500' or varnamesn[v] == 'Z30':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.5)  
        elif varnamesn[v] == 'SLP':
            cbar.set_label(r'\textbf{hPa}',fontsize=9,color='k')        
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)

plt.tight_layout()    
fig.subplots_adjust(wspace=-0.75,hspace=0)
    
plt.savefig(directoryfigure + 'SlicePV_Composites_SIT_Fu-SIT_Cu_%s.png' % period,dpi=300)

print('Completed: Script done!')

