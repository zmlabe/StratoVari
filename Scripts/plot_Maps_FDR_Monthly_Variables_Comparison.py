"""
Plot maps of PAMIP data for DJF data comparing different simulations. 
Statistical test uses the FDR method with alpha_FDR=0.05

Notes
-----
    Author : Zachary Labe
    Date   : 8 August 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
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
varnames = ['U10','Z50','U200','Z500','SLP','P','T2M','RNET']
experi = np.repeat([r'\textbf{$\bf{\Delta}$Pi}',r'\textbf{$\bf{\Delta}$Cu}',
          r'\textbf{$\bf{\Delta}$SIT}'],len(varnames))
letters = list(string.ascii_lowercase)
readallinfo = True
period = 'JF'

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Comparison/%s_Maps/' % period

######################
def readDataPeriods(varnames,simulations,period):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,simulations[0],'surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,simulations[1],'surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data [ensembles,months,lat,lon]
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Calculate over DJF
    if period == 'OND':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,-3:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,-3:,:,:],axis=1)
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
#    pruns = UT.calc_indttest(varfuturem[:,:,:],varpastm[:,:,:])[1] #ttest
    
    return zdiffruns,zclimo,pruns,lat,lon,lev

###########################################################################
###########################################################################
###########################################################################
### Read in data
if readallinfo == True:
    vari = np.empty((3,len(varnames),96,144)) # [variables,simulations,lat,lon]
    clim = np.empty((3,len(varnames),96,144)) # [variables,simulations,lat,lon]
    pval = np.empty((3,len(varnames),96,144)) # [variables,simulations,lat,lon]
    for v in range(len(varnames)):
        diffp,climop,pp,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['Future','Past'],
                                                         period)
        diffcu,climocu,pcu,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['Future','Current'],
                                                         period)
        diffsit,climosit,psit,lat,lon,lev = readDataPeriods(varnames[v],
                                                         ['SIT_Fu','SIT_Cu'],
                                                         period)
        
        vari[:,v,:,:] = np.asarray([diffp,diffcu,diffsit])
        clim[:,v,:,:] = np.asarray([climop,climocu,climosit])
        pval[:,v,:,:] = np.asarray([pp,pcu,psit])
    
    ### Reshape for subplot [subplots,lat,lon]
    var = np.reshape(vari,
                     (vari.shape[0]*vari.shape[1],vari.shape[2],vari.shape[3]))
    cli = np.reshape(clim,
                     (clim.shape[0]*clim.shape[1],clim.shape[2],clim.shape[3]))
    pva = np.reshape(pval,
                     (pval.shape[0]*pval.shape[1],pval.shape[2],pval.shape[3]))
    varnamesq = np.tile(varnames,3)
 
##########################################################################
##########################################################################
##########################################################################
### Plot settings
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure(figsize=(8,2.5))
for i in range(len(varnamesq)):
    print('Completed: Subplot for %s-%s!' % (varnamesq[i],i+1))
    variablessub = var[i,:,:]
    climossub = cli[i,:,:] 
    pvalues = pva[i,:,:]
    
    ### Plot settings
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnamesq[i] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,10)
    elif varnamesq[i] == 'SLP':
        limit = np.arange(-4,4.1,0.25)
        barlim = np.arange(-4,5,4)
    elif varnamesq[i] == 'Z500':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnamesq[i] == 'Z50':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,100) 
    elif varnamesq[i]=='U10':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnamesq[i]=='U200':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnamesq[i] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    elif varnamesq[i] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,2) 
    elif varnamesq[i] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,30)
    elif varnamesq[i] == 'EGR':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.2)
    elif varnamesq[i] == 'RNET':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50)
        
    ### Meshgrid lat and lon    
    lonq,latq = np.meshgrid(lon,lat)
    
    ax1 = plt.subplot(3,len(varnames),i+1)

    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    varn, lons_cyclic = addcyclic(var[i], lon)
    varn, lons_cyclic = shiftgrid(180., varn, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvarn,lons_cyclic = addcyclic(pvalues, lon)
    pvarn,lons_cyclic = shiftgrid(180.,pvarn,lons_cyclic,start=False)
    climoq,lons_cyclic = addcyclic(climossub, lon)
    climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
              
    circle = m.drawmapboundary(fill_color='white',
                               color='dimgrey',linewidth=0.7)
    circle.set_clip_on(False)
    
    if varnamesq[i] == 'RNET':
        varn = varn * -1. # change sign for upward fluxes as positive
    
    ### Plot contours
    cs = m.contourf(x,y,varn,limit,extend='both')
    cs1 = m.contourf(x,y,pvarn,colors='None',hatches=['.....'])
    if varnames[v] == 'Z50': # the interval is 250 m 
        cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                        colors='k',linewidths=1.1,zorder=10)
        
    ### Set map geography
    if varnamesq[i] == 'RNET':
        m.drawcoastlines(color='darkgrey',linewidth=0.15)
        m.fillcontinents(color='dimgrey')
    else:
        m.drawcoastlines(color='dimgrey',linewidth=0.5)
   
    ### Set colormap
    cs.set_cmap(cmocean.cm.balance)
    
    if i < len(varnames):
        ax1.annotate(r'\textbf{%s}' % varnamesq[i],
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=13,color='k',rotation=0,
                    ha='center',va='center')
    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
        xytext=(0.92,0.9),xycoords='axes fraction',
        color='dimgrey',fontsize=6)
    if any([i==0,i==8,i==16]):     
        plt.annotate(r'%s' % experi[i],
            xy=(0, 0),xytext=(-0.2,0.5),xycoords='axes fraction',
            fontsize=13,color='k',rotation=90,
            ha='center',va='center')  
    
    if i < 16:
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('bottom',size='0%',pad=0.15)   
        cax.axis('off')
    
    if i >= 16:      
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('bottom',size='0%',pad=0.15)   
        cbar = plt.colorbar(cs,ax=cax,orientation='horizontal',
                            extend='both',extendfrac=0.07,drawedges=False,
                            )
        cax.axis('off')
        
        if varnamesq[i] == 'T2M':
            cbar.set_label(r'\textbf{$^\circ$C}',fontsize=7,color='dimgray')  
        elif varnamesq[i] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray')  
        elif varnamesq[i] == 'Z50':
            cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray')  
        elif varnamesq[i] == 'SLP':
            cbar.set_label(r'\textbf{hPa}',fontsize=7,color='dimgray')  
        elif varnamesq[i] == 'U10' or varnamesq[i] == 'U200' or varnamesq[i] == 'U500':
            cbar.set_label(r'\textbf{m/s}',fontsize=7,color='dimgray')  
        elif varnamesq[i] == 'SWE':
            cbar.set_label(r'\textbf{mm}',fontsize=7,color='dimgray')
        elif varnamesq[i] == 'P':
            cbar.set_label(r'\textbf{mm/day}',fontsize=7,color='dimgray') 
        elif varnamesq[i] == 'THICK':
            cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray') 
        elif varnamesq[i] == 'EGR':
            cbar.set_label(r'\textbf{1/day}',fontsize=7,color='dimgray')
        elif varnamesq[i] == 'RNET':
            cbar.set_label(r'\textbf{W/m$\bf{^{2}}$',fontsize=7,color='dimgray')
            
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim))) 
        cbar.ax.tick_params(axis='x', size=.01)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        cbar.ax.tick_params(labelsize=5)

plt.subplots_adjust(hspace=-0.2)

plt.savefig(directoryfigure + 'variable_Comparison_FDR.png',dpi=300)
print('Completed: Script done!')


        