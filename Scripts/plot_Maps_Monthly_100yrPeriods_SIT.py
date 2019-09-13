"""
Plot maps of PAMIP data for each month from November to April using
the ensemble mean (300) for the thickness experiments

Notes
-----
    Author : Zachary Labe
    Date   : 13 September 2019
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
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Thickness/Periods100/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Maps- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','U700','U200','SLP','Z50','Z500','T2M']
varnames = ['U10']

######################
def readDataPeriods(varnames,sliceq):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'SIT_Fu','surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,'SIT_Cu','surface')
    
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
for v in range(len(varnames)):
    diffm,climom,pvalm,lat,lon,lev = readDataPeriods(varnames[v],'Mean')
    diffa,climoa,pvala,lat,lon,lev = readDataPeriods(varnames[v],'A')
    diffb,climob,pvalb,lat,lon,lev = readDataPeriods(varnames[v],'B')
    diffc,climoc,pvalc,lat,lon,lev = readDataPeriods(varnames[v],'C')

    var = list(itertools.chain(*[diffm,diffa,diffb,diffc]))
    climos = list(itertools.chain(*[climom,climoa,climob,climoc]))    
    pvar = list(itertools.chain(*[pvalm,pvala,pvalb,pvalc]))
    
    ### Plot Variables
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'SLP':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnames[v] == 'Z500':
        limit = np.arange(-60,60.1,5)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z50':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,50) 
    elif varnames[v] == 'U10' or varnames[v] == 'U200' or varnames[v] == 'U500' or varnames[v] == 'U700':
        limit = np.arange(-4,4.1,0.25)
        barlim = np.arange(-4,5,1)
    elif varnames[v] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    elif varnames[v] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,1) 
    elif varnames[v] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.2)
        
    lonq,latq = np.meshgrid(lon,lat)
    
    fig = plt.figure()
    for i in range(len(var)):
        ax1 = plt.subplot(4,6,i+1)
    
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        varn, lons_cyclic = addcyclic(var[i], lon)
        varn, lons_cyclic = shiftgrid(180., varn, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvarn,lons_cyclic = addcyclic(pvar[i], lon)
        pvarn,lons_cyclic = shiftgrid(180.,pvarn,lons_cyclic,start=False)
        climoq,lons_cyclic = addcyclic(climos[i], lon)
        climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
                  
        circle = m.drawmapboundary(fill_color='white',
                                   color='dimgrey',linewidth=0.7)
        circle.set_clip_on(False)
        
        if varnames[v] == 'RNET':
            varn = varn * -1. # change sign for upward fluxes as positive
        
        cs = m.contourf(x,y,varn,limit,extend='both')
        cs1 = m.contourf(x,y,pvarn,colors='None',hatches=['....'])
        if varnames[v] == 'Z30': # the interval is 250 m 
            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                            colors='k',linewidths=1.1,zorder=10)
        if varnames[v] == 'RNET':
            m.drawcoastlines(color='darkgray',linewidth=0.3)
            m.fillcontinents(color='dimgrey')
        else:
            m.drawcoastlines(color='dimgrey',linewidth=0.6)
        
        if any([i==0,i==6,i==12,i==18]):
            ax1.tick_params(labelleft='on')       
        else:
            ax1.tick_params(labelleft='off') 
        if i < 18:
            ax1.tick_params(labelbottom='off') 
        if any([i==0,i==6,i==12]):
            ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                            width=2,color='dimgrey')
            ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')  
        else:
            if i < 24 and i != 18:
                ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')
                if i < 18:
                    ax1.tick_params(axis='y',direction='out',which='major',
                                    pad=3,width=0,color='dimgrey')
                    ax1.tick_params(axis='x',direction='out',which='major',
                                    pad=3,width=0,color='dimgrey')    
   
        if varnames[v] == 'T2M':
            cmap = cmocean.cm.balance             
            cs.set_cmap(cmap)   
        elif varnames[v] == 'SLP':
            cmap = cmocean.cm.balance          
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z50':
            cmap = cmocean.cm.balance  
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U200' or varnames[v] == 'U500' or varnames[v] == 'U700':
            cmap = cmocean.cm.balance            
            cs.set_cmap(cmap)  
        elif varnames[v] == 'SWE':
            cmap = cmap = cmocean.cm.tarn
            cs.set_cmap(cmap)
        elif varnames[v] == 'P':
            cmap = cmocean.cm.tarn         
            cs.set_cmap(cmap) 
        elif varnames[v] == 'THICK':
            cmap = cmocean.cm.balance          
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.diff
            cs.set_cmap(cmap)
        elif varnames[v] == 'RNET':
            cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
    
        labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                        fontsize=13,color='dimgrey',rotation=0,
                        ha='center',va='center')
        if i==0:     
            plt.annotate(r'\textbf{Mean}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==6:     
            plt.annotate(r'\textbf{A}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==12:     
            plt.annotate(r'\textbf{B}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==18:     
            plt.annotate(r'\textbf{C}',
                xy=(0, 0),xytext=(-0.3,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
    

    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.02])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z50':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U200' or varnames[v] == 'U500' or varnames[v] == 'U700':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=6)

    plt.subplots_adjust(hspace=0.0,bottom=0.14,top=0.93,wspace=0.0)
    
    plt.savefig(directoryfigure + '%s_MonthlyProfiles_100yr.png' % varnames[v],
                dpi=300)
    print('Completed: Script done!')        