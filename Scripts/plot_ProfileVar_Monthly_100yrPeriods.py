"""
Plot vertical plots of PAMIP data for each month from November to April using
the ensemble mean (300)

Notes
-----
    Author : Zachary Labe
    Date   : 26 June 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
import cmocean
import itertools

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Vertical Profiles- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U','GEOP','TEMP','V','EGR']

######################
def readDataPeriods(varnames,sliceq):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'Future','profile')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,'Current','profile')
    
    ### Select ensemble mean period
    if sliceq == 'Mean':
        varfuture = varfuture[:,:,:,:,:]
        varpast = varpast[:,:,:,:,:]
    elif sliceq == 'A':
        varfuture = varfuture[:100,:,:,:,:]
        varpast = varpast[:100,:,:,:,:]
    elif sliceq == 'B':
        varfuture = varfuture[100:200,:,:,:,:]
        varpast = varpast[100:200,:,:,:,:]
    elif sliceq == 'C':
        varfuture = varfuture[200:,:,:,:,:]
        varpast = varpast[200:,:,:,:,:]
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Rearrange months (N,D,J,F,M,A)
    varfuturem = np.append(varfuture[:,-2:,:,:,:],varfuture[:,:4,:,:,:],
                           axis=1)
    varpastm = np.append(varpast[:,-2:,:,:,:],varpast[:,:4,:,:,:],axis=1)
    
    ### Calculate zonal means
    varfuturemz = np.nanmean(varfuturem,axis=4)
    varpastmz = np.nanmean(varpastm,axis=4)
    
    ### Calculate anomalies
    anompi = varfuturemz - varpastmz
    
    ### Calculate ensemble mean
    anompim = np.nanmean(anompi,axis=0)
    zdiffruns = anompim
    
    ### Calculate climatologies
    zclimo = np.nanmean(varpastmz,axis=0)
    
    ### Calculate significance for each month
    stat_past = np.empty((varpastmz.shape[1],len(lev),len(lat)))
    pvalue_past= np.empty((varpastmz.shape[1],len(lev),len(lat)))
    for i in range(varpastmz.shape[1]):
        stat_past[i],pvalue_past[i] = UT.calc_indttest(varfuturemz[:,i,:,:],
                                                       varpastmz[:,i,:,:])
    
    pruns = pvalue_past
    
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

    zdiffruns = list(itertools.chain(*[diffm,diffa,diffb,diffc]))
    zclimo = list(itertools.chain(*[climom,climoa,climob,climoc]))    
    pruns = list(itertools.chain(*[pvalm,pvala,pvalb,pvalc]))
    
    ### Plot Variables
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'U':
        limit = np.arange(-2,2.1,0.1)
        barlim = np.arange(-2,3,1)
    elif varnames[v] == 'TEMP':
        limit = np.arange(-4,4.1,0.2)
        barlim = np.arange(-4,5,1)
    elif varnames[v] == 'GEOP':
        limit = np.arange(-60,61,2)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'V':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.1)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.08,0.081,0.005)
        barlim = np.arange(-0.08,0.09,0.04)
        
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    latq,levq = np.meshgrid(lat,lev)
    
    fig = plt.figure()
    for i in range(len(zdiffruns)):
        ax1 = plt.subplot(4,6,i+1)
    
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
                
        cs = plt.contourf(lat,lev,zdiffruns[i],limit,extend='both')
        
        if varnames[v] == 'U': 
            cs2 = plt.contour(lat,lev,zclimo[i],np.arange(-20,101,5),
                              linewidths=0.5,colors='dimgrey')
            
        plt.contourf(latq,levq,pruns[i],colors='None',hatches=['//////'],
                     linewidth=5)   
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=5)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=5)
        plt.minorticks_off()
        
        plt.xlim([0,90])
        plt.ylim([1000,10])
        
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
   
        if varnames[v] == 'U':
            cmap = cmocean.cm.balance            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'TEMP':
            cmap = cmocean.cm.balance            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'GEOP':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'V':
            cmap = cmocean.cm.balance             
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.diff           
            cs.set_cmap(cmap) 
    
        labelmonths = [r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
        if i < 6:
            ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                        xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                        fontsize=13,color='dimgrey',rotation=0,
                        ha='center',va='center')
        if i==0:     
            plt.annotate(r'\textbf{Mean}',
                xy=(0, 0),xytext=(-0.6,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==6:     
            plt.annotate(r'\textbf{A}',
                xy=(0, 0),xytext=(-0.6,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==12:     
            plt.annotate(r'\textbf{B}',
                xy=(0, 0),xytext=(-0.6,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
        elif i==18:     
            plt.annotate(r'\textbf{C}',
                xy=(0, 0),xytext=(-0.6,0.5),xycoords='axes fraction',
                fontsize=15,color='k',rotation=90,
                ha='center',va='center')  
    
    cbar_ax = fig.add_axes([0.312,0.07,0.4,0.02])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{m/s}',fontsize=9,color='dimgray',
                       labelpad=0)
    elif varnames[v] == 'TEMP':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=9,color='dimgray',
                       labelpad=0)
    elif varnames[v] == 'GEOP':
        cbar.set_label(r'\textbf{m}',fontsize=9,color='dimgray',
                       labelpad=0)
    elif varnames[v] == 'V':
        cbar.set_label(r'\textbf{m/s}',fontsize=9,color='dimgray',
                       labelpad=0)
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=9,color='dimgray',
                       labelpad=0)
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=6)
        
    plt.annotate(r'\textbf{Latitude ($^{\circ}$N)',
        xy=(0, 0),xytext=(0.515,0.12),xycoords='figure fraction',
        fontsize=6,color='k',rotation=0,
        ha='center',va='center')  

    plt.subplots_adjust(hspace=0.1,bottom=0.17,top=0.93,wspace=0.1)
    
    plt.savefig(directoryfigure + '%s_MonthlyProfiles_100yr.png' % varnames[v],
                dpi=300)
    print('Completed: Script done!')


        