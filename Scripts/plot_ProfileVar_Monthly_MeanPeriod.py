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
for v in range(len(varnames)):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll('%s' % varnames[v],'Future',
                                               'profile')
    lat,lon,lev,varpast = MO.readExperiAll('%s' % varnames[v],'Past',
                                               'profile')
    lat,lon,lev,varcurrent = MO.readExperiAll('%s' % varnames[v],'Current',
                                           'profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    varcurrent[np.where(varcurrent <= -1e10)] = np.nan
    
    ### Rearrange months (N,D,J,F,M,A)
    varfuturem = np.append(varfuture[:,-2:,:,:,:],varfuture[:,:4,:,:,:],
                           axis=1)
    varpastm = np.append(varpast[:,-2:,:,:,:],varpast[:,:4,:,:,:],axis=1)
    varcurrentm = np.append(varcurrent[:,-2:,:,:,:],varcurrent[:,:4,:,:,:],
                            axis=1)
    
    ### Calculate zonal means
    varfuturemz = np.nanmean(varfuturem,axis=4)
    varpastmz = np.nanmean(varpastm,axis=4)
    varcurrentmz = np.nanmean(varcurrentm,axis=4)
    
    ### Calculate anomalies
    anompi = varfuturemz - varpastmz
    anomcu = varfuturemz - varcurrentmz
    
    ### Calculate ensemble mean
    anompim = np.nanmean(anompi,axis=0)
    anomcum = np.nanmean(anomcu,axis=0)
    zdiffruns = np.append(anompim,anomcum,axis=0)
    
    ### Calculate climatologies
    zclimo = np.append(np.nanmean(varpastmz,axis=0),
                           np.nanmean(varcurrentmz,axis=0),axis=0)
    
    ### Calculate significance for each month
    stat_past = np.empty((varpastmz.shape[1],len(lev),len(lat)))
    stat_current = np.empty((varpastmz.shape[1],len(lev),len(lat)))
    pvalue_past= np.empty((varpastmz.shape[1],len(lev),len(lat)))
    pvalue_current = np.empty((varpastmz.shape[1],len(lev),len(lat)))
    for i in range(varpastmz.shape[1]):
        stat_past[i],pvalue_past[i] = UT.calc_indttest(varfuturemz[:,i,:,:],
                                                       varpastmz[:,i,:,:])
        stat_current[i],pvalue_current[i] = UT.calc_indttest(
                                            varfuturemz[:,i,:,:],
                                            varcurrentmz[:,i,:,:])
    
    pruns = np.append(pvalue_past,pvalue_current,axis=0)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot U
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
    for i in range(12):
        ax1 = plt.subplot(2,6,i+1)
    
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
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
                
        cs = plt.contourf(lat,lev,zdiffruns[i],limit,extend='both')
        
        if varnames[v] == 'U': 
            cs2 = plt.contour(lat,lev,zclimo[i],np.arange(-20,101,5),
                              linewidths=0.5,colors='dimgrey')
            
        plt.contourf(latq,levq,pruns[i],colors='None',hatches=['//////'],
                     linewidth=5)   
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=7)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=7)
        plt.minorticks_off()
        
        plt.xlim([0,90])
        plt.ylim([1000,10])
        
        if i==1 or i==2 or i==3 or i==4 or i==5 or i==7 or i==8 or i==9 or i==10 or i==11:
            ax1.tick_params(labelleft='off') 
        if i < 6:
            ax1.tick_params(labelbottom='off') 
            
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
                        xy=(0, 0),xytext=(0.5,1.08),xycoords='axes fraction',
                        fontsize=13,color='dimgrey',rotation=0,
                        ha='center',va='center')
    
    cbar_ax = fig.add_axes([0.312,0.09,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    elif varnames[v] == 'TEMP':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
    elif varnames[v] == 'GEOP':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    elif varnames[v] == 'V':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=8)
          
    plt.annotate(r'\textbf{PAST}',
            xy=(0, 0),xytext=(0.055,0.73),xycoords='figure fraction',
            fontsize=15,color='k',rotation=90,
            ha='center',va='center')        
    plt.annotate(r'\textbf{CURRENT}',
        xy=(0, 0),xytext=(0.055,0.36),xycoords='figure fraction',
        fontsize=15,color='k',rotation=90,
        ha='center',va='center')  
    plt.annotate(r'\textbf{Latitude ($^{\circ}$N)',
        xy=(0, 0),xytext=(0.515,0.15),xycoords='figure fraction',
        fontsize=8,color='k',rotation=0,
        ha='center',va='center')  

    plt.subplots_adjust(hspace=0.1)
    plt.subplots_adjust(bottom=0.21)
    
    plt.savefig(directoryfigure + '%s_MonthlyProfiles.png' % varnames[v],
                dpi=300)
    print('Completed: Script done!')


        