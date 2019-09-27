"""
Plot vertical plots of PAMIP data for each month from January to April using
the ensemble mean RESPONSE by compositing on PV years

Notes
-----
    Author : Zachary Labe
    Date   : 27 September 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyData as MO
import calc_SlicePolarVortexYears as PV
import cmocean
import itertools

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/SortedPVYears/Composites/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Downward Coupling- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = 'GEOP'
period = 'J'
periodall = 'JFMA'

###############################################################################
def readDataPeriods(varnames,simulations,PVstrength,periodpv,periodall):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAveVar(varnames,simulations[0],
                                                'polar','profile')
    lat,lon,lev,varpast = MO.readExperiAveVar(varnames,simulations[1],
                                              'polar','profile')
    
    ### Remove missing data [ensembles,months,lat,lon]
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan

    ### Slice Polar Vortex Years
    if PVstrength == 'all':
        varf = varfuture
        varh = varpast
        print('\n---Number of NEW ensembles is %s!---\n' % (varfuture.shape[0]))
    else:        
        PVq = PV.polarVortexStats(simulations[0],simulations[1],'U10','surface',
                                  PVstrength,period)
        varf = varfuture[PVq,:,:]
        varh= varpast[PVq,:,:]
        print('\n---Number of NEW ensembles is %s!---\n' % (varfuture.shape[0]))
    
    ### Slice by time period
    if periodall == 'JFMA':
        varff = varf[:,0:4,:]
        varhh = varh[:,0:4,:]
        
    ### Calculate anomalies
    anom = np.nanmean(varff - varhh,axis=0)
    
    return anom,lev

### Call data
anomw_weak,lev = readDataPeriods(varnames,['Future','Past'],'<33',
                                       period,periodall)
anomw_strong,lev = readDataPeriods(varnames,['Future','Past'],'>66',
                                         period,periodall)
anomw_all,lev = readDataPeriods(varnames,['Future','Past'],'all',
                                         period,periodall)

anome_weak,leve = readDataPeriods(varnames,['E3SM_Fu','E3SM_Pi'],'<33',
                             period,periodall)
anome_strong,leve = readDataPeriods(varnames,['E3SM_Fu','E3SM_Pi'],'>66',
                             period,periodall)
anome_all,leve = readDataPeriods(varnames,['E3SM_Fu','E3SM_Pi'],'all',
                             period,periodall)

plotall = np.array([anomw_weak,anomw_strong,anome_weak[:,:-1],anome_strong[:,:-1]])

###########################################################################
###########################################################################
###########################################################################
### Read in data   
### Plot Variables
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
if varnames == 'U':
    limit = np.arange(-2,2.1,0.1)
    barlim = np.arange(-2,3,1)
elif varnames == 'TEMP':
    limit = np.arange(-4,4.1,0.2)
    barlim = np.arange(-4,5,1)
elif varnames == 'GEOP':
    limit = np.arange(-200,201,10)
    barlim = np.arange(-200,201,100)
elif varnames == 'V':
    limit = np.arange(-0.2,0.21,0.02)
    barlim = np.arange(-0.2,0.3,0.1)
elif varnames == 'EGR':
    limit = np.arange(-0.08,0.081,0.005)
    barlim = np.arange(-0.08,0.09,0.04)
    
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
time = np.arange(0,4,1)
latq,levq = np.meshgrid(time,lev)

fig = plt.figure()
for i in range(len(plotall)):
    ax1 = plt.subplot(2,2,i+1)

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
            
    cs = plt.contourf(time,lev,plotall[i].transpose(),limit,extend='both') 
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    periodxticks = [r'JAN',r'FEB',r'MAR',r'APR']
    plt.xticks(np.arange(0,4,1),periodxticks,fontsize=5)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=5)
    plt.minorticks_off()
    
    plt.xlim([0,3])
    plt.ylim([1000,10])
    
    if any([i==0]):
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')
        ax1.tick_params(labelbottom=False) 
    elif any([i==1,i==3,i==5,i==7]):
        ax1.tick_params(labelleft=False) 
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                                width=0,color='dimgrey')
    if any([i==1]):
        ax1.tick_params(labelbottom=False) 
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey')  
   
    if varnames == 'U':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames == 'TEMP':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames == 'GEOP':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap) 
    elif varnames == 'V':
        cmap = cmocean.cm.balance             
        cs.set_cmap(cmap) 
    elif varnames == 'EGR':
        cmap = cmocean.cm.diff           
        cs.set_cmap(cmap) 

    labelpv = [r'Weak PV',r'Strong PV']
    if i < 2:
        ax1.annotate(r'\textbf{%s}' % labelpv[i],
                    xy=(0, 0),xytext=(0.5,1.08),xycoords='axes fraction',
                    fontsize=13,color='dimgrey',rotation=0,
                    ha='center',va='center')
    if i==0:     
        plt.annotate(r'\textbf{$\bf{\Delta}$WACCM}',
            xy=(0, 0),xytext=(-0.17,0.5),xycoords='axes fraction',
            fontsize=12,color='k',rotation=90,
            ha='center',va='center')  
    elif i==2:     
        plt.annotate(r'\textbf{$\bf{\Delta}$E3SM}',
            xy=(0, 0),xytext=(-0.17,0.5),xycoords='axes fraction',
            fontsize=12,color='k',rotation=90,
            ha='center',va='center')  

cbar_ax = fig.add_axes([0.312,0.07,0.4,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)

if varnames == 'U':
    cbar.set_label(r'\textbf{m/s}',fontsize=9,color='dimgray',
                   labelpad=0)
elif varnames == 'TEMP':
    cbar.set_label(r'\textbf{$^\circ$C}',fontsize=9,color='dimgray',
                   labelpad=0)
elif varnames == 'GEOP':
    cbar.set_label(r'\textbf{m}',fontsize=9,color='dimgray',
                   labelpad=0)
elif varnames == 'V':
    cbar.set_label(r'\textbf{m/s}',fontsize=9,color='dimgray',
                   labelpad=0)
elif varnames == 'EGR':
    cbar.set_label(r'\textbf{1/day}',fontsize=9,color='dimgray',
                   labelpad=0)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.5)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(hspace=0.05,bottom=0.14,top=0.93,wspace=0.08)

plt.savefig(directoryfigure + '%s_PVCompositesResponse.png' % varnames,
            dpi=300)
print('Completed: Script done!')


    

