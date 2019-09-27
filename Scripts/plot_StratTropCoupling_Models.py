"""
Plot vertical plots of PAMIP data for each month from January to April using
the ensemble mean (300) by compositing on PV years

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
import calc_PVComposites as PV
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
def readDataPeriods(varnames,simu,PVstrength,periodpv,periodall):
    ### Call function for 4d variable data
    lat,lon,lev,varq = MO.readExperiAveVar(varnames,simu,'polar','profile')

    ### Slice by years    
    if PVstrength != 'all':
        argyears = PV.polarVortexComp(simu,'U10','surface',PVstrength,periodpv)
        varpv = varq[argyears,:,:]
    elif PVstrength == 'all':
        varpv = varq
    
    ### Slice by time period
    if periodall == 'JFMA':
        var = varpv[:,0:4,:]
    
    return var,lev
###############################################################################
###############################################################################
###############################################################################
### Read in WACCM data
#var_f_weak,lev = readDataPeriods(varnames,'Future','<33',period,periodall)
#var_f_strong,lev = readDataPeriods(varnames,'Future','>66',period,periodall)
#var_f_all,lev = readDataPeriods(varnames,'Future','all',period,periodall)
#var_h_weak,lev = readDataPeriods(varnames,'Past','<33',period,periodall)
#var_h_strong,lev = readDataPeriods(varnames,'Past','>66',period,periodall)
#var_h_all,lev = readDataPeriods(varnames,'Past','all',period,periodall)
#
#anom_f_weak = np.nanmean(var_f_weak,axis=0) - np.nanmean(var_f_all,axis=0)
#anom_h_weak = np.nanmean(var_h_weak,axis=0) - np.nanmean(var_h_all,axis=0)
#anom_f_strong = np.nanmean(var_f_strong,axis=0) - np.nanmean(var_f_all,axis=0)
#anom_h_strong = np.nanmean(var_h_strong,axis=0) - np.nanmean(var_h_all,axis=0)
#
#### Read in E3SM data
#var_ef_weak,leve = readDataPeriods(varnames,'E3SM_Fu','<33',period,periodall)
#var_ef_strong,leve = readDataPeriods(varnames,'E3SM_Fu','>66',period,periodall)
#var_ef_all,leve = readDataPeriods(varnames,'E3SM_Fu','all',period,periodall)
#var_eh_weak,leve = readDataPeriods(varnames,'E3SM_Pi','<33',period,periodall)
#var_eh_strong,leve = readDataPeriods(varnames,'E3SM_Pi','>66',period,periodall)
#var_eh_all,leve = readDataPeriods(varnames,'E3SM_Pi','all',period,periodall)
#
#anom_ef_weak = np.nanmean(var_ef_weak,axis=0) - np.nanmean(var_ef_all,axis=0)
#anom_eh_weak = np.nanmean(var_eh_weak,axis=0) - np.nanmean(var_eh_all,axis=0)
#anom_ef_weak = anom_ef_weak[:,:-1]
#anom_eh_weak = anom_eh_weak[:,:-1]
#anom_ef_strong = np.nanmean(var_ef_strong,axis=0) - np.nanmean(var_ef_all,axis=0)
#anom_eh_strong = np.nanmean(var_eh_strong,axis=0) - np.nanmean(var_eh_all,axis=0)
#anom_ef_strong = anom_ef_strong[:,:-1]
#anom_eh_strong = anom_eh_strong[:,:-1]
#
#### Create arguements for plots
#plotall = np.array([anom_f_weak,anom_f_strong,anom_h_weak,anom_h_strong,
#           anom_ef_weak,anom_ef_strong,anom_eh_weak,anom_eh_strong])
    
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
    ax1 = plt.subplot(4,2,i+1)

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
    
    if any([i==0,i==2,i==4]):
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                            width=0,color='dimgrey')
        ax1.tick_params(labelbottom=False) 
    elif any([i==1,i==3,i==5,i==7]):
        ax1.tick_params(labelleft=False) 
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                                width=0,color='dimgrey')
    if any([i==1,i==3,i==5]):
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
                    xy=(0, 0),xytext=(0.5,1.13),xycoords='axes fraction',
                    fontsize=13,color='dimgrey',rotation=0,
                    ha='center',va='center')
    if i==0:     
        plt.annotate(r'\textbf{WACCM-F}',
            xy=(0, 0),xytext=(-0.2,0.5),xycoords='axes fraction',
            fontsize=8,color='k',rotation=90,
            ha='center',va='center')  
    elif i==2:     
        plt.annotate(r'\textbf{WACCM-H}',
            xy=(0, 0),xytext=(-0.2,0.5),xycoords='axes fraction',
            fontsize=8,color='k',rotation=90,
            ha='center',va='center')  
    elif i==4:     
        plt.annotate(r'\textbf{E3SM-F}',
            xy=(0, 0),xytext=(-0.2,0.5),xycoords='axes fraction',
            fontsize=8,color='k',rotation=90,
            ha='center',va='center')  
    elif i==6:     
        plt.annotate(r'\textbf{E3SM-H}',
            xy=(0, 0),xytext=(-0.2,0.5),xycoords='axes fraction',
            fontsize=8,color='k',rotation=90,
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

plt.subplots_adjust(hspace=0.1,bottom=0.14,top=0.93,wspace=0.08)

plt.savefig(directoryfigure + '%s_PVComposites.png' % varnames,
            dpi=300)
print('Completed: Script done!')


    

