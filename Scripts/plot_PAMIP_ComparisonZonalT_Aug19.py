"""
This script compares all PAMIP simulations as of 8/28/2019 for the zonal mean
temperature response to historical and preindustrial forcings. All statistical
test uses the FDR method with alpha_FDR=0.05

Notes
-----
    Author : Zachary Labe
    Date   : 28 August 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyData as MO
import read_MonthlyData_ANT as MOA
import read_MonthlyData_SST as MOS
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
print('\n' '----Plotting Zonal Mean T - %s----' % titletime)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
letters = list(string.ascii_lowercase)
periodm = 'JFM'
periodma = 'JAS'
dataread = True
vari = 'TEMP'

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Aug19_PAMIPcomp/'

######################
def readDataPeriods(typesimu,varnames,simulations,period):
    ### Call function for 4d variable data
    if typesimu == 'arc_sea_ice':
        lat,lon,lev,varfuture = MO.readExperiAll(varnames,
                                                 simulations[0],'zonmean')
        lat,lon,lev,varpast = MO.readExperiAll(varnames,
                                               simulations[1],'zonmean')
    elif typesimu == 'sst':
        lat,lon,lev,varfuture = MOS.readExperiAll(varnames,
                                                 simulations[0],'zonmean')
        lat,lon,lev,varpast = MOS.readExperiAll(varnames,
                                               simulations[1],'zonmean')
    elif typesimu == 'ant_sea_ice':
        lat,lon,lev,varfuture = MOA.readExperiAll(varnames,
                                                 simulations[0],'zonmean')
        lat,lon,lev,varpast = MOA.readExperiAll(varnames,
                                               simulations[1],'zonmean')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data [ensembles,months,lat,lon]
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Calculate over DJF
    if period == 'JJA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,5:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,5:8:,:,:],axis=1)
    elif period == 'JAS':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:9,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:9:,:,:],axis=1)
    elif period == 'July':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:7,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:7:,:,:],axis=1)
    elif period == 'August':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,7:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,7:8:,:,:],axis=1)
    elif period == 'JA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:8:,:,:],axis=1)
    elif period == 'S':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,8:9,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,8:9:,:,:],axis=1)
    elif period == 'AMJ':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,3:6,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,3:6:,:,:],axis=1)
    elif period == 'OND':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,-3:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,-3:,:,:],axis=1)
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
    elif period == 'Annual':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,:,:,:],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))
    
    ### Calculate anomalies and ensemble means
    varfutureq = varfuturem.copy()
    varpastq = varpastm.copy()
    anommean = np.nanmean(varfutureq,axis=0) - np.nanmean(varpastq,axis=0)
    
    ### Calculate climatologies
    climo = np.nanmean(varpastq,axis=0)
    
    ### Calculate significance for each month (pick method)
    pruns = UT.calc_FDR_ttest(varfuturem[:,:,:],varpastm[:,:,:],0.05) #FDR
    
    return anommean,climo,pruns,lat,lev

###############################################################################
###############################################################################
###############################################################################
### Read in data
if dataread == True:
    ###############################################################################
    anom_sstcu,climo_sstcu,p_sstcu,lat,lev = readDataPeriods('sst',vari,
                                                             ['SST_Fu','SST_Cu'],
                                                             periodm)
    anom_sstpi,climo_sstpi,p_sstpi,lat,lev = readDataPeriods('sst',vari,
                                                             ['SST_Fu','SST_Pi'],
                                                             periodm)
    ###############################################################################
    anom_antcu,climo_antcu,p_antcu,lat,lev = readDataPeriods('ant_sea_ice',vari,
                                                             ['ANT_Fu','ANT_Cu'],
                                                             periodma)
    anom_antpi,climo_antpi,p_antpi,lat,lev = readDataPeriods('ant_sea_ice',vari,
                                                             ['ANT_Fu','ANT_Pi'],
                                                             periodma)
    ###############################################################################
    anom_BKcu,climo_BKcu,p_BKcu,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['BKsea_Fu','Current'],
                                                             'JFM')
    anom_BKpi,climo_BKpi,p_BKpi,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['BKsea_Fu','Past'],
                                                             periodm)
    ###############################################################################
    anom_Ocu,climo_Ocu,p_Ocu,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['Osea_Fu','Current'],
                                                             periodm)
    anom_Opi,climo_Opi,p_Opi,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['Osea_Fu','Past'],
                                                             periodm)
    ###############################################################################
    anom_ARCcu,climo_ARCcu,p_ARCcu,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['Future','Current'],
                                                             periodm)
    anom_ARCpi,climo_ARCpi,p_ARCpi,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['Future','Past'],
                                                             periodm)
    ###############################################################################
    anom_SITcu,climo_SITcu,p_SITcu,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['SIT_Fu','SIT_Cu'],
                                                             periodm)
    anom_SITpi,climo_SITpi,p_SITpi,lat,lev = readDataPeriods('arc_sea_ice',vari,
                                                             ['SIT_Fu','Past'],
                                                             periodm)
    ###############################################################################
    #anom_E3cu,climo_E3cu,p_E3cu,lat,lev = readDataPeriods('arc_sea_ice',vari,
    #                                                         ['E3SM_Fu','E3SM_Cu'],
    #                                                         periodm)
    #anom_E3pi,climo_E3pi,p_E3pi,lat,lev = readDataPeriods('arc_sea_ice',vari,
    #                                                         ['E3SM_Fu','E3SM_Pi'],
    #                                                         periodm)
###############################################################################
### Assign into lists for plotting
datas = np.array([anom_sstpi,anom_sstcu,
                 anom_ARCpi,anom_ARCcu,
                 anom_SITpi,anom_SITcu,
                 anom_BKpi,anom_BKcu,
                 anom_Opi,anom_Ocu])
climos = np.array([climo_sstpi,climo_sstcu,
                 climo_ARCpi,climo_ARCcu,
                 climo_SITpi,climo_SITcu,
                 climo_BKpi,climo_BKcu,
                 climo_Opi,climo_Ocu])
pvals = np.array([p_sstpi,p_sstcu,
                 p_ARCpi,p_ARCcu,
                 p_SITpi,p_SITcu,
                 p_BKpi,p_BKcu,
                 p_Opi,p_Ocu])

###########################################################################
###########################################################################
###########################################################################
#### Plot T
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-4,4.1,0.25)
barlim = np.arange(-4,5,2)
zscale = np.array([1000,700,500,300])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(len(datas)):
    anom = datas[i]
    climo = climos[i]
    pval = pvals[i]
    
    ax1 = plt.subplot(5,2,i+1)
    
    if i == 8:
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
    elif i == 9:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
    elif any([i==0,i==2,i==4,i==6]):
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey',labelleft=True)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='w',labelbottom=False)    
    elif any([i==1,i==3,i==5,i==7]):
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='w',labelbottom=False)    
    
    ### Set limits    
    cs = plt.contourf(lat,lev,anom,limit,extend='both')
    cs1 = plt.contour(lat,lev,climo,np.arange(-100,110,10),
                      linewidths=0.5,colors='dimgrey') 
    cs2 = plt.contourf(lat,lev,pval,
                       colors='None',hatches=['///////'])     
    
    cs.set_cmap(cmocean.cm.balance)
    
    plt.yscale('log',nonposy='clip')
    plt.xlim([-90,90])
    plt.ylim([1000,300])
    plt.xticks(np.arange(-90,96,15),map(str,np.arange(-90,91,15)),fontsize=4)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=4)
    plt.minorticks_off()
    plt.gca().invert_yaxis()
    
    if i==0:
        plt.title(r'\textbf{$\bf{\Delta}$Pi}',fontsize=11,color='k')    
        plt.ylabel(r'\textbf{SST}',fontsize=6,color='k')   
    if i==1:
        plt.title(r'\textbf{$\bf{\Delta}$Cu}',fontsize=11,color='k')    
    if i==2:
        plt.ylabel(r'\textbf{Arctic}',fontsize=6,color='k')   
    if i==4:
        plt.ylabel(r'\textbf{Thickness}',fontsize=6,color='k')   
    if i==6:
        plt.ylabel(r'\textbf{BK-Sea}',fontsize=6,color='k')   
    if i==8:
        plt.ylabel(r'\textbf{O-Sea}',fontsize=6,color='k')    
        
        
#    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
#            xytext=(0.01,0.84),xycoords='axes fraction',
#            color='k',fontsize=11)

cbar_ax = fig.add_axes([0.326,0.08,0.4,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{T [$\bf{^{\circ}}$C]}',fontsize=8,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15,hspace=0.15,wspace=0.05,top=0.94)

plt.savefig(directoryfigure + 'T_zonalmean_PAMIPcomp_Aug19.png',dpi=300)