"""
Plot vertical plots of PAMIP data for AMS presentation

Notes
-----
    Author : Zachary Labe
    Date   : 10 December 2019
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
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Vertical Profiles for PAMIP- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['TEMP']
period = 'DJFM' # Enter temporal period (DJF,DJFM,JFM,JFMA,ND,MA)
letters = [r'$\Delta$SIC',r'$\Delta$SIC+$\Delta$SIT']
######################
def readData(simuh,simuf,varnames,period):
    for v in range(len(varnames)):
        ### Call function for 4d variable data
        lat,lon,lev,varfuture = MO.readExperiAll('%s' % varnames[v],simuf,
                                                   'profile')
        lat,lon,lev,varpast = MO.readExperiAll('%s' % varnames[v],simuh,
                                                   'profile')
        
        ### Create 2d array of latitude and longitude
        lon2,lat2 = np.meshgrid(lon,lat)
        
        ### List of experiments
        runs = [varfuture,varpast]
        
        ### Separate per monthly periods
        if period == 'DJF':
            varmo = np.empty((len(runs),varpast.shape[0]-1,varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                      lon,'profile',17)  
        elif period == 'DJFM':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo1 = runs[i][:,:3,:,:,:]
                varmo2 = runs[i][:,-1,:,:,:]
                varmo2 = np.expand_dims(varmo2,axis=1)
                varmoall = np.append(varmo1,varmo2,axis=1)
                varmo[i] = np.nanmean(varmoall,axis=1)
        elif period == 'JFM':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = np.nanmean(runs[i][:,:3,:,:,:],axis=1)
        elif period == 'JFMA':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = np.nanmean(runs[i][:,:4,:,:,:],axis=1)
        elif period == 'ND':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = np.nanmean(runs[i][:,-2:,:,:,:],axis=1)
        elif period == 'MA':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = np.nanmean(runs[i][:,2:4,:,:,:],axis=1)
        else:
            ValueError('Wrong period selected! (DJF,JFM,JFMA,ND)')
            
        ### Remove missing data 
        varmo[np.where(varmo < -1e10)] = np.nan
            
        ### Take total ensemble mean (300)
        meanens = np.nanmean(varmo,axis=1)
        
        ### Calculate climos
        climo = np.nanmean(meanens[1,:,:,:],axis=2)
        
        ### Calculate anomaly
        anom = np.nanmean(meanens[0,:,:,:] - meanens[1,:,:,:],axis=2)
        
        ### Calculate significance testing at 95% confidence level
        meanx = np.nanmean(varmo[0,:,:,:,:],axis=3)
        meany = np.nanmean(varmo[1,:,:,:,:],axis=3)      
        pvalue = UT.calc_FDR_ttest(meanx,meany,0.05) #FDR
        
        pvalue[np.where(pvalue < 0.05)] = 1.
        pvalue[np.where(np.isnan(pvalue))] = 0.
        
        return anom,climo,pvalue,lev,lat
    
anom_sitq,climo_sitq,pvalue_sitq,lev,lat = readData('SIT_Cu','SIT_Fu',varnames,period)
anom_sic,climo_sic,pvalue_sic,lev,lat = readData('Past','Future',varnames,period)

### Adjust levels
anom_sit = anom_sitq[:-1,:]
climo_sit = climo_sitq[:-1,:]
pvalue_sit = pvalue_sitq[:-1,:]

varall = [anom_sic,anom_sit]
climoall = [climo_sic,climo_sit]
pvalall = [pvalue_sic,pvalue_sit]

###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black') 

### Set limits for contours and colorbars
if varnames[0] == 'U':
    limit = np.arange(-1,1.01,0.05)
    barlim = np.arange(-1,2,1)
elif varnames[0] == 'TEMP':
    limit = np.arange(-3,3.1,0.25)
    barlim = np.arange(-3,4,1)
elif varnames[0] == 'GEOP':
    limit = np.arange(-60,61,2)
    barlim = np.arange(-60,61,30)
elif varnames[0] == 'V':
    limit = np.arange(-0.2,0.21,0.02)
    barlim = np.arange(-0.2,0.3,0.1)
elif varnames[0] == 'EGR':
    limit = np.arange(-0.08,0.081,0.005)
    barlim = np.arange(-0.08,0.09,0.04)
    
zscale = np.array([1000,700,500,300])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(len(varall)):
    ax1 = plt.subplot(1,2,i+1)
    
    var = varall[i]
    pvar = pvalall[i]
    climovar = climoall[i]
    
    if i == 0:
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
                        width=2,color='dimgrey',labelbottom=True)    
        ax1.yaxis.set_ticks_position('left')
    elif i == 1:
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
                        width=2,color='dimgrey',labelbottom=True)    
       
    cs = plt.contourf(lat,lev,var*pvar,limit,extend='both')

    if varnames[0] == 'U': 
        cs2 = plt.contour(lat,lev,climovar,np.arange(-20,101,5),
                          linewidths=0.6,colors='dimgrey')  
    if varnames[0] == 'TEMP': 
        cs2 = plt.contour(lat,lev,climovar,np.arange(-70,101,5),
                          linewidths=0.3,colors='dimgrey')  
    
    plt.yscale('log',nonposy='clip')
    plt.xlim([0,90])
    plt.ylim([1000,300])
    plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
    plt.minorticks_off()
    plt.gca().invert_yaxis()
    
    if varnames[0] == 'U':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames[0] == 'TEMP':
        cmap = cmocean.cm.balance            
        cs.set_cmap(cmap) 
    elif varnames[0] == 'GEOP':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap) 
    elif varnames[0] == 'V':
        cmap = cmocean.cm.balance             
        cs.set_cmap(cmap) 
    elif varnames[0] == 'EGR':
        cmap = cmocean.cm.diff           
        cs.set_cmap(cmap) 
        
    ax1.set_aspect('equal')
    plt.title(r'\textbf{%s}' % letters[i],color='w',fontsize=30)
    
### Add experiment text to subplot
cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

if varnames[0] == 'U':
    cbar.set_label(r'\textbf{[U] m/s}',fontsize=8,color='darkgray',labelpad=0)
elif varnames[0] == 'TEMP':
    cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='darkgray',labelpad=0)
elif varnames[0] == 'GEOP':
    cbar.set_label(r'\textbf{m}',fontsize=11,color='darkgray',labelpad=0)
elif varnames[0] == 'V':
    cbar.set_label(r'\textbf{m/s}',fontsize=11,color='darkgray',labelpad=0)
elif varnames[0] == 'EGR':
    cbar.set_label(r'\textbf{1/day}',fontsize=11,color='darkgray',labelpad=0)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelsize=7)
cbar.outline.set_edgecolor('darkgrey')

plt.tight_layout()
plt.subplots_adjust(bottom=0.18,hspace=0.05,wspace=0.09)

plt.savefig(directoryfigure + 'AMS_MeanResponse_300ens_%s_%s_FDR2.png' % (
            varnames[0],period),dpi=300)
print('Completed: Script done!')

    