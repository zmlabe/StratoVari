"""
Plot vertical plots of PAMIP data segmented over (3) 100 year periods

Notes
-----
    Author : Zachary Labe
    Date   : 21 June 2019
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
print('\n' '----Plotting 100 Year Periods (Vertical)- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U','GEOP','TEMP','V','EGR']
period = 'ND' # Enter temporal period (DJF,JFM,JFMA,ND)
simuh = 'Current' # Enter simulation time (Current,Past)
######################
if simuh == 'Current':
    simuq = 'Cu'
elif simuh == 'Past':
    simuq = 'Pi'
else:
    print(ValueError('Wrong simulation selected!'))
######################
letters = [r'Mean',r'A',r'B',r'C']
for v in range(len(varnames)):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll('%s' % varnames[v],'Future',
                                               'profile')
    lat,lon,lev,varpast = MO.readExperiAll('%s' % varnames[v],'%s' % simuh,
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
    else:
        ValueError('Wrong period selected! (DJF,JFM,JFMA,ND)')
        
    ### Remove missing data 
    varmo[np.where(varmo < -1e10)] = np.nan
        
    ### Take total ensemble mean (300)
    meanens = np.nanmean(varmo,axis=1)
    
    ### Calculate climos
    climo = meanens[1,:,:,:]
    
    ### Calculate anomaly
    meananom = meanens[0,:,:,:] - meanens[1,:,:,:]
    
    ### Calculate zonal means
    meananomz = np.nanmean(meananom,axis=2)
    meanclimoz = np.nanmean(climo,axis=2)
    
    ### Calculate periods (3 - A/B/C)
    climoA = np.nanmean(varmo[1,:100,:,:,:],axis=0)
    climoB = np.nanmean(varmo[1,100:200,:,:,:],axis=0)
    climoC = np.nanmean(varmo[1,200:,:,:,:],axis=0)
    
    meanA = np.nanmean(varmo[0,:100,:,:,:],axis=0) - climoA
    meanB = np.nanmean(varmo[0,100:200,:,:,:],axis=0) - climoB 
    meanC = np.nanmean(varmo[0,200:,:,:,:],axis=0) - climoC
    
    climoAz = np.nanmean(climoA[:,:,:],axis=2)
    climoBz = np.nanmean(climoB[:,:,:],axis=2)
    climoCz = np.nanmean(climoC[:,:,:],axis=2)
    meanAz = np.nanmean(meanA[:,:,:],axis=2)
    meanBz = np.nanmean(meanB[:,:,:],axis=2)
    meanCz = np.nanmean(meanC[:,:,:],axis=2)
    
    ### Calculate significance testing at 95% confidence level
    meanx = np.nanmean(varmo[0,:,:,:,:],axis=3)
    meany = np.nanmean(varmo[1,:,:,:,:],axis=3)
    Ax = np.nanmean(varmo[0,:100,:,:,:],axis=3)
    Ay = np.nanmean(varmo[1,:100,:,:,:],axis=3)
    Bx = np.nanmean(varmo[0,100:200,:,:,:],axis=3)
    By = np.nanmean(varmo[1,100:200,:,:,:],axis=3)
    Cx = np.nanmean(varmo[0,200:,:,:,:],axis=3)
    Cy = np.nanmean(varmo[1,200:,:,:,:],axis=3)
    stat_mean,pvalue_mean = UT.calc_indttest(meanx,meany)
    stat_A,pvalue_A = UT.calc_indttest(Ax,Ay)
    stat_B,pvalue_B = UT.calc_indttest(Bx,By)
    stat_C,pvalue_C = UT.calc_indttest(Cx,Cy)

    ### Append lists of variables
    varall = [meananomz,meanAz,meanBz,meanCz]
    climoall = [meanclimoz,climoAz,climoBz,climoCz]
    pvalall = [pvalue_mean,pvalue_A,pvalue_B,pvalue_C]

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
    
    fig = plt.figure(figsize=(5,7))
    for i in range(len(varall)):
        ax1 = plt.subplot(2,2,i+1)
        
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
                            width=0,color='dimgrey',labelbottom=False)    
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
                            width=0,color='dimgrey',labelbottom=False)    
        elif i == 2:
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
        elif i == 3:
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
        
        
        cs = plt.contourf(lat,lev,var,limit,extend='both')
        
        if varnames[v] == 'U': 
            cs2 = plt.contour(lat,lev,climovar,np.arange(-20,101,5),
                              linewidths=0.6,colors='dimgrey')
        plt.contourf(latq,levq,pvar,colors='None',hatches=['////'],
                     linewidth=5)   
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xlim([0,90])
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
        plt.minorticks_off()
        
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
        
    ### Add experiment text to subplot
        ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
                xytext=(0.03,0.94),xycoords='axes fraction',
                color='k',fontsize=11)
    cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
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
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15,hspace=0.05,wspace=0.05)
    
    plt.savefig(directoryfigure +'%s/' % simuq + 'MeanResponse_300ens_%s_%s.png' % (
                varnames[v],period),dpi=300)
    print('Completed: Script done!')

    