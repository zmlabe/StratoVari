"""
Script calculates the minimum number of enembles to detect a statistically
significant ensemble mean difference.

Notes
-----
    Author : Zachary Labe
    Date   : 24 June 2019
"""

def readinData(varnames,simuh,period):
    ### Import modules
    import numpy as np
    import datetime
    import read_MonthlyData as MO
    import calc_Utilities as UT
    
    ### Define time           
    now = datetime.datetime.now()
    currentmn = str(now.month)
    currentdy = str(now.day)
    currentyr = str(now.year)
    titletime = currentmn + '/' + currentdy + '/' + currentyr
    print('\n' '----Calculate Minimum Ensembles- %s----' % titletime)
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### Read in data
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
        elif period == 'MA':
            varmo = np.empty((len(runs),varpast.shape[0],varpast.shape[2],
                              varpast.shape[3],varpast.shape[4]))
            for i in range(len(runs)):
                varmo[i] = np.nanmean(runs[i][:,2:4,:,:,:],axis=1)
        else:
            ValueError('Wrong period selected! (DJF,JFM,JFMA,ND)')
            
        ### Remove missing data 
        varmo[np.where(varmo < -1e10)] = np.nan
        
        ### Simulation
        control = varmo[1,:,:,:,:]
        future = varmo[0,:,:,:,:]
        
        ### Calculate anomaly
        anom = future - control
        
        
    return anom,future,control,lat,lon,lev

def computeMean(datain,typemean):
    """
    Compute simple means for different dimensions of array
    """
    ### Import modules
    import numpy as np
    
    ### Calculate various means
    if typemean == 'zonal': # zonal mean
        dataout = np.nanmean(datain,axis=datain.ndim-1)
    elif typemean == 'ensemble': # ensemble mean
        dataout = np.nanmean(datain,axis=0)
     
    return dataout

def computeSTD(datain,df):
    """
    Compute standard deviation of ensemble members
    """
    ### Import modules
    import numpy as np
    
    ### Compute standard deviation
    dataout = np.nanstd(datain,axis=0,ddof=df,
                        dtype=np.float64) # calculate for ensemble members
    
    return dataout

def computePooledSD(xstd,ystd,xn,yn):
    """
    Compute pooled standard deviation
    """
    ### Import modules
    import numpy as np
    
    ### Compute pooled standard deviation
    if xstd.ndim == 2:
        sp = np.empty((xstd.shape))
        for i in range(xstd.shape[0]):
            for j in range(xstd.shape[1]):
                sp[i,j] = np.sqrt(((xn - 1)*xstd[i,j]**2 + \
                          (yn - 1)*ystd[i,j]**2) \
                          /(xn + yn - 2))
    
    return sp

def computeMinEns(sp,xmean,ymean,alpha):
    """
    Compute minimum ensemble number using formula 4 from 
    Screen et al. 2013, Climate Dynamics
    """
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### Calculate t statistic for confidence level
    tc = sts.t.ppf(1-(alpha/2),len(futurem)-1) # two-tailed
    
    ### Compute minimum ensemble number
    nmin = (2*tc**2) * (sp/(xmean - ymean))**2
    
    nmin[np.where(nmin >= 300)] = np.nan
    
    return nmin
    
###############################################################################
###############################################################################
###############################################################################
### Calculate functions
simuh = 'Current' # Enter simulation time (Current,Past)
varnames = ['TEMP']
period = 'MA'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/'
######################
if simuh == 'Current':
    simuq = 'Cu'
elif simuh == 'Past':
    simuq = 'Pi'
else:
    print(ValueError('Wrong simulation selected!'))
anom,future,climo,lat,lon,lev = readinData(varnames,simuh,period)

futurez = computeMean(future,'zonal')
futurem = computeMean(futurez,'ensemble')
futurestd = computeSTD(futurez,1)

climoz = computeMean(climo,'zonal')
climom = computeMean(climoz,'ensemble')
climostd = computeSTD(climoz,1)

sp = computePooledSD(futurestd,climostd,len(futurez),len(climoz))
nmin = computeMinEns(sp,futurem,climom,0.05)

###############################################################################
###############################################################################
###############################################################################
#### Plot minimum number of ensembles
### Import modules
import numpy as np
import matplotlib.pyplot as plt
import palettable.cubehelix as cm
import cmocean

### Set parameters for matplotlib
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(0,300.1,20)
barlim = np.arange(0,301,50)
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

### Begin plot
fig = plt.figure(figsize=(5,7))
ax1 = plt.subplot(111)

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

cs = plt.contourf(lat,lev,nmin,limit) 

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')
plt.xlim([0,90])
plt.ylim([1000,10])
plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()

cmap = cm.cubehelix1_16.mpl_colormap     
cs.set_cmap(cmap) 
cbar_ax = fig.add_axes([0.312,0.07,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{Minimum Ensemble Size}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure +'%s/' % simuq + 'MinEns_%s_%s.png' % (
            varnames[0],period),dpi=300)
print('Completed: Script done!')