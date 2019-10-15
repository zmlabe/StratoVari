"""
Composite SSW events for E3SM 1.1

Notes
-----
    Author : Zachary Labe
    Date   : 15 October 2019
"""

### Import modules
import numpy as np
import cmocean
import matplotlib.pyplot as plt
import scipy.stats as sts

### Directories
directoryfigure = '/home/zlabe/Desktop/'

def readData(exp):
    """
    Read in data for E3SM experiments
    """
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    if exp == '1.1':
        experiment = 'PAMIP-1.1-E3SM'
        NENS = 100
        
    directorydata = '/seley/zlabe/simu/' + experiment + '/daily/'
    
    ### Read in geopotential polar cap data
    polarcap = np.empty((NENS-1,365*2,17))
    for i in range(1,NENS,1):
        filepath = directorydata + experiment + '_%s/' % i
        filename = filepath + 'GEOP_%s_Extended_polarmean.nc' % i
        
        data = Dataset(filename,'r')
        lev = data.variables['level'][:-1]
        polarcap[i-1,:,:] = data.variables['GEOP'][:,:-1].squeeze()
        data.close()
        
    print('Completed: Data read for polar cap!')
        
    ### Read in SSW Count
    count = np.empty((NENS-1))
    for i in range(1,NENS,1):
        filepath = directorydata + experiment + '_%s/' % i
        filename = filepath + 'sswcount_%s_11-3_winter.txt' % i
        count[i-1] = np.genfromtxt(filename,unpack=True,delimiter=',')
        
    print('Completed: Data read for ssw count!')
        
    ### Read in SSW Stats
    stat = np.empty((NENS-1,polarcap.shape[1]))
    for i in range(1,NENS,1):
        filepath = directorydata + experiment + '_%s/' % i
        filename = filepath + 'sswstats_%s_11-3_winter.txt' % i
        statq = np.genfromtxt(filename,usecols=[6])
        
        ### Add data for January 1 - October 31
        zeroJO = np.zeros((304))
        statqJO = np.append(zeroJO,statq)
        ### Add data for April 1 - December 31
        zeroAD = np.zeros((275))
        statqall = np.append(statqJO,zeroAD)
        
        stat[i-1,:] = statqall
        
    print('Completed: Data read for ssw statistics!')
    
    return polarcap,count,stat,lev

def calculateAnom(var):
    """
    Calculate polar cap anomalies
    """
    ### Import modules
    import numpy as np
    
    mean = np.nanmean(var,axis=0)
    anom = var - mean
    
    print('Completed: Calculate anomalies for polar cap!')
    return anom
    
### Call functions for data
polarcap,count,stat,lev = readData('1.1')
anomgeo = calculateAnom(polarcap)

### Locate SSW 
NENS = polarcap.shape[0]
days = polarcap.shape[1]
levq = lev.shape[0]

indexssw = []
for i in range(NENS):
    indexq = np.where(stat[i,:])[0]
    indexssw.append(indexq)
    
### Average 30 days prior and 90 days after
prior = 30
after = 90
time = np.arange(prior+after)
sswperiods = np.zeros((NENS,prior+after,levq))
for i in range(i):
    if int(count[i]) == 1:
        indexval = indexssw[i][0]
        slicegeomean = anomgeo[i,indexval-prior:indexval+after,:]
    elif int(count[i]) == 2:
        indexval1 = indexssw[i][0]
        slicegeo1 = anomgeo[i,indexval1-prior:indexval1+after,:]
        indexval2 = indexssw[i][1]
        slicegeo2 = anomgeo[i,indexval2-prior:indexval2+after,:]
        slicegeomean = (slicegeo1 + slicegeo2)/2
    elif int(count[i]) == 3:
        indexval1 = indexssw[i][0]
        slicegeo1 = anomgeo[i,indexval1-prior:indexval1+after,:]
        indexval2 = indexssw[i][1]
        slicegeo2 = anomgeo[i,indexval2-prior:indexval2+after,:]
        indexval3 = indexssw[i][1]
        slicegeo3 = anomgeo[i,indexval3-prior:indexval3+after,:]
        slicegeomean = (slicegeo1 + slicegeo2 + slicegeo3)/3
    else:
        empty = np.zeros((sswperiods.shape[1],sswperiods.shape[2]))
        empty[:] = np.nan
        slicegeomean = empty
        
    sswperiods[i,:,:] = slicegeomean
    
### Calculate ensemble mean climatology for SSW events
meanssw = np.nanmean(sswperiods,axis=0)

### Plot Variables
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-300,301,5)
barlim = np.arange(-300,301,100)
    
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
timeq,levq = np.meshgrid(time,lev)

fig = plt.figure()
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
        
cs = plt.contourf(time,lev,meanssw.transpose(),limit,extend='both')
cs1 = plt.contour(time,lev,meanssw.transpose(),np.arange(-900,901,50),
                  extend='both',linewidths=0.4,colors='dimgrey')

plt.axvline(x=31,linestyle='--',dashes=(1,0.3),color='k',linewidth=2)

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.xticks(np.arange(0,121,30),map(str,np.arange(0,121,30)),fontsize=5)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=5)
plt.minorticks_off()

plt.xlabel(r'\textbf{Days}',fontsize=8)
plt.ylabel(r'\textbf{Pressure (hPa)}',fontsize=8)
plt.title(r'\textbf{Climatology of SSW Events -- Polar Cap Heights -- E3SM-1.1}',fontsize=11)

plt.xlim([0,120])
plt.ylim([1000,10])

ax1.tick_params(axis='y',direction='out',which='major',pad=1,
                width=2,color='dimgrey')
ax1.tick_params(axis='x',direction='out',which='major',pad=1,
                width=2,color='dimgrey')  

cmap = cmocean.cm.balance             
cs.set_cmap(cmap) 

cbar_ax = fig.add_axes([0.312,0.06,0.4,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{m}',fontsize=9,color='k',
               labelpad=0)

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.5)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(bottom=0.17,top=0.93)

plt.savefig(directoryfigure + 'SSWEvents_E3SM1.1.png',dpi=300)
print('Completed: Script done!')

np.savetxt('/home/zlabe/Documents/Research/StratoVari/Data/E3SM_SSWClimo.txt',
           meanssw)
np.savetxt('/home/zlabe/Documents/Research/StratoVari/Data/Levels.txt',
           lev)