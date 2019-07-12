"""
Script calculates the 'average' minimum number of ensembles needed for
each variable

Notes
-----
    Author : Zachary Labe
    Date   : 12 July 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_MinEns_Maps as MENS
import palettable.cubehelix as cm
import calc_Utilities as UT

### Define directories
directorydata = '/seley/zlabe/simu/'
directorydataout = '/home/zlabe/Documents/Research/StratoVari/Data/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating minimum ensembles!- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','U30','U50','U300','U700','SLP','Z500','Z200','Z30','T2M',
            'THICK','P']
varnamesq = [r'\textbf{U10}',r'\textbf{U30}',r'\textbf{U50}',r'\textbf{U300}',
             r'\textbf{U700}',r'\textbf{SLP}',r'\textbf{Z500}',r'\textbf{Z200}',
             r'\textbf{Z30}',r'\textbf{T2M}',r'\textbf{THICK}',r'\textbf{P}']
months = [r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
monthsq = [r'\textbf{Nov}',r'\textbf{Dec}',r'\textbf{Jan}',
           r'\textbf{Feb}',r'\textbf{Mar}',r'\textbf{Apr}'] 

######################
def readDataPeriods(varnames,sliceq,simu):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,'Future','surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,simu,'surface')
    
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

    return varfuturem,varpastm,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
sliceq = 'Mean'
simu = 'Current'
weighted = False
allnmin = np.empty((len(varnames),len(months)))
for v in range(len(varnames)):
    future,climo,lat,lon = readDataPeriods(varnames[v],sliceq,simu)
    nmin = MENS.computeMinEns(future,climo,0.05)
    
    ### Calculate only for north of 30N
    latq = np.where(lat > 30)[0]
    lat = lat[latq]
    nmin = nmin[:,latq,:]
        
    ### Meshgrid of lat and lon
    lon2,lat2 = np.meshgrid(lon,lat)
    lat2mo = np.repeat(lat2[np.newaxis,:,:],6,axis=0)
    lon2mo = np.repeat(lon2[np.newaxis,:,:],6,axis=0)
    
    ### Mask values that are not significant
    nminmask = nmin.copy()
    nminmask[np.isnan(nminmask)] = 0.
    nminmask[np.where(nminmask > 0.)] = 1.
    
    ### Mask out lat and lon arrays
    ### These lines are not needed for function
    ###lat2mox = lat2mo * nminmask
    ###lon2mox = lon2mo * nminmask
    ###lat2mox[np.where(lat2mox == 0.)] = np.nan
    ###lon2mox[np.where(lon2mo == 0.0)] = np.nan
    
    ### Calculated weighted average for only significant points
    if weighted == True:
        avenmin = UT.calc_weightedAve(nmin,lat2)
    elif weighted == False:
        avenmin = np.empty((len(months)))
        for mo in range(len(months)):
            avenmin[mo] = np.nanmean(nmin[mo,:,:])
    
    ### Save to look at all values
    allnmin[v,:] = avenmin
    
    ### Save individual files
    np.savetxt(directorydataout + '%s_minens95_monthly_%s.txt' % (simu,
               varnames[v]),avenmin,delimiter=',',
               fmt='%3.1f',header='  '.join(months)+'\n',
               footer='\n Minimum number of ensembles needed to get' \
               '\n statistical significance at the 95% confidence\n' \
               ' level for November through April',newline='\n\n')
        
###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.get_xaxis().set_tick_params(direction='out', width=0,length=0,
            color='w')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on')

boundaries = np.arange(0,301,25)
cmap = cm.cubehelix1_16.mpl_colormap

cs = plt.pcolormesh(allnmin,shading='faceted',edgecolor='w',
                    linewidth=0.3,norm=mpl.colors.BoundaryNorm(
                            boundaries,ncolors=cmap.N,clip=False))
cs.set_cmap(cmap)

for i in range(allnmin.shape[0]):
    for j in range(allnmin.shape[1]):
        value = allnmin[i,j]
        if np.isnan(value) == True:
            value = '$>$300'
        else:
            value = int(value)
        plt.text(j+0.5,i+0.5,r'\textbf{%s}' % value,fontsize=6,
                 color='k',va='center',ha='center')

ylabels = varnamesq
plt.yticks(np.arange(0.5,13.5,1),ylabels,ha='right',color='dimgrey',
           va='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=0.7)

xlabels = monthsq
plt.xticks(np.arange(0.5,7.5,1),xlabels,ha='center',color='dimgrey',
           va='center')
xax = ax.get_xaxis()
xax.set_tick_params(pad=8)
plt.xlim([0,6])

cbar = plt.colorbar(cs,orientation='horizontal',aspect=50,drawedges=False,
                    extend='max')
ticks = np.arange(0,301,50)
cbar.set_ticks(ticks)
cbar.set_ticklabels(list(map(str,ticks))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(1.2)
cbar.ax.tick_params(labelsize=10)
cbar.set_label(r'\textbf{Minimum Ensemble Members}',
               fontsize=13,color='dimgray',labelpad=3)  

plt.tight_layout()

if weighted == True:
    plt.savefig(directoryfigure + '%s_minens.png' % simu,dpi=300)
elif weighted == False:
    plt.savefig(directoryfigure + '%s_minens.png' % simu,dpi=300)

        