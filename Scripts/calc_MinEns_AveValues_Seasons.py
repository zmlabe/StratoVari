"""
Script calculates the 'average' minimum number of ensembles needed for
each variable for seasonal means

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
    
    runs = [varfuture,varpast]
    
### Separate per monthly periods
    period = 'DJF'
    if period == 'DJF':
        varmo = np.empty((len(runs),varpast.shape[0]-1,varpast.shape[2],
                          varpast.shape[3]))
        for i in range(len(runs)):
            varmo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'surface',17) 
        varfuturem = varmo[0]
        varpastm = varmo[1]
    elif period == 'NDJFMA':
        varfuturem = np.nanmean(np.append(varfuture[:,-2:,:,:],
                                          varfuture[:,:4,:,:],
                                          axis=1),axis=1)
        varpastm = np.nanmean(np.append(varpast[:,-2:,:,:],
                                        varpast[:,:4,:,:],
                                        axis=1),axis=1)
    else:
        ValueError('Wrong period selected! (DJF,JFM,JFMA,ND)')

    return varfuturem,varpastm,lat,lon
    
###########################################################################
###########################################################################
###########################################################################
### Read in data
sliceq = 'A'
simu = 'Past'
weighted = True
allnmin = np.empty((len(varnames)))
for v in range(len(varnames)):
    future,climo,lat,lon = readDataPeriods(varnames[v],sliceq,simu)
    nmin = MENS.computeMinEns(future,climo,0.05)
    
    ### Calculate only for north of 30N
    latq = np.where(lat > 30)[0]
    lat = lat[latq]
    nmin = nmin[latq,:]
        
    ### Meshgrid of lat and lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Mask values that are not significant
    nminmask = nmin.copy()
    nminmask[np.isnan(nminmask)] = 0.
    nminmask[np.where(nminmask > 0.)] = 1.
    
    ### Calculated weighted average for only significant points
    if weighted == True:
        avenmin = UT.calc_weightedAve(nmin,lat2)
    elif weighted == False:
        avenmin = np.nanmean(nmin[:,:])
    
    ### Save to look at all values
    allnmin[v] = avenmin
    
#    ### Save individual files
#    np.savetxt(directorydataout + '%s_minens95_DJF_%s.txt' % (simu,
#               varnames[v]),np.array([avenmin]),delimiter=',',fmt='%3.1f',
#               footer='\n Minimum number of ensembles needed to get' \
#               '\n statistical significance at the 95% confidence\n' \
#               ' level for DJF',newline='\n\n')
        