"""
Compare SSW climatologies for daily geopotential polar cap height data

Notes
-----
    Author : Zachary Labe
    Date   : 15 October 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
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
print('\n' '----Plotting Comparison for SSW- %s----' % titletime)

### Call arguments
directorydata = '/home/zlabe/Documents/Research/StratoVari/Data/'
directoryfigure = '/home/zlabe/Desktop/'

### Read in data
ssw_waccm = np.genfromtxt(directorydata + 'WACCM_SSWClimo.txt')
ssw_e3sm = np.genfromtxt(directorydata + 'E3SM_SSWClimo.txt')
lev = np.genfromtxt(directorydata + 'Levels.txt')
diff = ssw_waccm - ssw_e3sm
time = np.arange(diff.shape[0])

dataall = [ssw_waccm,ssw_e3sm,diff]

###############################################################################
###############################################################################
###############################################################################
### Call arguments

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-100,101,5)
barlim = np.arange(-100,101,50)
    
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
        
cs = plt.contourf(time,lev,diff.transpose(),limit,extend='both')
cs1 = plt.contour(time,lev,diff.transpose(),np.arange(-900,901,25),
                  extend='both',linewidths=0.4,colors='dimgrey')

plt.axvline(x=31,linestyle='--',dashes=(1,0.3),color='k',linewidth=2)

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.xticks(np.arange(0,121,30),map(str,np.arange(0,121,30)),fontsize=5)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=5)
plt.minorticks_off()

plt.xlabel(r'\textbf{Days}',fontsize=8)
plt.ylabel(r'\textbf{Pressure (hPa)}',fontsize=8)
plt.title(r'\textbf{SC-WACCM4-1.1 minus E3SM-1.1}',fontsize=11)

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

plt.savefig(directoryfigure + 'SSWEvents_Difference.png',dpi=300)
print('Completed: Script done!')