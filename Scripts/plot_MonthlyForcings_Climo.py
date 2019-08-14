"""
Script plots monthly mean surface heat fluxes and sea ice extent in the Arctic

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2019
"""

### Import modules
import numpy as np
import datetime
import matplotlib.pyplot as plt
import cmocean

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/Forcing/'
directoryoutput = '/home/zlabe/Documents/Research/StratoVari/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plot Join SIE/RNET Figure- %s----' % titletime)

### Alott time series
months = np.arange(1,12+1,1)

### Read in sea ice extent data
sief = np.genfromtxt(directoryoutput + 'SIE/Arctic_SIE_PAMIP-1.6.txt')
siepd = np.genfromtxt(directoryoutput + 'SIE/Arctic_SIE_PAMIP-1.1.txt')
siepi = np.genfromtxt(directoryoutput + 'SIE/Arctic_SIE_PAMIP-1.5.txt')

### Calculate sea ice extent anomalies
icesat = sief - siepd
iceold = sief - siepi

### Read in heat flux data
heatsat = np.genfromtxt(directoryoutput + 'Arctic_RNETA_PD.txt')
heatold = np.genfromtxt(directoryoutput + 'Arctic_RNETA_PI.txt')

###############################################################################
###############################################################################
###############################################################################
### Create subplots of sea ice anomalies 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
    if 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        ax.yaxis.set_ticks([])
        
fig = plt.figure()
ax1 = plt.subplot(211)

#adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('w')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(0)
ax1.spines['bottom'].set_linewidth(2)
ax1.tick_params('both',length=4,width=2,which='major',color='dimgrey')

N = len(icesat)
ind = np.arange(N)
width = 0.9
rects = ax1.bar(ind,icesat*-1,zorder=1,edgecolor='deepskyblue',
                facecolor='deepskyblue')

plt.yticks(np.arange(0,12,1),list(map(str,np.arange(0,12,1))),fontsize=6)
plt.ylabel(r'\textbf{Sea Ice Extent [-1$\times$10$^6$ km$^2$]}',
                     color='k',fontsize=6,labelpad=5)
xlabels = [r'JAN',r'FEB',r'MAR',r'APR',r'MAY',r'JUN',
           r'JUL',r'AUG',r'SEP',r'OCT',r'NOV',r'DEC',]
plt.xticks(np.arange(0,30,1),xlabels,fontsize=8)
plt.xlim([0,11])
plt.ylim([0,6])

ax2 = ax1.twinx()
#adjust_spines(ax2, ['right','bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('dimgrey')
ax2.spines['left'].set_color('w')
ax2.spines['bottom'].set_color('dimgrey')
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(2)
ax2.spines['bottom'].set_linewidth(2)
ax2.tick_params('both',length=4,width=2,which='major',color='dimgrey')

ax2.plot(heatsat,color='crimson',linewidth=4,markersize=6,
         marker='o',clip_on=False,zorder=10)
plt.ylabel(r'\textbf{Net Surface Heat Flux [W/m$\bf{^{2}}$]}',
                     color='k',fontsize=6,labelpad=5)
plt.yticks(np.arange(-40,90,10),list(map(str,np.arange(-40,90,10))),fontsize=6)
plt.ylim([-10,40])

ax1.annotate(r'\textbf{[%s]}' % 'a',xy=(0,0),
        xytext=(0.02,0.94),xycoords='axes fraction',
        color='dimgrey',fontsize=9)

###############################################################################
ax1 = plt.subplot(212)

#adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('w')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(0)
ax1.spines['bottom'].set_linewidth(2)
ax1.tick_params('both',length=4,width=2,which='major',color='dimgrey')

N = len(icesat)
ind = np.arange(N)
width = 0.9
rects = ax1.bar(ind,iceold*-1,zorder=1,edgecolor='deepskyblue',
                facecolor='deepskyblue')

plt.yticks(np.arange(0,12,1),list(map(str,np.arange(0,12,1))),fontsize=6)
plt.ylabel(r'\textbf{Sea Ice Extent [-1$\times$10$^6$ km$^2$]}',
                     color='k',fontsize=6,labelpad=5)
xlabels = [r'JAN',r'FEB',r'MAR',r'APR',r'MAY',r'JUN',
           r'JUL',r'AUG',r'SEP',r'OCT',r'NOV',r'DEC',]
plt.xticks(np.arange(0,30,1),xlabels,fontsize=8)
plt.xlim([0,11])
plt.ylim([0,6])

ax2 = ax1.twinx()
#adjust_spines(ax2, ['right','bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('dimgrey')
ax2.spines['left'].set_color('w')
ax2.spines['bottom'].set_color('dimgrey')
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(2)
ax2.spines['bottom'].set_linewidth(2)
ax2.tick_params('both',length=4,width=2,which='major',color='dimgrey')

ax2.plot(heatold,color='crimson',linewidth=4,markersize=6,
         marker='o',clip_on=False,zorder=10)
plt.ylabel(r'\textbf{Net Surface Heat Flux [W/m$\bf{^{2}}$]}',
                     color='k',fontsize=6,labelpad=5)
plt.yticks(np.arange(-40,90,10),list(map(str,np.arange(-40,90,10))),fontsize=6)
plt.ylim([-10,40])

ax1.annotate(r'\textbf{[%s]}' % 'b',xy=(0,0),
        xytext=(0.02,0.94),xycoords='axes fraction',
        color='dimgrey',fontsize=9)

plt.subplots_adjust(hspace=0.3)
plt.tight_layout()

plt.savefig(directoryfigure + 'RNETSIE_Anomalies.png',dpi=900)

