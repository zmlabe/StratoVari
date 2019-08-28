"""
Scripts looks at resampling methods to understanding the number of ensembles
for statistically significant responses

Notes
-----
    Author : Zachary Labe
    Date   : 16 August 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyData as MO
import string

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Map Comparison- %s----' % titletime)

### Alott time series (300 ensemble members)
year1 = 1701
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
directoryfigure = '/home/zlabe/Desktop/STRATOVARI/'
varnames = ['U10','Z50','U200','Z500','SLP','THICK','T2M','RNET']
experi = np.array([r'\textbf{$\bf{\Delta}$Pi}',r'\textbf{$\bf{\Delta}$Cu}'])
letters = list(string.ascii_lowercase)

### Functions
def readPolarCapData(varnames,level,levq,sliceq,period):
    """
    Read in all data for polar cap (>65N)
    """
    if varnames == 'U10': # data needed to calculate polar vortex definition
        lat,lon,lev,varf = MO.readExperiAll(varnames,'Future','surface')
        lat,lon,lev,varcu = MO.readExperiAll(varnames,'Current','surface')
        lat,lon,lev,varpi = MO.readExperiAll(varnames,'Past','surface')
    else: # averaged over the polar cap
        lat,lon,lev,varf = MO.readExperiAveVar(varnames,'Future','polar',level)
        lat,lon,lev,varcu = MO.readExperiAveVar(varnames,'Current','polar',level)
        lat,lon,lev,varpi = MO.readExperiAveVar(varnames,'Past','polar',level)
    
    if sliceq == True: # data from the polar cap
        print('\nSlicing grid at level=%s!' % levq)
        levqq = np.where(lev == levq)[0]
        varf = varf[:,:,levqq].squeeze()
        varcu = varcu[:,:,levqq].squeeze()
        varpi = varpi[:,:,levqq].squeeze()
    elif varnames == 'U10':
        latq = np.where((lat >= 59.5) & (lat <= 60.5))[0]
        varf = np.nanmean(varf[:,:,latq,:].squeeze(),axis=2)
        varcu = np.nanmean(varcu[:,:,latq,:].squeeze(),axis=2)
        varpi = np.nanmean(varpi[:,:,latq,:].squeeze(),axis=2)
    else:
        print('\nSurface variable!')
        
    ### Calculate time mean [returns 1D array]
    print('-------Period of time: %s!-------\n' % period)
    if period == 'NDJFM':
        varwf = np.nanmean(np.append(varf[:,-2:],varf[:,:3],axis=1),axis=1)
        varwcu = np.nanmean(np.append(varcu[:,-2:],varcu[:,:3],axis=1),axis=1)
        varwpi = np.nanmean(np.append(varpi[:,-2:],varpi[:,:3],axis=1),axis=1)
    elif period == 'DJF':
        varwf = np.nanmean(np.append(varf[:,-1:],varf[:,:2],axis=1),axis=1)
        varwcu = np.nanmean(np.append(varcu[:,-1:],varcu[:,:2],axis=1),axis=1)
        varwpi = np.nanmean(np.append(varpi[:,-1:],varpi[:,:2],axis=1),axis=1)
    elif period == 'JFM':
        varwf = np.nanmean(varf[:,0:3],axis=1)
        varwcu = np.nanmean(varcu[:,0:3],axis=1)
        varwpi = np.nanmean(varpi[:,0:3],axis=1)
    elif period == 'JF':
        varwf = np.nanmean(varf[:,0:2],axis=1)
        varwcu = np.nanmean(varcu[:,0:2],axis=1)
        varwpi = np.nanmean(varpi[:,0:2],axis=1)
    elif period == 'OND':
        varwf = np.nanmean(varf[:,-3:],axis=1)
        varwcu = np.nanmean(varcu[:,-3:],axis=1)
        varwpi = np.nanmean(varpi[:,-3:],axis=1)
    elif period == 'ND':
        varwf = np.nanmean(varf[:,-2:],axis=1)
        varwcu = np.nanmean(varcu[:,-2:],axis=1)
        varwpi = np.nanmean(varpi[:,-2:],axis=1)
    elif period == 'MAM':
        varwf = np.nanmean(varf[:,2:5],axis=1)
        varwcu = np.nanmean(varcu[:,2:5],axis=1)
        varwpi = np.nanmean(varpi[:,2:5],axis=1)
    elif period == 'Annual':
        varwf = np.nanmean(varf[:,:],axis=1)
        varwcu = np.nanmean(varcu[:,:],axis=1)
        varwpi = np.nanmean(varpi[:,:],axis=1)
    elif period == 'Dec':
        varwf = varf[:,-1]
        varwcu = varcu[:,-1]
        varwpi = varpi[:,-1]
    elif period == 'Jan':
        varwf = varf[:,0]
        varwcu = varcu[:,0]
        varwpi = varpi[:,0]
    elif period == 'Feb':
        varwf = varf[:,1]
        varwcu = varcu[:,1]
        varwpi = varpi[:,1]
    elif period == 'Mar':
        varwf = varf[:,2]
        varwcu = varcu[:,2]
        varwpi = varpi[:,2]
    
    ### Check for missing data! this occurs between 100,200,300 ensembles
    varwf[np.where(varwf < -1e20)] = np.nan
    varwcu[np.where(varwcu < -1e20)] = np.nan
    varwpi[np.where(varwpi < -1e20)] = np.nan
    
    ### Calculate anomalies [1D arrays]
    anomcu = varwf - varwcu
    anompi = varwf - varwpi
    
    return lat,lon,anomcu,anompi

def calcSubEns(size,var):
    """ 
    Calculate random subsamples from large ensemble
    """
    subens = np.empty((var.shape[0]))
    subens[:] = np.nan
    for i in range(2,var.shape[0]): # min ensemble size of 2
        superens = np.random.choice(var,size=(size,i),replace=True)
        supermean = np.nanmean(superens,axis=1)
        subens[i] = abs(np.nanmax(supermean) - np.nanmin(supermean))
    print('Completed: calcSubEns function!')
    return subens

def readSubEns(variable,level,levq,levslice,period):
    """
    Process all data
    """
    lat,lon,anomcu,anompi = readPolarCapData(variable,level,levq,
                                             levslice,period)

    ### Perform combinations (100,000 random samples)
    subens_cu = calcSubEns(100000,anomcu)
    subens_pi = calcSubEns(100000,anompi)
    
    print('\n<<<<Completed: finished readings %s data combos!>>>>' % variable)
    return subens_cu,subens_pi

###############################################################################
###############################################################################
###############################################################################
### Read all functions
periodm = 'Annual'
pv_cu,pv_pi = readSubEns('U10','surface',None,False,periodm)
u10_cu,u10_pi = readSubEns('U','profile',10,True,periodm)     
z50_cu,z50_pi = readSubEns('GEOP','profile',50,True,periodm)     

slp_cu,slp_pi = readSubEns('SLP','surface',None,False,periodm)     
z500_cu,z500_pi = readSubEns('Z500','surface',None,False,periodm)     

t2_cu,t2_pi = readSubEns('T2M','surface',None,False,periodm)     
p_cu,p_pi = readSubEns('P','surface',None,False,periodm)     

### Create subplot files and arguement
suball = np.array([[u10_cu,u10_pi],[z50_cu,z50_pi],[slp_cu,slp_pi],
                   [z500_cu,z500_pi],[t2_cu,t2_pi],[p_cu,p_pi]])
xaxis = np.array([np.arange(0,301,50),np.arange(0,301,50),np.arange(0,301,50),
                  np.arange(0,301,50),np.arange(0,301,50),np.arange(0,301,50)])
yaxis = np.array([np.arange(0,51,10),np.arange(0,501,100),
                  np.arange(0,11,2),np.arange(0,201,25),
                  np.arange(0,8,1),np.round(np.arange(0,0.51,0.1),2)])
variablelist = np.array(['U10','Z50','SLP','Z500','T2M','P'])
ylabels = np.array([r'm/s',r'm',r'hPa',r'm',r'$\bf{^{\circ}}$C',r'mm/day'])
explist = np.array([[r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}'],
                    [r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}'],
                    [r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}'],
                    [r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}'],
                    [r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}'],
                    [r'{\textbf{$\bf{\Delta}$Cu}',r'\textbf{$\bf{\Delta}$Pi}']])
    
###############################################################################
###############################################################################
###############################################################################
### Create subplots of ensembles
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
        
fig = plt.figure(figsize=(8,5))
for i in range(len(suball)):
    ax = plt.subplot(2,3,i+1)
    
    if i < 3:
        adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('dimgrey')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(0)
        ax.tick_params('y',length=4,width=2,which='major',color='dimgrey',
                       pad=0.3)
        ax.tick_params('x',length=0,width=0,which='major',color='w',
                       pad=0,labelbottom=False)
    elif i >= 3:
        adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('dimgrey')
        ax.spines['bottom'].set_color('dimgrey')
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',
                       pad=0.3)
        
    plt.axvline(100,color='dimgrey',linestyle='--',dashes=(1,0.3),
                linewidth=2)    
    plt.hlines(suball[i][0][100],xmin=100,xmax=300,color='dimgrey',
              linestyle='--',linewidth=0.8)
    plt.plot(suball[i][0],linestyle='-',linewidth=0.8,color='rebeccapurple')
    plt.plot(suball[i][1],linestyle='-',linewidth=0.8,color='teal')
    
    if i == 0:
       plt.plot(pv_pi,linestyle='-',linewidth=0.8,color='darkorange')
       plt.plot(pv_cu,linestyle='-',linewidth=0.8,color='dodgerblue')
    
    plt.xticks(xaxis[i],list(map(str,xaxis[i])),fontsize=6)
    plt.xlim([0,xaxis[i].max()])
    plt.yticks(yaxis[i],list(map(str,yaxis[i])),fontsize=6)
    plt.ylim([0,yaxis[i].max()])
    
    if i == 4:
        plt.xlabel(r'\textbf{Number of Ensembles}',color='k',fontsize=10)
    ax.annotate(r'\textbf{[%s]}' % variablelist[i],xy=(0,0),
            xytext=(0.75,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=17)
    plt.ylabel(r'\textbf{Uncertainty [%s]}' % (ylabels[i]),
                         color='k',fontsize=10)
    
plt.tight_layout()
    
plt.savefig(directoryfigure + 'Ensemble_Subsampling_%s.png' % periodm,dpi=900)
    
    
    
    
    
    
    