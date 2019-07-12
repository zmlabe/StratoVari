"""
Functions to calculate minimum number of ensembles for 2d maps
 
Notes
-----
    Author : Zachary Labe
    Date   : 27 June 2019
    
Usage
-----
    [1] computeMean(datain,typemean)
    [2] computeSTD(datain,df)
    [3] computePooledSD(xstd,ystd,xn,yn)
    [4] computeMinEns(future,climo,alpha)
"""
def computeMean(datain,typemean):
    """
    Compute simple means for different dimensions of array
    """
    print('>>>>>>>> Using computeMean function!')
    
    ### Import modules
    import numpy as np
    
    ### Calculate various means
    if typemean == 'zonal': # zonal mean
        dataout = np.nanmean(datain,axis=datain.ndim-1)
        print('Completed: zonal mean!')
    elif typemean == 'ensemble': # ensemble mean
        dataout = np.nanmean(datain,axis=0)
        print('Completed: ensemble mean!')
     
    print('>>>>>>>> Ending computeMean function!\n')    
    return dataout
###############################################################################
###############################################################################
###############################################################################
def computeSTD(datain,df):
    """
    Compute standard deviation of ensemble members
    """
    print('>>>>>>>> Using computeSTD function!')
    
    ### Import modules
    import numpy as np
    
    ### Compute standard deviation
    dataout = np.nanstd(datain,axis=0,ddof=df,
                        dtype=np.float64) # calculate for ensemble members
    
    print('>>>> Ending computeSTD function!\n')
    return dataout
###############################################################################
###############################################################################
###############################################################################
def computePooledSD(xstd,ystd,xn,yn):
    """
    Compute pooled standard deviation
    """
    print('>>>>>>>> Using computePooledSD function!')
    
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
    elif xstd.ndim == 3:
        sp = np.empty((xstd.shape))
        for mo in range(xstd.shape[0]):
            for i in range(xstd.shape[1]):
                for j in range(xstd.shape[2]):
                    sp[mo,i,j] = np.sqrt(((xn - 1)*xstd[mo,i,j]**2 + \
                              (yn - 1)*ystd[mo,i,j]**2) \
                              /(xn + yn - 2))
    
    print('>>>>>>>> Ending computePooledSD function!\n')
    return sp
###############################################################################
###############################################################################
###############################################################################
def computeMinEns(future,climo,alpha):
    """
    Compute minimum ensemble number using formula 4 from 
    Screen et al. 2013, Climate Dynamics
    """
    print('>>>>>>>> Using computeMinEns function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    futurem = computeMean(future,'ensemble')
    futurestd = computeSTD(future,1)
    
    climom = computeMean(climo,'ensemble')
    climostd = computeSTD(climo,1)
    
    ### Calculate t statistic for confidence level
    tc = sts.t.ppf(1-(alpha/2),len(future)-1) # two-tailed
    
    ### Calculate pooled standard deviation
    sp = computePooledSD(futurestd,climostd,len(future),len(climo))
    
    ### Compute minimum ensemble number
    nmin = (2*(tc**2)) * (sp/(futurem - climom))**2
    
    nmin[np.where(nmin >= len(future))] = np.nan
    
    print('>>>>>>>> Ending computeMinEns function!\n')
    return nmin
###############################################################################
###############################################################################
###############################################################################