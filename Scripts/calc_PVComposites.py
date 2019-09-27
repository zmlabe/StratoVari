"""
Function to slice polar vortex years in each experiment

Notes
-----
    Author : Zachary Labe
    Date   : 27 September 2019
"""

def polarVortexComp(simu,varname,vartype,stat,time):
    """
    Calculate polar vortex statistics per year
    """
    
    ### Import modules
    import numpy as np
    import read_MonthlyData as MO

    ### Call function for 4d variable data
    lat,lon,lev,var= MO.readExperiAll(varname,simu,vartype)
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data 
    var[np.where(var < -1e10)] = np.nan
    
    ### Calculate polar vortex strength using 60N
    if varname == 'U10':
        latq = np.where((lat >= 59.5) & (lat <= 60.5))[0]
        varu = var[:,:,latq,:].squeeze()
        
    ### Calculate time month mean [ensemble,month,longitude]
    if time == 'annual':
        varf = np.nanmean(varu[:,:,:],axis=2)
    elif time == 'JFM':
        varf = np.nanmean(varu[:,0:3,:],axis=2)
    elif time == 'JFMA':
        varf = np.nanmean(varu[:,0:4,:],axis=2)
    elif time == 'JF':
        varf = np.nanmean(varu[:,0:2,:],axis=2)
    elif time == 'J':
        varf = np.nanmean(varu[:,0:1,:],axis=2)
    elif time == 'FM':
        varf = np.nanmean(varu[:,1:3,:],axis=2)
    elif time == 'D':
        varf = np.nanmean(varu[:,-1:,:],axis=2)
    elif time == 'NDJFMA':
        varf= np.append(varu[:,-2:,:],varu[:,:4,:],axis=1)

    ### Calculate zonal mean [ensemble,longitude]
    if time == 'NDJFMA':
        pvf = np.nanmean(varf,axis=2)
    else:
        pvf = np.nanmean(varf,axis=1)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Calculate statistics
    if time == 'NDJFMA':
        stats = []
        for i in range(pvf.shape[1]):
            if stat == '5-95':
                pvf5 = np.nanpercentile(pvf[:,i],5)
                pvf95 = np.nanpercentile(pvf[:,i],95)
                argyearsq = np.where((pvf[:,i] > pvf5) & (pvf[:,i] < pvf95))[0]
                stats.append(argyearsq)
            elif stat == '10-90':
                pvf10 = np.nanpercentile(pvf[:,i],10)
                pvf90 = np.nanpercentile(pvf[:,i],90)
                argyearsq = np.where((pvf[:,i] > pvf10) & (pvf[:,i] < pvf90))[0]
                stats.append(argyearsq)
        argyears = np.asarray(stats)
    else:
        if stat == '5-95':
            pvf5 = np.nanpercentile(pvf,5)
            pvf95 = np.nanpercentile(pvf,95)
            argyears = np.where((pvf > pvf5) & (pvf < pvf95))[0]
        elif stat == '10-90':
            pvf10 = np.nanpercentile(pvf,10)
            pvf90 = np.nanpercentile(pvf,90)
            argyears = np.where((pvf > pvf10) & (pvf < pvf90))[0]
        elif stat == '33-66':
            pvf33 = np.nanpercentile(pvf,33.333)
            pvf66 = np.nanpercentile(pvf,66.666)
            argyears = np.where((pvf > pvf33) & (pvf < pvf66))[0]
        elif stat == '>50':
            pvf50 = np.nanpercentile(pvf,50)
            argyears = np.where((pvf > pvf50))[0]
        elif stat == '<50':
            pvf50 = np.nanpercentile(pvf,50)
            argyears = np.where((pvf < pvf50))[0]
        elif stat == '>66':
            pvf66 = np.nanpercentile(pvf,66.666)
            argyears = np.where((pvf > pvf66))[0]
        elif stat == '<33':
            pvf33 = np.nanpercentile(pvf,33.333)
            argyears = np.where((pvf < pvf33))[0]
        elif stat == '1sigma':
            stdev = np.nanstd(pvf)
            stdl = np.nanmean(pvf) - stdev
            stdu = np.nanmean(pvf) + stdev
            argyears = np.where((pvf > stdl) & (pvf < stdu))[0]
    
    return argyears

###############################################################################
###############################################################################
###############################################################################
### Test functions (not needed!!!)
#argyears = polarVortexComp('Future','U10','surface','<50','J')