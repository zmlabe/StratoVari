"""
Function to slice polar vortex response years

Notes
-----
    Author : Zachary Labe
    Date   : 24 September 2019
"""

def polarVortexStats(simuF,simuP,varname,vartype,stat,time):
    """
    Calculate polar vortex statistics per year
    """
    
    ### Import modules
    import numpy as np
    import read_MonthlyData as MO

    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varname,simuF,vartype)
    lat,lon,lev,varpast = MO.readExperiAll(varname,simuP,vartype)
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data 
    varfuture[np.where(varfuture < -1e10)] = np.nan
    varpast[np.where(varpast < -1e10)] = np.nan
    
    ### Calculate polar vortex strength using 60N
    if varname == 'U10':
        latq = np.where((lat >= 59.5) & (lat <= 60.5))[0]
        latu = lat[latq].squeeze()
        varfutureu = varfuture[:,:,latq,:].squeeze()
        varpastu = varpast[:,:,latq,:].squeeze()
        
    ### Calculate time month mean [ensemble,month,longitude]
    if time == 'annual':
        varf = np.nanmean(varfutureu[:,:,:],axis=2)
        varp = np.nanmean(varpastu[:,:,:],axis=2)
    elif time == 'JFM':
        varf = np.nanmean(varfutureu[:,0:3,:],axis=2)
        varp = np.nanmean(varpastu[:,0:3,:],axis=2)
    elif time == 'JF':
        varf = np.nanmean(varfutureu[:,0:2,:],axis=2)
        varp = np.nanmean(varpastu[:,0:2,:],axis=2)
    elif time == 'J':
        varf = np.nanmean(varfutureu[:,0:1,:],axis=2)
        varp = np.nanmean(varpastu[:,0:1,:],axis=2)
    elif time == 'FM':
        varf = np.nanmean(varfutureu[:,1:3,:],axis=2)
        varp = np.nanmean(varpastu[:,1:3,:],axis=2)
    elif time == 'D':
        varf = np.nanmean(varfutureu[:,-1:,:],axis=2)
        varp = np.nanmean(varpastu[:,-1:,:],axis=2)
    elif time == 'NDJFMA':
        varf= np.append(varfutureu[:,-2:,:],varfutureu[:,:4,:],axis=1)
        varp = np.append(varpastu[:,-2:,:],varpastu[:,:4,:],axis=1)

    ### Calculate zonal mean [ensemble,longitude]
    if time == 'NDJFMA':
        pvf = np.nanmean(varf,axis=2)
        pvp = np.nanmean(varp,axis=2)
    else:
        pvf = np.nanmean(varf,axis=1)
        pvp = np.nanmean(varp,axis=1)
    
    ### Calculate anomalies
    anom = pvf - pvp 

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Calculate statistics
    if time == 'NDJFMA':
        stats = []
        for i in range(anom.shape[1]):
            if stat == '5-95':
                anom5 = np.nanpercentile(anom[:,i],5)
                anom95 = np.nanpercentile(anom[:,i],95)
                argyearsq = np.where((anom[:,i] > anom5) & (anom[:,i] < anom95))[0]
                stats.append(argyearsq)
            elif stat == '10-90':
                anom10 = np.nanpercentile(anom[:,i],10)
                anom90 = np.nanpercentile(anom[:,i],90)
                argyearsq = np.where((anom[:,i] > anom10) & (anom[:,i] < anom90))[0]
                stats.append(argyearsq)
        argyears = np.asarray(stats)
    else:
        if stat == '5-95':
            anom5 = np.nanpercentile(anom,5)
            anom95 = np.nanpercentile(anom,95)
            argyears = np.where((anom > anom5) & (anom < anom95))[0]
        elif stat == '10-90':
            anom10 = np.nanpercentile(anom,10)
            anom90 = np.nanpercentile(anom,90)
            argyears = np.where((anom > anom10) & (anom < anom90))[0]
        elif stat == '33-66':
            anom33 = np.nanpercentile(anom,33.333)
            anom66 = np.nanpercentile(anom,66.666)
            argyears = np.where((anom > anom33) & (anom < anom66))[0]
        elif stat == '>50':
            anom50 = np.nanpercentile(anom,50)
            argyears = np.where((anom > anom50))[0]
        elif stat == '<50':
            anom50 = np.nanpercentile(anom,50)
            argyears = np.where((anom < anom50))[0]
        elif stat == '>66':
            anom66 = np.nanpercentile(anom,66.66)
            argyears = np.where((anom > anom66))[0]
        elif stat == '<33':
            anom33 = np.nanpercentile(anom,33.333)
            argyears = np.where((anom < anom33))[0]
        elif stat == '1sigma':
            stdev = np.nanstd(anom)
            stdl = np.nanmean(anom) - stdev
            stdu = np.nanmean(anom) + stdev
            argyears = np.where((anom > stdl) & (anom < stdu))[0]
    
    return argyears

###############################################################################
###############################################################################
###############################################################################
### Test functions (not needed!!!)
#argyears = polarVortexStats('Future','Past','U10','surface','10/90','NDJFMA')