"""
Script reads in monthly data from PAMIP experiments using SC-WACCM4 
for all 300 ensemble members! 
 
Notes
-----
    Author : Zachary Labe
    Date   : 21 June 2019
    
Usage
-----
    [1] readExperiAll(varid,timeperiod,level)
    [2] readExperiAveVar(varid,timeperiod,region,level)
"""

def readExperiAll(varid,timeperiod,level):
    """
    Function reads monthly data from PAMIP simulations for 300 members

    Parameters
    ----------
    varid : string
        variable name to read
    timeperiod: string
        time of analysis (Future, Current, Past)
    level : string
        Height of variable (surface or profile)
        

    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    var : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]

    Usage
    -----
    lat,lon,lev,var = readExperiAll(varid,timeperiod,level)
    """
    print('\n>>>>>>>>>> Using readExperiAll function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory (1-300 members)
    if any([timeperiod=='Future',timeperiod=='Current',timeperiod=='Past']):
        ### Directory for all ensemble members (remote server - seley)
        directorydata = '/seley/zlabe/simu/'
        
        if timeperiod == 'Future':
            experi = 'PAMIP_Fu'
        elif timeperiod == 'Current':
            experi = 'PAMIP_Cu'
        elif timeperiod == 'Past':
            experi = 'PAMIP_Pi'
        else:
            print(ValueError('Selected wrong time period (Future,Current,Past!')) 
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1701-2000.nc'
    
    ### Directories for thickness experiments (1-100 members)
    elif any([timeperiod=='SIT_Fu',timeperiod=='SIT_Cu']):
        if timeperiod == 'SIT_Fu':
            experi = 'PAMIP-1.10-300yr'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1700-2000.nc'
            print('-----------USING THICKNESS EXPERIMENTS (Future)!-----------')
        elif timeperiod == 'SIT_Cu':
            experi = 'PAMIP-1.9-300yr'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1700-2000.nc'
            print('-----------USING THICKNESS EXPERIMENTS (Present-Day)!-----------')
        else:
            print(ValueError('Selected wrong time period (SIT_Fu,SIT_Cu!')) 
    elif any([timeperiod=='E3SM_Fu',timeperiod=='E3SM_Cu',timeperiod=='E3SM_Pi']):
        if timeperiod == 'E3SM_Fu':
            experi = 'PAMIP-1.6-E3SM'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        elif timeperiod == 'E3SM_Cu':
            experi = 'PAMIP-1.1-E3SM'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        elif timeperiod == 'E3SM_Pi':
            experi = 'PAMIP-1.5-E3SM'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        else:
            print(ValueError('Selected wrong time period (E3SM_Fu,E3SM_Cu!')) 
    elif any([timeperiod=='Osea_Fu',timeperiod=='BKsea_Fu']):
        if timeperiod == 'Osea_Fu':
            experi = 'PAMIP-3.1'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('-----------USING WACCM OSeaIce EXPERIMENTS!-----------')
        elif timeperiod == 'BKsea_Fu':
            experi = 'PAMIP-3.2'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('-----------USING WACCM BKSeaIce EXPERIMENTS!-----------')
        else:
            print(ValueError('Selected wrong time period (Osea_Fu,BKsea_Fu!')) 
    else:
        print(ValueError('Selected wrong experiment name!'))
    
    if varid == 'EGR' and level == 'surface': # integrated from 500-850 hPa
        filename = totaldirectory + varid + '_500_850.nc'

    ### Read in Data
    if level == 'surface': # 3d variables
        data = Dataset(filename,'r')
        lev = 'surface'
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 4d variables
        data = Dataset(filename,'r')
        lev = data.variables['level'][:]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'zonmean': # 3d variables (zonal mean!)
        varidz = varid + '_' + level
        filename = totaldirectory + varidz + '_1701-2000.nc'
        if any([timeperiod=='BKsea_Fu',timeperiod=='Osea_Fu',
                timeperiod=='E3SM_Cu',timeperiod=='E3SM_Pi',
                timeperiod=='SIT_Fu',timeperiod=='SIT_Cu']):
            filename = totaldirectory + varidz + '_1900-2000.nc'
        elif timeperiod=='E3SM_Fu':
            filename = totaldirectory + varidz + '_1910-2000.nc'
        data = Dataset(filename,'r')
        lev = data.variables['level'][:]
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        varq = data.variables['%s' % varid][:].squeeze()
        data.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:],varid))

    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(varq.shape[0]//months,months,
                              int(lat.shape[0]),int(lon.shape[0])))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(varq.shape[0]//months,months,int(lev.shape[0]),
                      int(lat.shape[0]),int(lon.shape[0])))
    elif level == 'zonmean': # 3d variables (zonal mean!)
        var = np.reshape(varq,(varq.shape[0]//months,months,int(lev.shape[0]),
                      int(lat.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif varid == 'SWE':
        var = var*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')
        
    print('Completed: Read members 1-300!')

    print('>>>>>>>>>> Completed: Finished readExperiAll function!')
    return lat,lon,lev,var

###############################################################################
###############################################################################

def readExperiAveVar(varid,timeperiod,region,level):
    """
    Function reads weighted averaged monthly data from PAMIP simulations
    for 300 members

    Parameters
    ----------
    varid : string
        variable name to read
    timeperiod: string
        time of analysis (Future, Current, Past)
    region : string
        Area where the weighted average took place
    level : string
        Height of variable (surface or profile)
        

    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    var : 2d numpy array or 3d numpy array 
        [year,month] or [year,month,level]

    Usage
    -----
    lat,lon,lev,var = readExperiAveVar(varid,timeperiod,region)
    """
    print('\n>>>>>>>>>> Using readExperiAveVar function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory (1-300 members)
    if any([timeperiod=='Future',timeperiod=='Current',timeperiod=='Past']):
        ### Directory for all ensemble members (remote server - seley)
        directorydata = '/seley/zlabe/simu/'
        
        if timeperiod == 'Future':
            experi = 'PAMIP_Fu'
        elif timeperiod == 'Current':
            experi = 'PAMIP_Cu'
        elif timeperiod == 'Past':
            experi = 'PAMIP_Pi'
        else:
            print(ValueError('Selected wrong time period (Future,Current,Past!')) 
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_' + region + 'mean' + '_1701-2000.nc'
    
    ### Directories for thickness experiments (1-100 members)
    elif any([timeperiod=='SIT_Fu',timeperiod=='SIT_Cu']):
        if timeperiod == 'SIT_Fu':
            experi = 'PAMIP-1.10'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING THICKNESS EXPERIMENTS!-----------')
        elif timeperiod == 'SIT_Cu':
            experi = 'PAMIP-1.9'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING THICKNESS EXPERIMENTS!-----------')
        else:
            print(ValueError('Selected wrong time period (SIT_Fu,SIT_Cu!')) 
    elif any([timeperiod=='E3SM_Fu',timeperiod=='E3SM_Cu',timeperiod=='E3SM_Pi']):
        if timeperiod == 'E3SM_Fu':
            experi = 'PAMIP-1.6-E3SM'
            directorydata = '/seley/zlabe/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        elif timeperiod == 'E3SM_Cu':
            experi = 'PAMIP-1.1-E3SM'
            directorydata = '/seley/zlabe/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        elif timeperiod == 'E3SM_Pi':
            experi = 'PAMIP-1.5-E3SM'
            directorydata = '/seley/zlabe/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING DOE E3SM EXPERIMENTS!-----------')
        else:
            print(ValueError('Selected wrong time period (E3SM_Fu,E3SM_Cu!')) 
    elif any([timeperiod=='Osea_Fu',timeperiod=='BKsea_Fu']):
        if timeperiod == 'Osea_Fu':
            experi = 'PAMIP-3.1'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING WACCM OSeaIce EXPERIMENTS!-----------')
        elif timeperiod == 'BKsea_Fu':
            experi = 'PAMIP-3.2'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_' + region + 'mean' + '_1900-2000.nc'
            print('-----------USING WACCM BKSeaIce EXPERIMENTS!-----------')
        else:
            print(ValueError('Selected wrong time period (Osea_Fu,BKsea_Fu!')) 
    else:
        print(ValueError('Selected wrong experiment name!'))

    ### Read in Data
    if level == 'surface': # 1d variables
        data = Dataset(filename,'r')
        lev = 'surface'
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 2d variables
        data = Dataset(filename,'r')
        lev = data.variables['level'][:]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:,:] 
        data.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:],varid))

    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(varq.shape[0]//months,months))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(varq.shape[0]//months,months,int(lev.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped *%s_%smean* array!' % (varid,region))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif varid == 'SWE':
        var = var*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')
        
    print('Completed: Read members 1-300!')

    print('>>>>>>>>>> Completed: Finished readExperiAveVar function!')
    return lat,lon,lev,var

#### Test function -- no need to use    
#varid = 'T2M'
#timeperiod = 'E3SM_Pi'
#level = 'surface'
#lat,lon,lev,var = readExperiAll(varid,timeperiod,level)
    
#varid = 'TEMP'
#timeperiod = 'Past'
#level = 'profile'
#region = 'polar'
#lat,lon,lev,var = readExperiAveVar(varid,timeperiod,region,level)