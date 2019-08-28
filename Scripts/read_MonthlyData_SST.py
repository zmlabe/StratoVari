"""
Script reads in monthly data from PAMIP experiments using SC-WACCM4 
for all 100 ensemble members! 
 
Notes
-----
    Author : Zachary Labe
    Date   : 12 August 2019
    
Usage
-----
    readExperiAll(directory,varid,timeperiod,level)
"""

def readExperiAll(varid,timeperiod,level):
    """
    Function reads monthly data from PAMIP simulations for 100 members

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
    ### Directories for SST experiments (1-100 members)
    if any([timeperiod=='SST_Fu',timeperiod=='SST_Cu',timeperiod=='SST_Pi']):
        if timeperiod == 'SST_Fu':
            experi = 'PAMIP-1.4'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('Reading in Future Sea Surface Temperatures!')
        elif timeperiod == 'SST_Cu':
            experi = 'PAMIP-1.1-QBO'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('Reading in Present-Day Sea Surface Temperatures!')
        elif timeperiod == 'SST_Pi':
            experi = 'PAMIP-1.3'
            directorydata = '/seley/ypeings/simu/'
            totaldirectory = directorydata + experi + '/monthly/'
            filename = totaldirectory + varid + '_1900-2000.nc'
            print('Reading in Pre-Industrial Sea Surface Temperatures!')
        else:
            print(ValueError('Selected wrong time period!')) 
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
        filename = totaldirectory + varidz + '_1900-2000.nc'
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
        
    print('Completed: Read members 1-100!')

    print('>>>>>>>>>> Completed: Finished readExperiAll function!')
    return lat,lon,lev,var

#### Test function -- no need to use    
#varid = 'U'
#timeperiod = 'SST_Fu'
#level = 'profile'
#lat,lon,lev,var = readExperiAll(varid,timeperiod,level)