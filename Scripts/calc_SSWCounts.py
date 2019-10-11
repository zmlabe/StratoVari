"""
Count total number of SSW in NDJFM

Notes
-----
    Author : Zachary Labe
    Date   : 11 October 2019
"""

def calc_SSWCount(experi,time):
    """
    Calculate total sudden stratospheric warming event count
    """
    import numpy as np
    
    ### Select experiment path
    if experi == 'E3SM_Cu':
        directorydata = '/seley/zlabe/simu/PAMIP-1.1-E3SM/daily/'
        directorydata2 = '/home/zlabe/Documents/Research/StratoVari/Data/'
        expn = 100
        
        if time == 'NDJFM':
            ND = np.empty((expn))
            JFM = np.empty((expn))
            for i in range(expn):
                filenameND = 'PAMIP-1.1-E3SM_%s/sswcount_%s_11-12.txt' % (i+1,i+1)
                filenameJFM = 'PAMIP-1.1-E3SM_%s/sswcount_%s_1-3.txt' % (i+1,i+1)
                ND[i] = np.genfromtxt(directorydata + filenameND,unpack=True)
                JFM[i] = np.genfromtxt(directorydata + filenameJFM,unpack=True)
                
            ### Tally SSW totals
            count = ND + JFM
            freq = np.nansum(count)/10.
            
            ### Save output
            np.savetxt(directorydata2 + 'NDJFM_SSW_%s' % experi,count)
###############################################################################
###############################################################################
###############################################################################            
    elif experi == 'PAMIP_Cu':
        directorydata = '/seley/zlabe/simu/PAMIP_Cu/daily/'
        directorydata2 = '/home/zlabe/Documents/Research/StratoVari/Data/'
        expn = 300
        
        if time == 'NDJFM':
            ND = np.empty((expn))
            JFM = np.empty((expn))
            for i in range(1,expn,1):
                filenameND = 'PAMIP-1.1-QBO_%s/sswcount_%s_11-12.txt' % (i+1,i+1)
                filenameJFM = 'PAMIP-1.1-QBO_%s/sswcount_%s_1-3.txt' % (i+1,i+1)
                ND[i] = np.genfromtxt(directorydata + filenameND,unpack=True)
                JFM[i] = np.genfromtxt(directorydata + filenameJFM,unpack=True)
                
            ### Tally SSW totals
            count = ND + JFM
            freq = np.nansum(count)/10.
            
            ### Save output
            np.savetxt(directorydata2 + 'NDJFM_SSW_%s' % experi,count)
###############################################################################
###############################################################################
###############################################################################
    print('SSW Frequency - %s per decade!' % freq)        
    return count,freq
            
### Test functions (do not use!)
#count,freq = calc_SSWCount('E3SM_Cu','NDJFM')
count,freq = calc_SSWCount('PAMIP_Cu','NDJFM')
            
        
