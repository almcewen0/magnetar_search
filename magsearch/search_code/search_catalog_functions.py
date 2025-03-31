from search_code.search_config import config
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from search_code.search_convenience_functions import *
import numpy as np



def read_catalog(c_cut):
    """
    Read in the catalog and make cuts as described below.
    """
    cf = config()
    cat = Table.read(cf.catalog_file)
    s_cat = SkyCoord(np.array(cat['LII'])*u.deg, np.array(cat['BII'])*u.deg,frame='galactic')
    
    b_cut = 5*u.deg
    
    
    # abridged version has the following changes:
    #    - omitted detections outside of |b| <= 5 deg
    #    - omitted the SmallWindow observations
    #    - omitted sources with PN_8_CTS < c_cut (including sources without a PN_8_CTS measurement
    #    - added column for Hui's IDs
    #    - added column for SIMBAD associations (position matching within 2*uncertainty of 4XMM position)
    #    - excluded sources with SUM_FLAG >= 3
 
    
    cat.add_column(cat['MJD_STOP'] - cat['MJD_START'],name='TOBS')
    abr_cnd = np.array(abs(s_cat.b.to('deg')) <= b_cut) \
            & np.array(~cat['PN_8_CTS'].mask) \
            & np.array((cat['PN_8_CTS'] >= c_cut).data) \
            & np.array(cat['PN_SUBMODE'] != 'PrimeSmallWindow') \
            & np.array(cat['SUM_FLAG'] <= 3)
    
    cat_abr = cat[abr_cnd]
    return cat,cat_abr
    
def add_SN_columns(subcat):
    """
    Function to calculate S/N for different combinations of bands in 4XMM.
    Calculates S/N for each band individually, then (optionally) steps through each *contiguous* combination.
    """
    bands = [1,2,3,4,5,8]
    if 'SN_1' not in subcat.columns:
        [subcat.add_column(calc_joint_sn(subcat[f'PN_{i}_RATE'],
                                         subcat[f'PN_{i}_RATE_ERR'],
                                         subcat['TOBS']
                                        ),name=f"SN_{i}") for i in bands]
    colnames = np.array(subcat.colnames)[['SN_' in v for v in list(subcat.colnames)]]
    return subcat,colnames

def calc_joint_sn(rates,rate_errs,tobs):
    if len(np.shape(rates)) > 1:
        return np.sum(np.array(rates)*np.sqrt(tobs),axis=0)/np.sqrt(np.sum(np.array(rate_errs)),axis=0)
    else:
        return np.sum([np.array(rates)*np.sqrt(tobs)],axis=0)/np.sqrt(np.sum([np.array(rate_errs)],axis=0))

def optimal_bands(ln):
    """
    Function that utilizes combine_factors() to return the best energy ranges for the search. This will (hopefully)
    optimize the H power of a periodic signal in the data.
    EDIT 041724: No longer utilizes combine_factors, and just uses band with highest S/N.
    """

    joint_bands = np.array(["1","2","3","4","5",'1-2','2-3','3-4','4-5',
                            '1-2-3','2-3-4','3-4-5','1-2-3-4','2-3-4-5','1-2-3-4-5'])
    exposure = ln['TOBS'] * u.day.to('s')
    snrs = np.zeros(len(joint_bands))
    for i,combination in enumerate(joint_bands):
        bands_to_check = combination.split('-')
        if len(bands_to_check) == 1:
            rates  = ln[f'PN_{bands_to_check[0]}_RATE']
            drates = ln[f'PN_{bands_to_check[0]}_RATE_ERR']
        else:
            rates  = np.sum([ln[f'PN_{b}_RATE'] for b in bands_to_check])
            drates = np.sqrt(np.sum([ln[f'PN_{b}_RATE_ERR']**2 for b in bands_to_check]))
        snrs[i] = np.sum(rates)*np.sqrt(exposure)/np.sqrt(np.sum(drates))
    best_bandints = joint_bands[snrs == snrs.max()]
    if '-' in best_bandints[0]:
        bs = np.array(best_bandints[0].split('-')).astype(int)
        if bs.min() <= 3 and bs.max() > 3:
            b1 = '-'.join([str(b) for b in range(bs.min(),4)])
            b2 = '-'.join([str(b) for b in range(4,bs.max()+1)])
            best_bandints = [b1,b2]

    return [ints_to_band(bb) for bb in best_bandints]

