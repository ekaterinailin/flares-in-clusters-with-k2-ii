import pandas as pd
import numpy as np
import astropy.units as u
import pickle
import datetime
from warnings import warn

from .cmd import correct_for_extinction

from .__init__ import logger

def timestamp():
    '''
    Get immediate timestamp string in
    YYYY_MM_DD_hh format.
    '''
    return datetime.datetime.now().strftime("%Y_%m_%d_%H")

def _calculate_observation_time(s, f, n_duplicates=False):
    '''
    Input hours return yr.

    Parameters:
    ------------
    s : DataFrame
         OC.stars type DataFrame
    f : DataFrame
         OC.flares type DataFrame
    n_duplicates : False or bool
        return number of dropped duplicate rows

    Return:
    --------
    astropy Quantity : total observation time in df
    '''
    #safety check, that there are no duplicates of light curves for same target in the same campaign
    fd = f.drop_duplicates(subset=['EPIC','Campaign'])
    sd = s.drop_duplicates(subset=['EPIC','Campaign'])
    campaignobs = fd.groupby('Campaign').total_obs_time.mean() # mean observing time in each campaign
    sd['t'] = np.nan
    for c in list(set(sd.Campaign.values)):
        try:
            sd.loc[sd.Campaign==c,'t'] = campaignobs[c] #map mean obs. times on stars where light curves were not inspected
        except KeyError:
            logger.warning('This is just a workaround until all LCs are filed in the flare table, using NaNs if no flares occurred.')
            sd.loc[sd.Campaign==c,'t'] = np.mean(campaignobs.values)#use the mean of all othe available campaigns as a workaround
    t = sd['t'].sum()#sum up the all the light curves
    if np.isnan(fd.total_obs_time.values).all():
        logger.warn('No observation times given for flare targets.')
    if n_duplicates == True:
        return t/24./365.25, f.shape[0]-fd.shape[0]
    else:
        return t/24./365.25

def set_plotlimits(ax, a, freq, xlim=None, ylim=None):
    '''
    Set plot limits according to the data or given
    limits.
    '''
    if xlim == None:
        ax.set_xlim([a.min()*.5, a.max()*1.5])
    else:
        ax.set_xlim(xlim)
    if ylim == None:
        ax.set_ylim([freq.min()*.5, freq.max()*1.5])
    else:
        ax.set_ylim(ylim)
    return

def generate_random_power_law_distribution(a, b, g, size=1, seed=None):
    """
    Power-law generator for pdf(x)\propto x^{g-1} for a<=x<=b
    """
    if seed is not None:
        np.random.seed(seed)
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

def mag_err_to_rel_flux_err(dm):
    '''
    Convert uncertainty given in mag to
    uncertainty given in rel. flux units.
    Note: this is asymmetric!
    '''
    return np.power(100,dm/5.)-1.

def rel_flux_err_to_mag_err(rf):
    '''
    Convert uncertainty given in rel. flux units
    to uncertainty given in mag.
    Note: this is asymmetric!
    '''
    return 2.5*np.log10(rf+1.)


def rename_columns_for_consistency(df):
    df = df.rename(index=str, columns={'phot_g_mean_flux_Gaia' : 'Gflux_Gaia',
                                       'phot_g_mean_flux_error_Gaia' : 'e_Gflux_Gaia',
                                       'phot_g_mean_mag_Gaia' : 'G_Gaia',
                                       'phot_bp_mean_flux_Gaia' : 'BPflux_Gaia',
                                       'phot_bp_mean_flux_error_Gaia' : 'e_BPflux_Gaia',
                                       'phot_bp_mean_mag_Gaia' : 'BP_Gaia',
                                       'phot_rp_mean_flux_Gaia' : 'RPflux_Gaia',
                                       'phot_rp_mean_flux_error_Gaia' : 'e_RPflux_Gaia',
                                       'phot_rp_mean_mag_Gaia' : 'RP_Gaia',
                                       'bp_rp_Gaia' : 'BPRP_Gaia',
                                       'e_bp_min_rp_val_Gaia' : 'Ext_BPRP_Gaia',
                                       'e_bp_min_rp_percentile_lower' : 'elo_Ext_BPRP_Gaia',
                                       'e_bp_min_rp_percentile_upper' : 'eup_Ext_BPRP_Gaia',
                                       'Jmag_2MASS' : 'J_2MASS',
                                       'Hmag_2MASS' : 'H_2MASS',
                                       'Kmag_2MASS' : 'K_2MASS',
                                       'e_Jmag_2MASS' : 'e_J_2MASS',
                                       'e_Hmag_2MASS' : 'e_H_2MASS',
                                       'e_Kmag_2MASS' : 'e_K_2MASS',
                                       'gmag_PS1' : 'g_PS1',
                                       'e_gmag_PS1' : 'e_g_PS1',
                                       'rmag_PS1' : 'r_PS1',
                                       'e_rmag_PS1' : 'e_r_PS1',
                                       'imag_PS1' : 'i_PS1',
                                       'e_imag_PS1' : 'e_i_PS1',
                                       'zmag_PS1' : 'z_PS1',
                                       'e_zmag_PS1' : 'e_z_PS1',
                                       'ymag_PS1' : 'y_PS1',
                                       'e_ymag_PS1' : 'e_y_PS1',})
    todropcols = ['gKmag_PS1', 'e_gKmag_PS1', 'rKmag_PS1', 'e_rKmag_PS1',
                 'iKmag_PS1', 'e_iKmag_PS1','zKmag_PS1', 'e_zKmag_PS1',
                 'yKmag_PS1', 'e_yKmag_PS1']
    todropleft = list(set(list(df.columns.values)) & set(todropcols))
    df = df.drop(todropleft, axis=1)
    return df

def drop_if_no_photometry_source(df):
    '''
    Select only those that have entries with at least one of the following ids
    '''
    suff = ['2MASS_2MASS', 'source_id_Gaia', 'objID_PS1']
    condition = ((df[suff[0]].isnull()) & df[suff[1]].isnull() & df[suff[2]].isnull())
    try:
        df.insert(1, 'todrop', np.nan)
    except ValueError:
        logger.info('\"todrop\" columns exists.')
        pass
    df.loc[condition,'todrop'] = 'no photometry'
    logger.info('No photometry available for {} targets.\n'
                'Affected targets marked with \"no photometry\" in \"todrop\" column'.format(df[condition].shape[0]))
    return df

def remove_bad_2MASS_photometry(df):
    tempcols = ['tJ','tH','tK']
    df[tempcols] = df['Qfl_2MASS'].apply(lambda x: np.array(list(str(x)))=='A').apply(pd.Series)
    df.loc[df.tJ==False,'J_2MASS']= np.nan
    df.loc[df.tH==False,'H_2MASS']= np.nan
    df.loc[df.tK==False,'K_2MASS']= np.nan
    df = df.drop(tempcols,axis=1)
    return df

def remove_bad_PS1_photometry(df):

    def helper(x):
        if np.isnan(x):
            return "".join(['0']*25)
        else:
            return bin(int(x))[2:].zfill(25)

      # deprecated version

#     bands = list('grizy')
#     for b in bands:
#         f = df['{}Flags_PS1'.format(b)].apply(helper)
#         condition = ((f.str[-5]=='1') | (f.str[-6]=='1') | (f.str[-7]=='1'))# one flag must be set at least
#         df.loc[(~condition), ['{}_PS1'.format(b),'e_{}_PS1'.format(b)]] = np.nan # remove value if condition is not fulfilled

    # Sarah's version:

    bands = list('grizy')
    f = df['Qual_PS1'].apply(helper)
    print(f.str[-3].values)
    condition = (f.str[-3]=='1')# QF_OBJ_GOOD flag must be set at least
    print("Good PS1 photometry in {} case(s).".format(np.where(condition)[0].shape[0]))
    for b in bands:
        df.loc[(~condition), ['{}_PS1'.format(b),'e_{}_PS1'.format(b)]] = np.nan # remove value if condition is not fulfilled

    # another version

#     bands = list('grizy')
#     for b in bands:
#         f = df['{}Flags_PS1'.format(b)].apply(helper)
#         condition = ((f.str[-9]=='1') | (f.str[-10]=='1') | (f.str[-11]=='1') | (f.str[-12]=='1'))# one flag must be set at least
#         df.loc[(~condition), ['{}_PS1'.format(b),'e_{}_PS1'.format(b)]] = np.nan # remove value if condition is not fulfilled

    return df

def remove_bad_Gaia_photometry(df):
    '''
    Removes measurements in all three bands,
    where flux/uncertainty < 10, or where no
    flux or uncertainty is given.

    Parameters:
    ------------
    df : DataFrame
        Contains "Xflux_Gaia", "e_Xflux_Gaia", and "X_Gaia"
        columns.

    Return:
    --------------
    Input DataFrame with low S/N values removed.
    '''
    GaiaBands = ['G', 'BP', 'RP']
    for B in GaiaBands:
        df.loc[df['{}flux_Gaia'.format(B)].isnull(),'{}_Gaia'.format(B)] = np.nan
        df.loc[df['e_{}flux_Gaia'.format(B)].isnull(),'{}_Gaia'.format(B)] = np.nan
        #if there is no magnitude error calculated yet, skip but warn:
        if np.where(df.columns.values == 'e_{}_Gaia'.format(B))[0].shape[0] > 0:
            df.loc[df['{}flux_Gaia'.format(B)]/df['e_{}flux_Gaia'.format(B)] < 10.,['e_{}_Gaia'.format(B), '{}_Gaia'.format(B)]] = np.nan
        else:
            logger.info('e_{}_Gaia '.format(B), 'not in index. FYI.')
            df.loc[df['{}flux_Gaia'.format(B)]/df['e_{}flux_Gaia'.format(B)] < 10.,'e_{}_Gaia'.format(B)] = np.nan
    return df

def calculate_e_mag_for_Gaia(df):
    '''
    Gaia archive does not contain magnitude errors, because
    they would be somewhat asymmetric. We compute them anyway
    as a comparison to 2MASS and PS1, as these give magnitude
    errors nonentheless.
    '''
    GaiaBands = ['G', 'BP', 'RP']
    for B in GaiaBands:
        df[ 'e_{}_Gaia'.format(B)] = (df['e_{}flux_Gaia'.format(B)] / df['{}flux_Gaia'.format(B)]).apply(rel_flux_err_to_mag_err)
    return df

def correct_BP_RP_for_extinction(df):
    """Use internal E(BP-RP) values from Gaia
    and their uncertainties. Calculate a rough
    estimate for the uncertainty on the corrected
    value. If value or uncertainty on extinction are
    not available we just take the average value
    from the cluster as best estimate.
    """
    df.eup_Ext_BPRP_Gaia = df.eup_Ext_BPRP_Gaia.fillna(df.eup_Ext_BPRP_Gaia.mean())
    df.elo_Ext_BPRP_Gaia = df.elo_Ext_BPRP_Gaia.fillna(df.elo_Ext_BPRP_Gaia.mean())
    df.Ext_BPRP_Gaia = df.Ext_BPRP_Gaia.fillna(df.Ext_BPRP_Gaia.mean())
    df['BPRP_Gaia_corr'] = df.BPRP_Gaia - df.Ext_BPRP_Gaia
    df['e_BPRP_Gaia_corr'] = np.sqrt(df.e_BP_Gaia**2 +
                                     df.e_RP_Gaia**2 +
                                     ((df.eup_Ext_BPRP_Gaia +
                                       df.elo_Ext_BPRP_Gaia)/2)**2)
    return df


def convert_PS1_to_Sloan(df, **kwargs):
    """
    Apply PanSTARRS1 to SDSS grizy band magnitude
    conversion to all bands and keep the
    uncertainties.

    Parameters:
    ------------
    kwargs : dict
        Keyword arguments to pass to :func:`PS_to_SDSS`

    Return:
    --------
    DataFrame with extra columns for SDSS grizy band
    magnitudes and the uncertainties.
    """
    bands = list('grizy')
    for b in bands:
        try:
            sdss, sdsserr = PS_to_SDSS(df,b, **kwargs)
            df['{}_SDSS'.format(b)] = sdss
            df['e_{}_SDSS'.format(b)] = sdsserr
        except KeyError:
            message = ("\nPanSTARRS DR1 {0} band is missing in data."
                       "\nSkipped for conversion to SDSS {0}.".format(b))
            warn(message)
    return df


def PS_to_SDSS(params, band):
    '''
    Converts PS1 grizy photometry to SDSS grizy.
    Tested implicitly via :func:`convert_PS1_to_Sloan`

    Input:
    -----------
    params - dataframe: Target list with at least g and i magnitude
                        in PS1 given with gmag and imag as columns.
    band - str: One of the following: u, g, r, i, z, y, found in
                params with band+mag column.

    Return:
    ------------
    sdss - series: same length as params containing the converted SDSS
                   band magnitude
    '''
    x = params.g_PS1-params.i_PS1

    df = pd.read_csv('opencluster/static/SDSS_to_PS.csv')
    df = df[df.Band==band].iloc[0]
    band+='_PS1'

    # Convert PS1 to SDSS
    sdss = params[band] - df.a_0 - df.a_1 * x - df.a_2 * (x**2) - df.a_3 * (x**3) #sdss conversion

    # Propagate uncertainties:
    gierr = np.sqrt(params.e_g_PS1**2 + params.e_i_PS1**2) #uncertainty in g-i
    xerr_2 = ((df.a_1 - 2. * df.a_2 * x + 3. * df.a_3 * (x**2)) * gierr)**2 # variance in the part of sdss conversion term w/ g-i
    sdsserr = np.sqrt(params["e_"+band]**2 +  xerr_2) # quadratic propagation

    return sdss, sdsserr


def prepare_and_validate_stars(df, extinction=True):
    '''
    All preparatory steps to make table
    useful for stellar parameters determination.
    '''
    #check if all entries have K2 light curves
    #assert df[~df.RAJ2000_K2.isnull()].shape == df.shape #some have extra sources like SF and N for PSF de-trended light curves

    #drop target if no photometry is available
    logger.info('First check for lack of photometry: ')
    df = drop_if_no_photometry_source(df)

    #Rename columns so that one can easily iterate over multiple photometry sources
    df = rename_columns_for_consistency(df)

    # EPIC ids and Campaign should bei ints
    df[['EPIC','Campaign']] = df[['EPIC','Campaign']].fillna('0').astype(int)
    df.loc[df.EPIC==0, ['EPIC','Campaign']] = np.nan

    #calculate magnitude uncert for Gaia, ignore asymmetry
    df = calculate_e_mag_for_Gaia(df)

    # flux over flux_error for Gaia must be > 10
    df = remove_bad_Gaia_photometry(df)

    # Fill in the median cluster distance for targets without Gaia parallax:
    #df.loc[df['distance'].isnull(),"distance"] = df['distance'].median()

    # Use internal Gaia extinction correction to correct BP-RP
    if extinction == True:
        df = correct_BP_RP_for_extinction(df)
    else:
        df["BPRP_Gaia_corr"] = df["BPRP_Gaia"]
        df["e_BPRP_Gaia_corr"] = np.sqrt(df.e_BP_Gaia**2 + df.e_RP_Gaia**2)

    # measurements for each band to be A for 2MASS
    df = remove_bad_2MASS_photometry(df)

    # QF_OBJ_GOOD should be set in PS1 photometry
    df = remove_bad_PS1_photometry(df)

    #Convert PS1 photometry to SDSS:
    df = convert_PS1_to_Sloan(df)

    # Correct the transformed photometry for extinction using 3D dust maps
    if extinction == True:
        df = correct_for_extinction(df)

    #Drop target if no  g o o d  photometry is available
    logger.info('Second check for lack of photometry: ')
    df = drop_if_no_photometry_source(df)

    return df

def load_pickle(path):
    '''
    Load a pickled file, FFD, or OpenCluster object.
    '''
    #rb = read + binary
    obj = pickle.load(open(path, "rb" ) )
    logger.info('Load pickle: {}'.format(obj.__repr__()))
    return obj
