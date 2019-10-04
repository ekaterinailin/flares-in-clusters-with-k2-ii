import pandas as pd
import numpy as np

from .__init__ import NORM, NORMG

def calculate_distance_from_parallax(df, check_GoF=False):
    '''Calculte distance in pc and error
    from Gaia parallax.

    Parameters:
    -----------
    df : DataFrame
        table from Gaia query, needs parallax, and
        parallax_error
    check_GoF : bool
        if True, will check unit weight error.
        df then also needs astrometric_chi2_al,
        astrometric_n_good_obs_al, phot_g_mean_mag,
        and optionally bp_rp.

    Return:
    -------
    DataFrame with a distance and a distance_error column.
    If UWE flag is set - more columns, see :func:check_UWE
    '''
    df['distance'] = 1. / df.parallax * 1e3
    df['distance_error'] = 1. / (df.parallax**2) * df.parallax_error * 1e3

    if check_GoF == True:
        df = check_UWE(df)
    return df

def UWE(df):
    '''Calculate the unit weight error.

    Parameters:
    ------------
    df : DataFrame
        table from Gaia query, needs astrometric_chi2_al,
         astrometric_n_good_obs_al, phot_g_mean_mag, and
         optionally bp_rp.

    Return:
    -------
    DataFrame with a UWE and UWE_norm column that
    indicates the unit weight error and the normalized UWE
    '''
    df['UWE'] = np.sqrt(df.astrometric_chi2_al/(df.astrometric_n_good_obs_al-5))
    df['UWE_norm'] = df.UWE/UWEnorm(df)
    return df

def UWEnorm(df):
    '''Calculate the normalization factor for the
    unit weight error.

    Parameters:
    ------------
    df : DataFrame
        table from Gaia query, needs phot_g_mean_mag, and
         optionally bp_rp.

    Return:
    --------
    Series of normalization factors u0
    '''
    args=[]
    for i, r in df.iterrows():
        if np.isnan(r.bp_rp):
            arg = ((r.phot_g_mean_mag - NORMG.g_mag)**2).idxmin()
            try:
                args.append(NORMG.u0.iloc[arg])
            except TypeError:
                args.append(np.nan)
        else:
            arg = ((r.phot_g_mean_mag - NORM.g_mag)**2 + (r.bp_rp - NORM.bp_rp)**2).idxmin()
            try:
                args.append(NORM.u0.iloc[arg])
            except TypeError:
                args.append(np.nan)

    df['u0'] = args
    return df.u0

def check_UWE(df):
    '''Check if the solution's goodness of fit
    unit normalized unit weight errors (UWE). Works better
    for fainter targets than for bright ones. See Gaia docs
    for more.

    Parameters:
    ------------
    df : DataFrame
        table from Gaia query, needs astrometric_chi2_al,
         astrometric_n_good_obs_al, phot_g_mean_mag, and
         optionally bp_rp.

    Return:
    ----------
    DataFrame with a boolean good_distance column that indicates
    the goodness of fit
    '''
    df = UWE(df)
    df['good_distance'] = UWE(df).UWE_norm < 1.4
    return df