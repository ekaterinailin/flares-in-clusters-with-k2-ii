import numpy as np
import pandas as pd

from scipy import interpolate

from opencluster.lum import SME_find_library_match
from gaia.gaia import calculate_distance_from_parallax

def radius_specmatchemp(df, lib, Teffthresh=3000):
    '''Define stellar radius as a best fit
    to empirical stars from the SpecMatch-Emp
    empirical library. Use Mann+15 Teff-R
    relation if Teff < Teffthresh.

    Parameters:
    -----------
    df : DataFrame
        stellar characteristics including
        Teff and [Fe/H]
    lib : library
        SpecMatch-Emp library
    Teffthresh : float
        Teff below which to use Mann et al. (2015) CTRs, too

    Return:
    -------
    df : DataFrame
        same as input but now with radius
    '''

    def radius(d, c):
            if np.isnan(d.Teff_median) | (d.Teff_median > 7000.):
                return np.nan, np.nan
            else:
                r, rerr = SME_find_library_match(c, d.Teff_median, d.Teff_std, d.FeH, interpolate=True)

                if np.isnan(r):
                    return np.nan, np.nan
                else:

                    return r, rerr


    cut = lib.library_params.query('(logg > 4.43 and Teff < 6000) or  (logg > 4.2 and Teff > 6000)')# only select MS, luminosity class V, > F4
    df['res'] = df.apply(radius, args=(cut,), axis=1)
    df[['Rstar', 'e_Rstar']] = df['res'].apply(pd.Series)
    df = df.drop('res',axis=1)

    #bc Specmatch-Emp does not cover <3000K, use Mann15 here:
    df.loc[df.Teff_median < Teffthresh,:] = radius_mann(df[df.Teff_median < Teffthresh])

    return df

def radius_mann(df):
    '''Apply Table 1 eqn 4 and 3 to determine
    R_* from Teff in Mann+2015 with and without
    [Fe/H], respectively.

    Parameters:
    -----------
    df : DataFrame
        Stellar parameters table

    Return:
    -------
    DataFrame with 'R_Mann'
    '''
    def eqn3(T, Terr=None, e=False):

        T = T / 3500.
        R = 10.5440 - 33.7546 * T + 35.1909 * T**2 -11.59280 * T**3
        if e == False:
            return R
        elif e == True:
            Terr = Terr / 3500.
            sig2_T = (- 33.7546 + 35.1909 * 2 * T -11.59280 * 3 * T**2)**2 * Terr**2
            sig2_R = (R * 0.134)**2
            return np.sqrt(sig2_T + sig2_R)

    def eqn4(T, feh, Terr=None, feherr=None, e=False):

        T = T / 3500.
        R = (16.7700 - 54.3210 * T + 57.6627 * T**2 -19.69940 * T**3) * (1. + 0.45650 * feh)
        if e == False:
            return R
        elif e == True:
            Terr = Terr / 3500.
            sig2_T = ((-54.3210 + 57.6627 * 2 * T -19.69940 * 3 * T**2) *
                      (1. + 0.45650 * feh) *
                      Terr)**2
            sig2_R = (R * .093)**2
            sig2_FeH = (((16.7700 - 54.3210 * T + 57.6627 * T**2 -19.69940 * T**3) * 0.45650) *
                        feherr)**2
            return np.sqrt(sig2_T + sig2_R + sig2_FeH)

    #conditions for applying each function:
    condition_eqn3 = ((df.FeH.isnull() | df.e_FeH.isnull()) & (~df.Teff_median.isnull())  & (df.Teff_median > 2700.) & (df.Teff_median < 4100.))
    condition_eqn4 = ((~df.FeH.isnull()) & (~df.e_FeH.isnull()) & (~df.Teff_median.isnull()) & (df.Teff_median > 2700.) & (df.Teff_median < 4100.))

    df['Rstar'] = np.nan
    df['e_Rstar'] = np.nan
    df.loc[condition_eqn3,'Rstar'] = df.loc[condition_eqn3,:].apply(lambda x: eqn3(x.Teff_median), axis=1)
    df.loc[condition_eqn3,'e_Rstar'] = df.loc[condition_eqn3,:].apply(lambda x: eqn3(x.Teff_median, Terr=x.Teff_std,
                                                                                     e=True), axis=1)
    df.loc[condition_eqn4,'Rstar'] = df.loc[condition_eqn4,:].apply(lambda x: eqn4(x.Teff_median, x.FeH), axis=1)
    df.loc[condition_eqn4,'e_Rstar'] = df.loc[condition_eqn4,:].apply(lambda x: eqn4(x.Teff_median, x.FeH,
                                                                                   Terr=x.Teff_std, feherr=x.e_FeH,
                                                                                   e=True), axis=1)
    return df


def calculate_double_check_radii_mann(df):
    """Calculate stellar radii using Gaia distances
    and 2MASS Ks band. Use these values as double
    check for Teff derived radii.

    Parameters:
    ------------
    df : DataFrame
        needs parallax, parallax_error, K_2MASS,
        and optionally e_K_2MASS, FeH, and e_FeH

    Return:
    -------
    DataFrame with MKs, MKs_err, Rstar_double_check,
             distance, distance_err, and
             e_Rstar_double_check extra columns.
    """
    def distance_modulus(d, derr):
        dm = 5. * np.log10(d/10.)
        dmerr = 5. / np.log(10.) * derr/d
        return dm, derr

    def eqn1(MKs, MKs_err=None, e=False):
        """Equation 1 from Table in Mann et al. 2016
        Erratum."""
        R = 1.9515 - 0.3520 * MKs + 0.01680 * MKs**2
        if e==False:
            return R
        elif e==True:
            sig2_R = (0.0289 * R)**2
            sig2_MKs = (0.3520 + 0.01680 * 2 * MKs)**2 * MKs_err**2
            return np.sqrt(sig2_R + sig2_MKs)

    def eqn2(MKs, FeH, MKs_err=None, FeH_err=None, e=False):
        """Equation 2 from Table in Mann et al. 2016
        Erratum."""
        R = (1.9305 - 0.3466 * MKs + 0.01647 * MKs**2) * (1. + 0.04458 * FeH)
        if e==False:
            return R
        elif e==True:
            sig2_R = (0.027 * R)**2
            sig2_MKs = ((0.3466 + 0.01647 * 2 * MKs) * (1. + 0.04458 * FeH) * MKs_err)**2
            sig2_FeH = ((1.9305 - 0.3466 * MKs + 0.01647 * MKs**2) * 0.04458 * FeH_err)**2
            return np.sqrt(sig2_R + sig2_MKs + sig2_FeH)

    # calculate the distance from Gaia parallax (no goodness of fit test):
    df["parallax"] = df.parallax_Gaia
    df["parallax_error"] = df.parallax_error_Gaia
    df = calculate_distance_from_parallax(df)

    # Calculate distance modulus and convert Ks band magnitude to MKs
    DM, DMerr = distance_modulus(df.distance, df.distance_error)
    df["MKs"] = df.K_2MASS - DM
    df["MKs_err"] = np.sqrt(df.e_K_2MASS**2 + DMerr**2)

    #conditions for applying equation 1 or 2:
    condition_eqn1 = ((df.FeH.isnull() | df.e_FeH.isnull()) & (~df.MKs.isnull())  & (df.MKs > 4.6) & (df.MKs < 9.8))
    condition_eqn2 = ((~df.FeH.isnull()) & (~df.e_FeH.isnull()) & (~df.MKs.isnull())  & (df.MKs > 4.6) & (df.MKs < 9.8))

    # Init new columns
    df['Rstar_double_check'] = np.nan
    df['e_Rstar_double_check'] = np.nan

    #Apply formula whenever conditions apply
    if (condition_eqn1 == True).any():
        df.loc[condition_eqn1,'Rstar_double_check'] = df.loc[condition_eqn1,:].apply(lambda x: eqn1(x.MKs), axis=1)
        df.loc[condition_eqn1,'e_Rstar_double_check'] = df.loc[condition_eqn1,:].apply(lambda x: eqn1(x.MKs, MKs_err=x.MKs_err,
                                                                                                      e=True), axis=1)
    if (condition_eqn2 == True).any():
        df.loc[condition_eqn2,'Rstar_double_check'] = df.loc[condition_eqn2,:].apply(lambda x: eqn2(x.MKs, x.FeH), axis=1)
        df.loc[condition_eqn2,'e_Rstar_double_check'] = df.loc[condition_eqn2,:].apply(lambda x: eqn2(x.MKs, x.FeH,
                                                                                   MKs_err=x.MKs_err, FeH_err=x.e_FeH,
                                                                                   e=True), axis=1)
    # delete extra columns
    del df["parallax"]
    del df["parallax_error"]
    del df["MKs"]
    del df["MKs_err"]
    del df["distance"]
    del df["distance_error"]


    return df