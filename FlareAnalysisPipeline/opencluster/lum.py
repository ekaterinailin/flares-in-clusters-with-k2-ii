import pandas as pd
import numpy as np
from astropy import units as u
import copy
from astropy import constants as const
from astropy.constants import sigma_sb
from opencluster import PACKAGEDIR

def remove_outlier_Teffs(df):
    """Remove Teff results that use colors
    that are outliers in CMDs.
    """
    for id_ in range(df.shape[0]):
        df.todrop = df.todrop.astype(str)
        colors = [x.replace(" outlier", "") for x in df.todrop.iloc[id_].split(", ")]
        avail = df.columns[df.columns.str.contains("Teff")].values
        availidx = np.where(df.columns.str.contains("Teff"))[0]
        bands = np.array([x.strip("Teff_").strip("Mann_").split("_")[:2] for x in avail])
        newcolnames = [x for x in bands if len(x)==2]
        newcolsnameidx = [i for i, x in enumerate(bands) if len(x)==2]
        newavailidx = availidx[newcolsnameidx]
        assert len(newavailidx) == len(newcolnames)
        newcolnames = np.array(["{}-{}".format(x[0],x[1]) for x in newcolnames])
        a = [list(newavailidx[newcolnames == c]) for c in colors]
        a = [it for sl in a for it in sl]
        df.iloc[id_, a] = np.nan
    return df

def assign_SpT_to_Teff(df):
    '''Use Mamajek table to find a good SED match
    to what Sarah's library provides. If Teff falls
    right inbetween to values we take the hotter one.

    Parameters:
    -----------
    df : DataFrame
        must contain Teff_median

    Return:
    ---------
    '''
    def matchSpT(Teff, conv):
        spt = conv.loc[(conv.Teff - Teff).abs().argmin(),"SpT"]
        if spt == "K9":
            return ""
        else:
            return spt
    conv = pd.read_csv("{}/static/mamajek_teff_spt.csv".format(PACKAGEDIR))
    df.loc[~df.Teff_median.isnull(),"SpT"] = df[~df.Teff_median.isnull()].apply(lambda x: matchSpT(x.Teff_median, conv), axis=1)
    return df

def SME_find_library_match(cut, T, Terr, feh, interpolate=False):
    '''
    Search the best match to given Teff and metallicity
    in a cleaned library from SME.

    Parameters:
    ------------
    cut : specmatchemp.library
        cut-out from the entire library
    T : float
        effective temperature
    feh : float
        [Fe/H]

    Return:
    --------
    minind : int
        index into library for best match
    res : Series
        attributes of best matching target
    '''
    c = cut[(feh - 0.3 <= cut.feh) &
            (feh + 0.3 >= cut.feh)] #&
            #(cut.Teff >= T - Terr*3) &
          #  (cut.Teff <= T + Terr*3) ]

    #minimize distance in T-feh space
    if c.empty:
        #print('No matching entry in library for {}+-100 K and {}+-0.2.'.format(T,feh))
        return np.nan, np.nan
    else:
        if interpolate==False:

            minind = ((c.Teff - T)**2 + (c.feh - feh)**2).idxmin()
            res = c.loc[minind,:]
            return minind, res

        elif interpolate==True:

            c["delt"] = c.Teff-T
            if (c.delt.values > 0).all():
                l = c.sort_values(by="delt", ascending=True).iloc[0,:]
                h = c.sort_values(by="delt", ascending=True).iloc[1,:]
            elif (c.delt.values <= 0).all():
                h = c.sort_values(by="delt", ascending=False).iloc[0,:]
                l = c.sort_values(by="delt", ascending=False).iloc[1,:]
            else:
                l = c[c.delt < 0].sort_values(by="delt", ascending=False).iloc[0,:]
                h = c[c.delt >= 0].sort_values(by="delt", ascending=True).iloc[0,:]

            if h.Teff-l.Teff < Terr: # if the values are very close to each other take average value and take scatter into account
                cc = c[(c.Teff > T - Terr) & (c.Teff < T + Terr)]
                mr = cc.radius.mean()
                mrerr = np.sqrt(cc.radius.std()**2 + (cc.u_radius**2).sum()/cc.shape[0])

            else: #interpolate
                dx = (h.Teff - l.Teff)
                dy = (h.radius - l.radius)
                m = dy / dx
                c = (l.radius * h.Teff - h.radius * l.Teff) / dx

                mr = m * T + c
                errs = []
                errs.append(dy * (T - h.Teff) / dx**2 * l.u_Teff)
                errs.append(dy * (l.Teff - T) / dx**2 * h.u_Teff)
                errs.append((h.Teff - T) / dx * l.u_radius)
                errs.append((T - l.Teff) / dx * h.u_radius)
                errs.append(m * Terr)

                mrerr = np.sqrt(np.sum(np.array([errs])**2))

            return mr, mrerr

def read_Kepler_response():
    Kp = pd.read_csv('{}/static/Kepler_response.txt'.format(PACKAGEDIR),
                                      skiprows=9,
                                      header=None,
                                      delimiter='\t',
                                      names=['wav','resp'])
    Kp.wav = Kp.wav.astype(np.float)*10.*u.angstrom #convert to angström for spectrum function
    return Kp


def planck_full_blackbody(T):
    '''
    Calculate a full black body curve at a
    given T.

    Parameters:
    -----------
    T : float
        effective temperature in K
    Return:
    -------
    wav : array of floats
        wavelength array in cm
    planck : array of floats
        flux array in erg/(cm^3 * s)
    '''
    h = const.h.to('erg*s') # Planck constant in erg*s
    c = const.c.to('cm/s') # light speed in cm/s
    k = const.k_B.to('erg/K') # Boltzmann constant in erg/K
    wav = np.linspace(0.001,80000,10000)*1e-8*u.cm
    T = T * u.Kelvin
    planck = np.pi * 2.* h * (c**2) / (wav**5) / (np.exp( h * c / ( wav * k * T ) ) - 1. )
    assert planck.unit == u.Unit('erg * cm**-3 * s^-1')
    return wav, planck.value

def planck(w, T, e_T=None, deriv=False):
    '''
    Calculate the planck curve for a given array
    of wavelengths.

    Parameters:
    ------------
    w : numpy arrays of floats
        wavelength in cm
    T : float
        effective temperature in K

    Return:
    --------
    planck_curve : Planck curve
    '''
    planckunit = u.erg / u.s / (u.cm)**3
    h = const.h.to('erg*s')#*u.erg*u.s # Planck constant in erg*s
    c = const.c.to('cm/s')#*u.cm/u.s # light speed in cm/s
    k = const.k_B.to('erg/K')#*u.erg/u.Kelvin # Boltzmann constant in erg/K
    w = w.to('cm')
    T = T*u.Kelvin
    if deriv == False:
        val = np.pi* 2.* h * (c**2) / (w**5) / (np.exp( h * c / ( w * k * T ) ) - 1. )

    elif deriv == True:
        e_T = e_T*u.Kelvin
        e_ = np.exp( h * c / ( w * k * T ) )
        derivative = np.pi * 2. * h**2 * c**3 / w**6 / T**2 / k / (e_ - 1. )**2 * e_
        val = derivative * e_T

    assert val.unit == planckunit
    return val.value

def interpolate_spectrum(target_wav, spectrum_wav, spectrum,):
    '''

    Returns the one-dimensional piecewise linear
    interpolant to a function with given
    discrete data points (spectrum_wav, spectrum),
    evaluated at target_wav.

    Parameters:
    ------------
    target_wav : np.array
        wavelength array we want to work with
    spectrum_wav : np.array
        wavelength array we want to evaluate at the target
        wavelengths
    spectrum : np.array
        normalized flux values

    Return:
    ---------
        spectrum evaluated at target wavelengths
    '''

    return np.interp(target_wav, spectrum_wav, spectrum, left=1., right=1.)

def Kepler_mean_wavelength(Kp):
    '''
    Return with units.
    Determine the mean values between each two
    wavelengths given in Kp.wav.

    Parameters:
    ------------
    Kp : DataFrame
        Kepler response function with wavelength in angstrom
        but w/o astropy units.

    Return:
    -------
    Wavelength array with a 2-point window rolling mean for
    later avaluation.
    '''
    return Kp.rolling(2, min_periods=1).wav.mean().values*u.angstrom.to('cm')*u.cm

def projected_luminosity_Stefan_Boltzmann(T, R, e_T=None, e_R=None,
                                          deriv=False):
    '''
    Calculate projected luminosity from
    effective temperature and stellar radius.
    Takes values with astropy units only.

    Parameters:
    -----------
    T : float with astropy unit
        effective temperature
    R : float with astropy unit
        stellar radius

    Return:
    -------
    Stellar luminosity in erg/s.
    '''
    R = (R * u.solRad).to('cm')
    T = T * u.Kelvin
    assert R.unit == u.cm
    assert T.unit == u.Kelvin
    sigma_stefan = sigma_sb.to('erg * s**-1 * K**-4 * cm**-2')
    if deriv==False:
        val = (sigma_stefan * np.pi * R**2 * T**4).to('erg/s')
        return val.value
    elif deriv==True:
        e_R = (e_R * u.solRad).to('cm')
        e_T = e_T * u.Kelvin
        val = (2. * sigma_stefan * np.pi * R * T**3 *
               np.sqrt(T**2 * e_R**2 + 4. * R**2 * e_T**2)).to('erg/s')
        return val.value

def projected_luminosity_SED(planck, spectrum, wavelengths, R,
                             e_R=None, e_planck=None, deriv=False):
    '''
    Calculate the projected stellar luminosity, given the
    blackbody curve, a normalized spectrum and the stellar radius.

    Parameters:
    ------------
    planck : array with astropy units
        Blackbody SED
    spectrum : array, unitless
        Normalized stellar spectrum
    wavelengths : array with astropy units
        Wavelength array corresponding to planck and spectrum
    R : float with astropy units
        Stellar radius

    Return:
    -------
    Projected luminosity from the stellar SED in erg/s.
    '''
    R = (R * u.solRad).to('cm').value

    wavelengths = wavelengths.to('cm').value
    planck = planck
    flux_per_area = np.sum(np.diff(wavelengths) / 2. * (spectrum[1:] * planck[1:] + spectrum[:-1] * planck[:-1])) #trapezoidal rule
    if deriv==False:
        val = np.pi * R**2 * flux_per_area#).to('erg/s')
        return val
    elif deriv==True:
        e_R = (e_R * u.solRad).to('cm').value
        e_planck = e_planck
        e_fi = spectrum * e_planck
        e_flux_per_area = (.5 * np.sqrt(np.sum(np.diff(wavelengths)**2 * (e_fi[1:]**2 + e_fi[:-1]**2))))
        val = (np.pi * R * np.sqrt(4. * flux_per_area**2 * e_R**2 + R**2 * e_flux_per_area**2))
        return val


def projected_luminosity_Kepler(planck, spectrum, wavelengths, Kpresp, R,
                                e_planck=None, e_R=None, deriv=False):
    '''
    Calculate the projected stellar luminosity in the Kepler
    band, given the blackbody curve, a normalized spectrum,
    and the stellar radius.

    Parameters:
    ------------
    planck : array with astropy units
        Blackbody SED
    spectrum : array, unitless
        Normalized stellar spectrum
    wavelengths : array with astropy units
        Wavelength array corresponding to planck and spectrum
    Kpresp : array
        Kepler response function
    R : float with astropy units
        Stellar radius

    Return:
    -------
    Projected luminosity from the stellar SED in erg/s within the
    Kepler band.
    '''
    assert Kpresp.shape[0] == spectrum.shape[0]
    assert Kpresp.shape[0] == planck.shape[0]
    assert Kpresp.shape[0] == wavelengths.shape[0]
    R = (R * u.solRad).to('cm').value
    wavelengths = wavelengths.to('cm').value
    planck = planck
    flux_per_area = np.trapz(spectrum * planck * Kpresp, x=wavelengths)
    if deriv == False:
        val =  (np.pi * R**2 * flux_per_area) 
        return val
    elif deriv == True:
        e_R = (e_R * u.solRad ).to('cm').value
        e_fi = spectrum * Kpresp.values * (e_planck )
        e_flux_per_area = (.5 * np.sqrt(np.sum(np.diff(wavelengths)**2 * (e_fi[1:]**2 + e_fi[:-1]**2))))
        val = (np.pi * R * np.sqrt(4. * flux_per_area**2 * e_R**2 + R**2 * e_flux_per_area**2))
        return val
