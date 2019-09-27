import numpy as np
import pandas as pd

from astropy import units as u
from pytest import approx

import copy
import pickle
import logging

from opencluster.utils import (prepare_and_validate_stars,
                               _calculate_observation_time)

from opencluster.cmd import (Teff_Boyajian,
                             Teff_Mann,
                             Teff_Apsis,
                             prioritize_Teff_Mann,
                             color_2MASS_to_Johnson,
                             cut_Teff_Mann_on_Teff,
                             cut_Teff_Mann_on_Mks,
                             correct_for_extinction,
                             calculate_weighted_mean)

from opencluster.lum import (read_Kepler_response,
                             Kepler_mean_wavelength,
                             interpolate_spectrum,
                             planck,
                             planck_full_blackbody,
                             projected_luminosity_Stefan_Boltzmann,
                             projected_luminosity_SED,
                             projected_luminosity_Kepler,
                             SME_find_library_match,
                             assign_SpT_to_Teff,
                             remove_outlier_Teffs,
                             )

from opencluster.radius import (radius_specmatchemp,
                                radius_mann,
                                calculate_double_check_radii_mann)

from opencluster.ffd import FFD

from opencluster import PACKAGEDIR

TEFFBINS = pd.read_csv('{}/clusters/Teff_bins_merged.csv'.format(PACKAGEDIR), dtype={'Teff_min':int,'Teff_max':int,})
UNITTOMODE = {'erg':'energy','s':'ED', None: None}

from .__init__ import logger

class OpenCluster(object):
    '''
    Implements a class for an open cluster.

    Attributes:'
    ------------
    cluster : str
        Data extension suitable cluster name.
    h_cluster : str
        Plot legend suitable cluster name.
    age : float
        Cluster age.
    u_age_high : float
        Upper cluster age uncertainty
    u_age_low : float
        Lower cluster age uncertainty
    age_unit : astropy.units.Quantity
        Age unit
    unit : 'erg' or 's'
        Either work in the ED or the flare energy
        domain
    feh : float
        [Fe/H] in dex
    u_feh : float
        [Fe/H] uncertainty in dex
    csym : str
        Plotting symbol for matplotlib.
    ccol : str
        Plotting color for matplotlib.
    flares : DataFrame
        Table of individual flares tagged by start/stop
        times and EPIC ID.
    stars : DataFrame
        Table of stellar properties.
    flares_on_stars : DataFrame
        Merged table of stars and flares.
    ffds : dict
        Dictionary of FFD object with Teff as keys.

    '''
    def __init__(self, cluster=None, h_cluster=None, age=None, u_age_high=0.,
                 u_age_low=0., age_unit=None, feh=None, unit=None, distance=None,
                 u_feh=None, csym='o', ccol='k', flares=None, stars=None,
                 flares_on_stars=None, ffds=dict(), prep=True, **kwargs):
        logger.info('Creating an OpenCluster object for {}'.format(cluster))
        self.distance = distance
        self.cluster = cluster
        self.h_cluster = h_cluster
        self.age = age
        self.u_age_high = u_age_high
        self.u_age_low = u_age_low
        self.age_unit = age_unit
        self.feh = feh
        self.u_feh = u_feh
        self.csym = csym
        self.ccol = ccol
        self.flares = flares
        self.stars = stars
        self.flares_on_stars = flares_on_stars
        self.ffds = ffds
        self.unit = unit

        if self.stars is not None:
            if (prep==True):
                self.stars = prepare_and_validate_stars(self.stars, **kwargs)
            self.stars['FeH'] = self.feh
            self.stars['e_FeH'] = self.u_feh
        logger.info('Created an OpenCluster object for {}'.format(cluster))

    def determine_Teff(self):
        '''
        Uses photometry information stored in ``stars`` to
        calculate effective temperature with Boyajian et al. (2013)
        Table 8 and Mann et al. (2015) Table 2 color-temperature
        relations. Adds a "Teff" column to ``stars`` if there's none.
        Keeps of colors that were used to determine Teff, and
        if Boyajian and/or Mann was used in the process.

        Additionally we include Gaia Apsis effective temperatures
        for FGK stars (Andrae et al. 2018).

        Return:
        --------
        Table of stellar parameters with additional tables with Teff
        info.
        '''
        df = copy.copy(self.stars)

        # Boyajian CTRs
        df = Teff_Boyajian(df)

        # Gaia Apsis Teffs
        df = Teff_Apsis(df)

        # Mann CTRs for low mass stars
        df = Teff_Mann(df)
        df = prioritize_Teff_Mann(df)
        df = cut_Teff_Mann_on_Mks(df)
        df = cut_Teff_Mann_on_Teff(df)
       # teff = df.loc[:,df.columns.str.startswith('Teff_')].columns.values
       # eteff = df.loc[:,df.columns.str.startswith('e_Teff_')].columns.values
        df = remove_outlier_Teffs(df)

        teff = df.loc[:,df.columns.str.startswith('Teff_')].columns.values
        std = df[teff].std(axis=1, numeric_only=True).fillna(0)
        df = calculate_weighted_mean(df)

        # Add std of individual values to it

        df.Teff_std = np.sqrt(df.Teff_std**2 + std**2)
        df.loc[df.Teff_std/df.Teff_median > .10, "Teff_median"] = np.nan
        df.loc[df.Teff_std/df.Teff_median > .10, teff] = np.nan

        self.stars = df
        print("Number of Teff_median: ",df.Teff_median.count())

        return df[['EPIC','Teff_median','Teff_std']]

    def find_Rstar(self, lib, mode="specmatchemp"):
        '''
        Use SpecMatch-Emp and/or Mann et al. (2015) Teff-R relations
        to match temperature and metallicity to stellar radius. Also
        gives uncertainty on radius and add all the relevant info from
        SpecMatch-Emp to ``stars`` if not yet done.

        Parameters:
        --------------
        lib : SpecMatch-Emp library

        Return:
        --------
        DataFrame with stellar radius and uncertainties
        '''

        df = copy.copy(self.stars)
        Tt= 3500
        cols = df.columns.values
        if (('Teff_median' in cols) & ('FeH' in cols)):

            if mode=="specmatchemp":
                df = radius_specmatchemp(df, lib, Teffthresh=Tt)
            elif mode=="mann":
                df = radius_mann(df)
            df = calculate_double_check_radii_mann(df)

            # Do the consistency check with MKs derived radii whenever there is a meaninggful
            # value R from MKs available (not NaN, Teff < 3800 K, uncertainty < 50%)
            # Do they overlap within uncertainties?
            df["Rstar_consistent"] = np.nan
            g = lambda x: x.Rstar_double_check == approx(x.Rstar, abs=x.e_Rstar + x.e_Rstar_double_check)
            df.loc[df.Teff_median < Tt, "Rstar_consistent"] = df[(df.Teff_median < Tt) &
                                                                    (df.e_Rstar_double_check * 2 <
                                                                     df.Rstar_double_check)].apply(g, axis=1)
            print("Number of Rstar: ", df.Rstar.count())
            self.stars = df
            return df
        else:
            raise AttributeError('Teff_median and/or FeH not assigned in to any star.\n'
                                 'Run determine_Teff first.')


    def find_SED(self, lib):
        '''
        Use SpecMatch-Emp to match temperature and metallicity to an
        empirical spectrum. Also gives uncertainty on spectrum and adds
        the relevant info from SpecMatch-Emp to ``stars`` if not yet
        done.

        Use Sarah J. Schmidt's spectral library for M and L dwarfs.

        Parameters:
        --------------
        lib : SpecMatch-Emp library

        Return:
        --------
        SED as in normalized flux (wavelength)
        '''

        def SMEspectrum(d, c):
            if np.isnan(d.Teff_median):
                return np.nan, np.nan, np.nan
            else:
                minind, res = SME_find_library_match(c, d.Teff_median, d.Teff_std, d.FeH)
                if np.isnan(minind):
                    return np.nan, np.nan, np.nan
                else:
                    spec = lib.get_spectrum(minind)
                    #set masked values (telluric lines) to 1. with no uncertainty
                    spec.s[np.isnan(spec.s)] = 1.
                    spec.serr[np.isnan(spec.serr)] = 0.
                    return spec.s, spec.serr, spec.w*u.angstrom

        def SJSspectrum(d, pathtosed="{}/static/ML_SEDs/".format(PACKAGEDIR)):
            """Only M0-L9 wihout non-integer SpTs. No uncertainties given."""
            try:
                sed = pd.read_csv("{}{}_SED.txt".format(pathtosed, d.SpT), delim_whitespace=True, names=["wav","val","e_val"])
                return sed.val.values, sed.e_val, sed.wav.values*u.angstrom
            except FileNotFoundError: #no spectrum found in library
                return np.nan, np.nan, np.nan



        df = copy.copy(self.stars)
        cols = df.columns.values
        if 'Teff_median' in cols:

            # use SpecMatch-Emp
            if 'FeH' in cols:
                cut = lib.library_params.query('logg > 4.2')# only select MS, luminosity class V
                df['res'] = df.apply(SMEspectrum, args=(cut,), axis=1)
                df[['SED_SME', 'e_SED_SME', 'w_SED_SME']] = df['res'].apply(pd.Series)
                df = df.drop('res', axis=1)
            # use Sarah J. Schmidt's spectral library for M and L dwarfs
            df = assign_SpT_to_Teff(df)
            df["res"] = df.apply(SJSspectrum, axis=1)
            df[['SED_SJS', 'e_SED_SJS', 'w_SED_SJS']] = df['res'].apply(pd.Series)
            df = df.drop("res", axis=1)

            cols = df.columns.values
            if (("SED_SJS" in cols) & ("SED_SME" in cols)):

                def compare(x):
                    if isinstance(x.SED_SJS, float):
                        x.SED_SJS = []
                    if isinstance(x.SED_SME, float):
                        x.SED_SME = []

                    return len(x.SED_SJS) > len(x.SED_SME)

                pref = df.apply(compare, axis=1)
                df["SED"] = np.nan
                df["e_SED"] = np.nan
                df["w_SED"] = np.nan
                df.loc[pref == True, "SED"] = df.SED_SJS
                df.loc[pref == False, "SED"] = df.SED_SME
                df.loc[pref == True, "e_SED"] = df.e_SED_SJS
                df.loc[pref == False, "e_SED"] = df.e_SED_SME
                df.loc[pref == True, "w_SED"] = df.w_SED_SJS
                df.loc[pref == False, "w_SED"] = df.w_SED_SME

            elif (("SED_SJS" in cols) & ("SED_SME" not in cols)):
                df["SED"] = df["SED_SJS"]
                df["e_SED"] = df["e_SED_SJS"]
                df["w_SED"] = df["w_SED_SJS"]

            self.stars = df
            print("number of SEDs: ",df.SED.count())
            return df[['EPIC','SED', 'e_SED', 'w_SED']]

        else:
            raise AttributeError('Teff_median and/or FeH not assigned to any star.\n'
                                 'Run determine_Teff first.')
            return df

    def calculate_luminosities(self):
        '''Script that calculates all luminosities needed
        '''


        C = copy.deepcopy(self.stars)
        C.Rstar = C.Rstar
        C.e_Rstar = C.e_Rstar
        C.Teff_median = C.Teff_median
        Kp = read_Kepler_response()
        Kepler_wav = Kepler_mean_wavelength(Kp)

        nonan = (C.w_SED.notnull() & C.SED.notnull()).values
        #interpolate spectrum to fit Kepler wavelength array
        spectra = C.loc[nonan,:].apply(lambda x: interpolate_spectrum(Kp.wav, x.w_SED, x.SED), axis=1)
        C['spectrum'] = np.nan
        C.loc[nonan,'spectrum'] = spectra

        #calculate planck curve at Kepler wavelength at given temperature
        Planck_curve = C.loc[nonan,:].apply(lambda x: planck(Kepler_wav, x.Teff_median), axis=1)
        e_Planck_curve = C.loc[nonan,:].apply(lambda x: planck(Kepler_wav, x.Teff_median, e_T=x.Teff_std, deriv=True), axis=1)

        C['planck_curve'] = np.nan
        C['e_planck_curve'] = np.nan
        C.loc[nonan,'planck_curve'] = Planck_curve
        C.loc[nonan,'e_planck_curve'] = e_Planck_curve


        #calculate projected luminosity from Planck curve
        Lum_SED = C.loc[nonan,:].apply(lambda x: projected_luminosity_SED(x.planck_curve, x.spectrum, Kepler_wav, x.Rstar), axis=1)
        e_Lum_SED = C.loc[nonan,:].apply(lambda x: projected_luminosity_SED(x.planck_curve, x.spectrum, Kepler_wav, x.Rstar,
                                                                          e_planck=x.e_planck_curve, e_R=x.e_Rstar, deriv=True), axis=1)
        C['Lum_SED'] = np.nan
        C.loc[nonan,'Lum_SED'] = Lum_SED
        C['e_Lum_SED'] = np.nan
        C.loc[nonan,'e_Lum_SED'] = e_Lum_SED

        #calculate projected luminosity from Planck curve and Kepler response function
        Lum_Kepler = C.loc[nonan,:].apply(lambda x: projected_luminosity_Kepler(x.planck_curve, x.spectrum, Kepler_wav, Kp.resp, x.Rstar), axis=1)
        e_Lum_Kepler = C.loc[nonan,:].apply(lambda x: projected_luminosity_Kepler(x.planck_curve, x.spectrum, Kepler_wav, Kp.resp, x.Rstar,
                                                                                e_planck=x.e_planck_curve, e_R=x.e_Rstar, deriv=True), axis=1)
        C['Lum_Kepler'] = np.nan
        C.loc[nonan,'Lum_Kepler'] = Lum_Kepler
        C['e_Lum_Kepler'] = np.nan
        C.loc[nonan,'e_Lum_Kepler'] = e_Lum_Kepler

        #calculate projected bolometric luminosity from stefan boltzmann law as a sanity check with the previous two luminosities
        Lum_SB = C.loc[nonan,:].apply(lambda x: projected_luminosity_Stefan_Boltzmann(x.Teff_median, x.Rstar), axis=1)
        e_Lum_SB = C.loc[nonan,:].apply(lambda x: projected_luminosity_Stefan_Boltzmann(x.Teff_median, x.Rstar,
                                                                                      deriv=True, e_T=x.Teff_std, e_R=x.e_Rstar), axis=1)

        C['Lum_SB'] = np.nan
        C.loc[nonan,'Lum_SB'] = Lum_SB
        C['e_Lum_SB'] = np.nan
        C.loc[nonan,'e_Lum_SB'] = e_Lum_SB

        self.stars = C
        return


    def validate_luminosities(self):
        '''
        Integration test to check if
        Lbol > LKep,
        also if
        Teff+Rstar+SED+Planck < Teff+Rstar+Stefan-Boltzmann.
        Maybe also if Lbol_phot using photometry + BC from Mann et al. (2015)
        is approximately Lbol from ``determine_Lbol``.
        '''
        return

    def validate_params(self, param='Teff'):
        '''
        Check your Teff, Lbol, or Rstar results with 'teff_val_Gaia',
        'radius_val_Gaia', or 'lum_val_Gaia'. Is 'lum_val_Gaia' really
        Lbol? << CHECK first:

        lum_val : stellar luminosity (float, Luminosity[Solar Luminosity])

        Estimate of luminosity from Apsis-FLAME
        '''
        return

    def merge_flares_stars(self):
        '''
        Match stellar parameters to detected flares.
        Return a flare table with luminosities, Teff,
        and other parameters supplemented.
        '''
        self.flares_on_stars = self.flares.merge(self.stars,how='left',on=['EPIC','Campaign'])
        return

    def generate_FFD(self, Tlims=None):
        '''
        Generate a FFD object from an entire cluster.
        Optionally cut on temperature.

        Parameters:
        ------------
        Tlims : tuple of floats
            (Tmin,Tmax)

        Return:
        ---------
        ffd : FFD object
        '''
        #some conversions
        df = copy.deepcopy(self.flares_on_stars)
        stars = copy.deepcopy(self.stars)
        df.loc[:,'Teff_median'] = df.Teff_median.apply(lambda x: u.Quantity(x).value)
        stars.loc[:,'Teff_median'] = stars.Teff_median.apply(lambda x: u.Quantity(x).value)
        if Tlims == None:
            Tlims = [df.Teff_median.min(), df.Teff_median.max()]
        # confine sample to Teff bins
        df = df[((df.Teff_median >= Tlims[0].value) & (df.Teff_median <= Tlims[1].value))]
        stars = stars[((stars.Teff_median >= Tlims[0].value) & (stars.Teff_median <= Tlims[1].value))]
        if ((stars.shape[0] == 0) | (df.shape[0] == 0)):
            #create empty FFD
            self.ffds[Tlims[0].value] = FFD(Tmax=Tlims[1], Tmin=Tlims[0],
                                            age=self.age, mode = UNITTOMODE[self.unit],)
            return
        else:
            df.Lum_Kepler = df.Lum_Kepler.apply(lambda x: u.Quantity(x).value)
            stars.Lum_Kepler = stars.Lum_Kepler.apply(lambda x: u.Quantity(x).value)
            #throw out targets that cannot be analysed in ED or energy space:
            if self.unit == 's':
                # at least a Teff must be known
                stars = stars.loc[~stars.Teff_median.isnull(),:]
                df = df.loc[~df.Teff_median.isnull(),:]
            elif self.unit == 'erg':
                # Kepler luminosity needed to find E_kp,flare
                stars = stars.loc[~stars.Lum_Kepler.isnull(),:]
                df = df.loc[~df.Lum_Kepler.isnull(),:]
            else:
                logger.error('No unit specified for OC {}'.format(self.h_cluster))
            #number of targets contains both flaring and non-flaring stars with Lum_Kep or Teff_median
            nstars = stars.groupby('EPIC').first().values.shape[0]
            ffd = FFD(name=self.cluster,
                      h_name=self.h_cluster,
                      ED=df.ed_rec.values,
                      energy=(df.ed_rec*df.Lum_Kepler).values,
                      Tmax=Tlims[1],
                      Tmin=Tlims[0],
                      flares=df,
                      age=self.age,
                      agemax=self.age + self.u_age_high,
                      agemin=self.age - self.u_age_low,
                      mode = UNITTOMODE[self.unit],
                      nstars=nstars,
                      tot_obs_time = _calculate_observation_time(stars, df))
            self.ffds[Tlims[0].value] = ffd
            return

    def dump_pickle(self, add='general'):
        '''
        Save OpenCluster object as a pickled file.

        Parameters:
        ------------
        add : 'general' or string
            specify why you saved this object

        Return:
        --------
        path to pickled file
         '''
        #wb = write + binary
        obj = 'OpenCluster'
        path = 'saved_objects/{}/{}_{}.p'.format(obj, self.name, add)
        pickle.dump(self, open(path, "wb" ) )
        return path

    def generate_all_FFDs(self, teffbins=TEFFBINS):
        '''
        Generate a dictionary of FFDs according
        to a table of Teff bins in the teffbins
        table to be found under OC.ffds.

        Parameters:
        -----------
        OC : OpenCluster

        teffbins : DataFrame
            table with Teff bin edges
        '''
        OC = copy.deepcopy(self)
        [OC.generate_FFD(Tlims=[row.Teff_min * u.Kelvin, (row.Teff_max - 1.) * u.Kelvin]) for i, row in teffbins.iterrows()]
        return OC
