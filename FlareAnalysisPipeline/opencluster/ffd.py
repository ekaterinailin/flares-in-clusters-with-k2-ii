import pandas as pd
import numpy as np
import astropy.units as u
import pickle
import warnings
import logging
import operator

from copy import deepcopy
from scipy.optimize import fmin

from .linfit import linfit, define_linfit_arrays
from .powerlawfit import powerlawfit_Bauke2007 # for Bauke (2007)
from .powerlawfit import (ML_powerlaw_estimator,
                          de_biased_upper_limit,
                          de_bias_alpha,) #for Maschberger and Kroupa (2009)
from .statistics import (calculate_average_number_of_exceeding_values,
                         stabilised_KS_statistic,
                         calculate_KS_acceptance_limit,
                         calculate_cumulative_powerlaw_distribution)#for Maschberger and Kroupa (2009)
from .utils import (set_plotlimits,
                    _calculate_observation_time,
                    timestamp,
                    generate_random_power_law_distribution)

from .__init__ import logger

import matplotlib.pyplot as plt



class FFD(object):
    '''
    Attributes:
    -----------
    nstars : int
        Number of unique targets in Teff bin,
        including non-flaring targets
    nflares : int
        Number of detected flares in sample
    ED : array
        EDs in sample (s)
    ED_err : array
        ED uncertainties in sample (s)
    energy: array
        energies in samples (erg)
    energy_err: array
        energy error in samples (erg)
    tot_obs_time: astropy Quantity
        total observation time of all targets
        in Teff bin in years, including those
        where no flares were detected.
    Tmin : astropy Quantity
        minimum Teff of stars in sample
    Tmax : astropy Quantity
        maximum Teff of stars in sample
    alpha : float
        powerlaw exponent
    alpha_err : float
        powerlaw exponent uncertainty
    alpha_mode: str
        How was alpha determined?
    beta : float
        power law intercept
    beta_err : float
        power law intercept uncertainty
    mode : str
        ED or energy
    cutoff : astropy Quantity
        cutoff that was set last
    cutoff_energy_lower: astropy Quantity
        lower energy cutoff for the FFD in erg
    cutoff_energy_upper: astropy Quantity
        upper energy cutoff for the FFD in erg
    cutoff_ED_lower: astropy Quantity
        lower ED cutoff for the FFD in s
    cutoff_ED_upper: astropy Quantity
        upper ED cutoff for the FFD in s
    FR : float
        Flaring rate in 1/yr above cutoff
    FR_err : float
        Poissonian FR uncertainty
    flarelum : float
        flaring luminosity in erg/s, using Kepler flare energy
    KS : bool
        Is the fitted model consistens with a power law?
    truncated : bool
        Is the fitted powerlaw truncated, judging by an exceedance test?
    '''
    def __init__(self, ED=None, energy=None, tot_obs_time=None,
                 energy_err=None, ED_err=None, name=None,
                 h_name=None, Tmax=None, Tmin=None, flares=None,
                 alpha=None, beta=None, alpha_err=None, beta_err=None,
                  age=None, agemax=None, agemin=None, nstars=None,
                 multiplicity=None, alpha_mode=None,
                 truncated=None, KS=None, mode=None, cutoff=None,
                 cutoff_energy_lower=None, cutoff_ED_lower=None,
                 cutoff_energy_upper=None, cutoff_ED_upper=None,
                 FA=None, FR=None, FR_err=None, flarelum=None):
        #if data are give in the flares table: infer ED and energy
        if (flares is not None):
            if ('ed_rec' in flares.columns):
                self._ED = flares.ed_rec.values
                self.nflares = len(flares.ed_rec.values)
                logger.info('Flare sample taken from flares attribute.')
        else:
            self._ED = ED
            if ED is not None:
                self.nflares = len(ED)
        #if data are NOT given in a flares table: init a flares table from ED
        #do not trust energies that are given without context
        if ((flares is None) & (ED is not None)):
            self.flares = pd.DataFrame({'ed_rec' : ED})
        elif ((flares is None) & (energy is not None)):
            warnings.warn('Where did you find these energies?')
        else:
            self.flares = flares

        self._energy = energy

        self.tot_obs_time = tot_obs_time
        self.energy_err = energy_err
        self.ED_err = ED_err
        self.name = name
        self.h_name = h_name
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nstars = nstars
        self.alpha = alpha
        self.beta = beta
        self.alpha_err = alpha_err
        self.beta_err = beta_err
        self.age = age
        self.agemax = agemax
        self.agemin = agemin
        self.multiplicity = multiplicity
        self.alpha_mode = alpha_mode
        self.truncated = truncated
        self.KS = KS
        self._cutoff = cutoff
        self._cutoff_energy_lower = cutoff_energy_lower
        self._cutoff_ED_lower = cutoff_ED_lower
        self._cutoff_energy_upper = cutoff_energy_upper
        self._cutoff_ED_upper = cutoff_ED_upper
        self._mode = mode
        self.sol_alpha_mean = None
        self.sol_alpha_err_mean = None
        self.FA = FR
        self.FR = FR_err
        self.FR = FA
        self.flarelum = flarelum

        #Logging output modifications:
        if self.ED is not None:
            _nflares = len(self.ED)
        else:
            _nflares = 'no'
        if self.flares is not None:
            if 'EPIC' in self.flares:
                _flarestars = len(list(set(self.flares.EPIC)))
            else:
                _flarestars = 'any'
        else:
            _flarestars = 'any'
        if self.nstars is not None:
            _nstars = self.nstars
        else:
            _nstars = 'No'

        #Actual logging output:
        logger.info('Creating an instance of FFD for {}'
                    ' within {}-{}. {} stars with {}'
                    ' flares on {} of them.'.format(self.name, self.Tmin, self.Tmax,
                                                    _nstars, _nflares, _flarestars))
        #clean
        self.clean_FFD()

    def clean_FFD(self):
        '''
        Clean sample from invalid data and write to log about it.
        '''
        if self.flares is not None:
            if self.flares.ed_rec.values.shape[0] < 5:
                logger.info('Less than 5 flare candidates for {}-{}'.format(self.Tmin, self.Tmax))
            if (self.flares.ed_rec.values <= 0.).any():
                locn = self.flares[self.flares.ed_rec <= 0.]
                self.flares = self.flares[self.flares.ed_rec >0.]
                logger.warning('{} invalid values encountered in ED, and dropped.'.format(locn.shape[0]))
                logger.info('Dropped values: {}'.format(locn.ed_rec))
                self.ED = self.flares.ed_rec.values #automatically updates nflares too
            if (self.flares.ed_rec.isnull()).any():
                locn = self.flares[self.flares.ed_rec.isnull()]
                self.flares = self.flares[~self.flares.ed_rec.isnull()]
                logger.warning('{} invalid values encountered in ED, and dropped.'.format(locn.shape[0]))
                logger.info('Dropped values: {}'.format(locn.ed_rec))
                self.ED = self.flares.ed_rec.values #automatically updates nflares too
        return

    @property
    def ED(self):
        return self._ED

    @ED.setter
    def ED(self, ED):
        self.nflares = len(ED)
        self._ED = ED

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, energy):
        self.nflares = len(energy)
        self._energy = energy

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, m):
        self._mode = m
        message = 'A different mode was used for this FFD.'
        if m == 'ED':
            if self.cutoff_ED_lower is None:
                self.cutoff = 0. * u.s
            assert self.cutoff_energy_lower is None, message
            assert self.cutoff_energy_upper is None, message

        if m == 'energy':
            if self.cutoff_energy_lower is None:
                self.cutoff = 0. * u.erg
            assert self.cutoff_ED_lower is None, message
            assert self.cutoff_ED_upper is None, message


    #INFER CUTOFF
    @property
    def cutoff(self):
        return self._cutoff

    @cutoff.setter
    def cutoff(self, cutoff):
        self._cutoff = cutoff
        if (cutoff.unit == u.erg):
            self.cutoff_energy_lower = cutoff
        elif (cutoff.unit == u.s):
            self.cutoff_ED_lower = cutoff


    #LOWER CUTOFF
    @property
    def cutoff_ED_lower(self):
        return self._cutoff_ED_lower

    @cutoff_ED_lower.setter
    def cutoff_ED_lower(self, cutoff, mode='ED'):
        self._get_threshold_warning(cutoff, mode=mode, edge='lower')
        self._cutoff_ED_lower = u.Quantity(cutoff, 's')
        self._cut_threshold(mode, edge='lower')
        self.mode = mode

    @property
    def cutoff_energy_lower(self):
        return self._cutoff_energy_lower

    @cutoff_energy_lower.setter
    def cutoff_energy_lower(self, cutoff, mode='energy'):
        self._get_threshold_warning(cutoff, mode=mode, edge='lower')
        self._cutoff_energy_lower = u.Quantity(cutoff, 'erg')
        self._cut_threshold('energy', edge='lower')
        self.mode = 'energy'

    #UPPER CUTOFF
    @property
    def cutoff_ED_upper(self):
        return self._cutoff_ED_upper

    @cutoff_ED_upper.setter
    def cutoff_ED_upper(self, cutoff, mode='ED'):
        self._get_threshold_warning(cutoff, mode=mode, edge='upper')
        self._cutoff_ED_upper = u.Quantity(cutoff, 's')
        self._cut_threshold(mode, edge='upper')
        self.mode == mode

    @property
    def cutoff_energy_upper(self):
        return self._cutoff_energy_upper

    @cutoff_energy_upper.setter
    def cutoff_energy_upper(self, cutoff, mode='energy'):
        self._get_threshold_warning(cutoff, mode=mode, edge='upper')
        self._cutoff_energy_upper = u.Quantity(cutoff, 'erg')
        self._cut_threshold(mode, edge='upper')
        self.mode = mode

    def _get_threshold_warning(self, cutoff, mode='ED', edge='upper'):
        '''Trigger a warning if you set an energy
        or ED threshold twice'''
        assert isinstance(cutoff, u.Quantity), 'You need an astropy quantity.'
        message = 'Proceeding on an incomplete dataset.'
        if ((self.cutoff_ED_lower is not None)  & (edge=='lower')):
            logger.warning('Lower ED cutoff has been set before to {}'.format(self._cutoff_ED_lower))
            if cutoff < self.cutoff_ED_lower:
                logger.warning(message)
        if ((self.cutoff_energy_lower is not None) & (edge=='lower')):
            logger.warning('Lower ED or energy cutoff has been set before to {}'.format(self._cutoff_energy_lower))
            if cutoff < self.cutoff_energy_lower:
                logger.warning(message)
        if ((self.cutoff_ED_upper is not None)  & (edge=='upper')):
            logger.warning('upper ED cutoff has been set before to {}'.format(self._cutoff_ED_upper))
            if cutoff > self.cutoff_ED_upper:
                logger.warning(message)
        if ((self.cutoff_energy_upper is not None) & (edge=='upper')):
            logger.warning('upper ED or energy cutoff has been set before to {}'.format(self._cutoff_energy_upper))
            if cutoff > self.cutoff_energy_upper:
                logger.warning(message)

    def _cut_threshold(self, mode, edge='lower'):
        '''Internal function that cuts ffd.flares at given
        cutoff value in diff. modes.

        Only flares with energy/ED above cutoff_value are kept.

        Parameters:
        -----------
        mode : 'ED' or 'energy'
            different modes depending on the output FFD

        Return:
        --------
        FFD object with `flares` attribute without the
        flares at or below the energy threshold.
        '''
        if edge == 'lower':
            op = operator.gt
            threshED = self.cutoff_ED_lower
            threshenerg = self.cutoff_energy_lower
        elif edge == 'upper':
            op = operator.lt
            threshED = self.cutoff_ED_upper
            threshenerg = self.cutoff_energy_upper

        if np.isnan(self.flares.ed_rec.values).any():
            logger.warning('NaN encountered in ed_rec while cutting '
                           'ED thresholds. These are dropped automatically. '
                           'May cause problems later.')
        if mode == 'ED':
            df = self.flares.loc[op(self.flares.ed_rec, threshED.value), :]#ed
            df['en'] = df.ed_rec
            self.ED = df.ed_rec.values #automatically updates nflares too
        elif mode == 'energy':
            if self.flares.Lum_Kepler.values.shape[0]==0:
                raise ValueError('No luminosity data available.')
            if np.isnan(self.flares['Lum_Kepler']).any():
                logger.warning('NaN encountered in Lum_Kep while cutting '
                               'energy thresholds. These are dropped automatically. '
                               'May cause problems later.')
            df = self.flares.loc[op(self.flares.ed_rec * self.flares.Lum_Kepler, threshenerg.value), :]#energ
            df['en']  = df.ed_rec * df.Lum_Kepler
            self.energy = (df.ed_rec * df.Lum_Kepler).values #automatically updates nflares too

        self.flares = df

    def plot_FFD(self, ax, setplotlims=True, xlim=None, ylim=None, **kwargs):
        '''
        Plot and save a FFD of EDs.

        Parameters:
        -------------
        ax : Axes object
            matplotlib Axes to plot to
        setplotlims : True or bool
            set limits to the axes by default using the sample min/max,
            or deactivate this option with False
        xlim : None or tuple of floats
            limits to the plot's x-axis
        ylim : None or tuple of floats
            limits to the plot's y-axis
        kwargs : dict
            Keyowrd arguments to pass to
            :func:plt.scatter
        '''
        df = deepcopy(self.flares)

        if self.mode=='ED':
            energy = df.ed_rec
            ax.set_xlabel(r'$ED$ (s)', fontsize=12)
            ax.set_title('ED FFD {}'.format(self.h_name), fontsize=12)
        elif self.mode == 'energy':
            energy = df.ed_rec*df.Lum_Kepler
            ax.set_xlabel(r'$E_\mathrm{Kp}$' + ' (erg)', fontsize=12)
            ax.set_title('FFD {}'.format(self.h_name), fontsize=12)
        if len(energy) == 0:
            warnings.warn('Empty FFD. No plots.')
            return
        a, freq = define_linfit_arrays(df, mode=self.mode)
        collection = ax.scatter(a, freq/self.tot_obs_time, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        if setplotlims==True:
            set_plotlimits(ax, a, freq/self.tot_obs_time, xlim=xlim, ylim=ylim)
        ax.set_ylabel('Flare Frequency (yr' + r'$^{-1}$)', fontsize=12)
        return collection

    def plot_powerlaw(self, ax, custom_xlim=None, **kwargs):
        '''
        Plot the power law fit to the FFD.

        Parameters:
        -----------
        ax : matplotlibe Axes object
            plot to insert the power law in to
        kwargs : dict
            Keyword arguments to pass to plt.plot()

        Return:
        --------
        3 power law points to construct a line
        in log-log representation.
        '''
        if custom_xlim is None:
            if self.mode=='ED':
                data = self.ED
            elif self.mode=='energy':
                data = self.energy
            x = np.linspace(np.nanmin(data),np.nanmax(data),3)
        else:
            mi, ma = custom_xlim
            x = np.linspace(mi, ma, 3)
        y = self.beta / np.abs(self.alpha - 1.) * np.power(x,-self.alpha+1.)
        ax.plot(x, y,  **kwargs)
        return


    def plot_percentile_percentile(self, ax, sig_level=0.05, **kwargs):
        '''
        Plot the percentile-percentile, or
        probability-probability distribution, as
        suggested by Maschberger and Kroupa 2009.

        Parameters:
        --------------
        ax : Axes object
            panel to plot to
        sig_level : 0 < float < 1
            significance level for acceptance region
        '''
        if self.cutoff_ED_lower is not None:
            data = self.ED
        if self.cutoff_energy_lower is not None:
            data = self.energy
        alpha = self.alpha
        if alpha is None:
            raise ValueError('Compute power law exponent first.')
        sorted_data = np.sort(data)
        pp = calculate_cumulative_powerlaw_distribution(sorted_data, alpha)
        y = (np.arange(1, len(pp) + 1) - .5) / len(pp)
        limit = calculate_KS_acceptance_limit(len(data), sig_level=sig_level)
        ax.plot(pp, y, **kwargs)
        ax.plot(pp, pp + limit, c='k', label='$p = {}$'.format(sig_level))
        ax.plot(pp, pp - limit, c='k')
        ax.set_xlabel(r'$P_i$')
        ax.set_ylabel(r'S')
        ax.set_title('Percentile-percentile plot with acceptance region')
        return

    def fit_beta_to_powerlaw(self):
        '''
        Fit beta via non-linear least squares to a power
        law with given alpha using the cumulative
        FFD. Generate uncertainty using jackknife algorithm.
        '''
        def LSQ(x0,a,freq,alpha):
            zw = ((x0 / (np.power(a,alpha-1.) * (alpha-1.))-freq)**2).sum()
            return np.sqrt(zw)
        ffd = deepcopy(self)
        df = ffd.flares
        a, freq = define_linfit_arrays(df, mode=ffd.mode)
        N = len(a)
        if N==0:
            raise ValueError('No data.')
        #jackknife uncertainty
        x0starts = {'ED' : 100, 'energy' : 1e25}
        beta = np.array([fmin(LSQ,x0=x0starts[ffd.mode],
                              args=(np.delete(a,i),np.delete(freq,i),ffd.alpha),
                              disp=0)[0] for i in range(N)])
        #cumulative beta = beta_cum
        ffd.beta = beta.mean()/ffd.tot_obs_time
        ffd.beta_err = np.sqrt( (N-1) / N * ( (beta/ffd.tot_obs_time - ffd.beta)**2 ).sum() )
        #power law beta = beta_cum * |alpha-1|
        ffd.beta = ffd.beta * np.abs(ffd.alpha - 1.)
        #propagate errors on alpha to beta
        ffd.beta_err = np.sqrt(ffd.beta_err**2 * (ffd.alpha - 1.)**2 + ffd.beta**2 * ffd.alpha_err**2)
        return ffd

    def powerlaw_Maschberger_and_Kroupa(self, alims=[1.01,3.]):
        '''
        Calculate the un-biased ML power law estimator
        from Maschberger and Kroupa (2009), sections
        3.1.4. and 3.1.5.

        Parameters:

        ffd: FFD
            FFD for which to fit the power law
        alims:
            parameter range for power law exponent
        '''
        ffd = deepcopy(self)
        if 'en' not in ffd.flares.columns:
            ffd.cutoff_ED_lower = 0.*u.s
        data = ffd.flares.en.values
        a = np.linspace(alims[0], alims[1], 2000)
        N = len(data)
        alphaids = [np.argmin(np.abs(ML_powerlaw_estimator(np.delete(data,i), a))) for i in range(N)]
        alpha_ = a[alphaids]
        mean_alpha = np.mean(alpha_)
        sig_alpha = np.sqrt( (N-1) / N * ( (alpha_ - mean_alpha)**2 ).sum() )
        ffd.alpha = mean_alpha
        ffd.alpha_err = sig_alpha
        ffd.alpha_mode = 'MK'
        return ffd


    def linefit(self):
        '''
        Fit a line to a distribution of points.
        Calculate uncertainties on line parameters
        using the jacknife algorithm.

        TENTATIVE:
        Uncertainties in y are generated as poissonian.

        Parameters:
        ------------

        Return:
        --------
        m : tuple
            linear slope and its uncertainty
        b : tuple
            intercept and its uncertainty
        '''
        ffd = deepcopy(self)
        mode = self.mode
        if 'en' not in ffd.flares.columns:
            ffd.cutoff_ED = 0.*u.s
        df = ffd.flares
        if len(df.en.values)==0:
            raise ValueError('No data.')
        energy, freq = define_linfit_arrays(df, mode)
        df2 = pd.DataFrame({'x':np.log10(energy),
                            'y':np.log10(freq),
                            'sigy':np.log10(np.sqrt(freq)),
                            'sigx':0,
                            'rho': 0,
                            'norm':np.log10(ffd.tot_obs_time)})
        df2 = df2.iloc[:-1]
        #plt.scatter(df2.x,np.log10(freq[cutoff:-1]),c='k')
        m, b = linfit(df2)
        ffd.alpha_mode = 'linear fit'
        ffd.alpha = -m[0]+1
        ffd.beta = np.power(10,b[0])/ffd.tot_obs_time
        ffd.alpha_err = m[1]
        #need to be adapted to np.power(10,b[0]-/+b[1])
        ffd.beta_err = b[1]/ffd.tot_obs_time
        return ffd

    def powerlawfit(self):
        '''
        Determine the power law exponent using Bauke 2007
        solution. But use jacknife to determine uncertainty
        on alpha.

        May check later if uncertainty on alpha as provided
        by Bauke 2007 is worth considering (how does it compare in
        size to jacknife?) Random sampling from a power law gives
        broad variance in alpha, jackknife underestimates
        uncertainties.

        Return:
        ---------
        alpha and its uncertainty
        '''
        ffd = deepcopy(self)
        mode = self.mode
        #data = np.array([int(np.rint(i)) for i in data.en.values])
        if 'en' not in ffd.flares.columns:
            ffd.cutoff_ED_lower = 0.*u.s
        data = ffd.flares.en.values
        if ((data < 0) | np.isnan(data)).any():
            raise ValueError('Caught NaN or negative value in ed_rec data.\n'
                             'Check your inputs.')
        if len(data)==0:
            raise ValueError('No data.')
        else:
            alpha_ = [powerlawfit_Bauke2007(np.delete(data,i)) for i in range(len(data))]
            mean_alpha = np.mean(alpha_)
            N = len(data)
            sig_alpha = np.sqrt( (N-1) / N * ( (alpha_ - mean_alpha)**2 ).sum() )
            ffd.alpha_mode = 'Bauke'
            ffd.alpha = mean_alpha
            ffd.alpha_err = sig_alpha
            return ffd

    def calculate_R1s(self, sqr=False):
        '''
        Calculate the R1s measure - the cumulative
        flare frequency at ED=1s, extrapolated from
        the fitted power law.

        Parameters:
        ------------
        sqr : False or bool
            Calculate uncertainty on R1s either
            in quadrature (True) or as linear
            estimate (False).

        Return:
        -------
        R1s, sig_R1s : astropy Quantity
            R1s and its uncertainty
        '''
        ffd = deepcopy(self)
        R1s = ffd.beta / (ffd.alpha - 1.) #or equivalently, this is the cumulative beta value

        if ffd.alpha_err is None:
            logger.warning('No error for alpha given. Fill with 0.')
            ffd.alpha_err = 0.
            ae = 0.
        else:
            ae = ffd.alpha_err
        if ffd.beta_err is None:
            logger.warning('No error for beta given. Fill with 0.')
            ffd.beta_err = 0.
            be = 0.
        else:
            be = ffd.beta_err
        if sqr == False:
            sig_R1s = be / (ffd.alpha - 1.) + ae * ffd.beta / (ffd.alpha - 1.)**2
        elif sqr == True:
            sig_R1s = np.sqrt((be / (ffd.alpha - 1.))**2 + (ae * ffd.beta / (ffd.alpha - 1.)**2)**2)
        else:
            raise NameError('The sqr variable must be either True or False.')
        ffd.R1s = R1s
        ffd.R1s_err = sig_R1s
        return ffd

    def is_powerlaw_truncated(self, rejection=(.15, .05), nthresh=100):
        '''
        Apply the exceedance test recommended by
        Maschberger and Kroupa 2009. Update the
        `truncated` attribute in the FFD with either
        True or False.

        Parameters:
        ------------
        rejection : tuple of floats < 1.
            above these thresholds the distribution
            can be suspected to be truncated
        nthresh : int
            Number at which to use the more permissive
            or more restrictive truncation rejection
            limit, i.e. value 0 or 1 in `rejection`

        Return:
        ---------
        updated FFD object with `truncated` attribute
        '''
        ffd = deepcopy(self)
        mean, std = calculate_average_number_of_exceeding_values(ffd, 500)
        if ffd.alpha > 2.:
            warnings.warn('Power law exponent is steep. '
                          'Power of statistical tests decreases '
                          'according to Maschberger and Kroupa 2009.')
        if len(ffd.ED) >= nthresh:
            truncation_limit = rejection[1]
        else:
            truncation_limit = rejection[0]
        ffd.truncated = (mean/len(ffd.ED) > truncation_limit)
        if ffd.truncated == True:
            warnings.warn('This power law may be truncated.')
        return ffd

    def is_powerlaw_by_KS(self, sig_level=0.05):
        '''
        Test if we must reject the power law hypothesis
        judging by the stabilised Kolmogorov-Smirnov
        statistic, suggested by Maschberger and Kroupa
        2009.

        Parameters:
        -----------
        sig_level : float < 1.
            significance level for the hypothesis test

        Returns:
        ---------
        True if we cannot not reject the power law hypothesis.
        False if we must reject the power law hypothesis.
        '''
        ffd = deepcopy(self)
        if self.mode=='ED':
            data = ffd.ED
        elif self.mode =='energy':
            data = ffd.energy
        else:
            raise AttributeError('mode attribute not set to either energy or ED')
        ffd = ffd.is_powerlaw_truncated()
        KS = stabilised_KS_statistic(data, alpha=ffd.alpha, truncated=ffd.truncated)
        limit = calculate_KS_acceptance_limit(len(data), sig_level=sig_level)
        ffd.KS = KS < limit
        if ~ffd.KS:
            logger.warning('Kolmogorov-Smirnov tells us to reject'
                           r' the power law hypothesis at p={}.'
                           ' KS={}, limit={}'.format(sig_level, KS, limit))
        return ffd

    def calculate_FA(self):
        '''
        Calculate FA following Eqn. (11) in Ilin+2019.
        Flare luminosity of stars where no flares were found
        is set to zero.

        Return:
        -------
        FA : float
        '''
        def obstime(df):
            df = df.drop_duplicates(subset=['Campaign'])#in case multiple flares were observed in a campaign
            t = df.total_obs_time.sum() # sum over all campaigns
            return t/24./365.25 #convert to yr
        ffd = deepcopy(self)
        if 'en' not in ffd.flares.columns:
            ffd.cutoff_energy_lower = 0.*u.erg
            logger.warning('No mode or energy cutoff set. Default to 0.0 erg.')
        if ffd.flares.en.values.shape[0]==0:
            ffd.FA = np.nan
            logger.warning('Could not calculate FA. No values with energy available.')
            return ffd
        else:
            df = ffd.flares
            #total energy released in flares by each target
            fl = df.groupby('EPIC').en.sum()
            #bolometric luminosity of each target
            lum = df.groupby('EPIC').Lum_SED.first().values[0]
            lum = u.Quantity(lum).value
            #total observation time for each target
            #Mind that some targets are observed in multiple campaigns!
            time = df.groupby('EPIC').apply(obstime) #return yr
            #Eqn. (11):
            FA = (fl / lum / time).sum() / ffd.nstars * u.erg / u.erg * u.s / u.yr
            ffd.FA = FA.value
            return ffd

    def calculate_FR(self):
        '''
        Calculate the flaring rate FR in
        1/yr.
        '''
        ffd = deepcopy(self)
        ffd.FR = ffd.nflares / ffd.tot_obs_time
        ffd.FR_err = np.sqrt(ffd.nflares) / ffd.tot_obs_time
        return ffd

    def calculate_flaring_luminosity(self):
        '''
        Calculate the total flaring energy shedded
        per time in erg/s.
        '''
        ffd = deepcopy(self)
        if 'en' not in ffd.flares.columns:
            ffd.cutoff_energy_lower = 0.*u.erg
            logger.warning('No mode or energy cutoff set. Default to 0.0 erg.')
        if ffd.flares.en.values.shape[0]==0:
            ffd.flarelum = np.nan
            logger.warning('Could not calculate flaring luminosity. No values with energy available.')
            return ffd
        else:
            ffd.flarelum = ffd.flares.en.sum() / (ffd.tot_obs_time * 365.25 * 24. * 60 * 60)
        return ffd

    def calculate_flaring_fraction(self):
        '''
        Calculate the fraction of all stars in the sample
        that do flare.

        What's the uncertainty?
        '''
        ffd = deepcopy(self)
        nflaring = ffd.flares.EPIC.drop_duplicates().shape[0]
        ffd.flare_frac = float(nflaring)/float(ffd.nstars)
        return ffd

    def dump_pickle(self, add='general'):
        '''
        Save FFD object as a pickled file.

        Parameters:
        ------------
        add : 'general' or string
            specify why you saved this FFD

        Return:
        --------
        path to pickled file
         '''
        #wb = write + binary
        obj = 'FFD'
        if self.cutoff_energy is not None:
            if self.cutoff_energy.value > 0.:
                add = self.cutoff_energy
        path = 'saved_objects/{0}/{1}_{2:04d}_{3:04d}_{4}.p'.format(obj, self.name, int(self.Tmin.value), int(self.Tmax.value), add)
        pickle.dump(self, open(path, "wb" ) )
        return path

    def sample_ffd_solution(self, n=100, method='Bauke'):
        '''
        Sample random power laws with the same
        best-fit power law exponent, same minimum
        and maximum value, and of the same
        sample size.

        Parameters:
        -----------
        n : int
            Number of iterations
        Return:
        -------
        amean, aerr : float, float
            mean alpha value and mean uncertainty
            from n itearations passed to sol_alpha_mean
            and sol_alpha_mean_err
        '''
        methodmap = {'Bauke':FFD.powerlawfit,
                     'MK':FFD.powerlaw_Maschberger_and_Kroupa,
                     'fix2':FFD.fit_beta_to_powerlaw}
        assert self.flares is not None
        assert 'ed_rec' in self.flares
        assert self.alpha is not None
        alphas, alpha_errs = [], []
        fit = methodmap[method]
        size = self.flares.ed_rec.shape[0]
        minval = self.flares.ed_rec.min()
        maxval = self.flares.ed_rec.max()
        for i in range(n):
            data = generate_random_power_law_distribution(minval,
                                                          maxval,
                                                          -self.alpha+1.,
                                                          size=size)
            ffdt = FFD(ED=data)
            try:
                ffdt_b = fit(ffdt)
                alphas.append(ffdt_b.alpha)
                alpha_errs.append(ffdt_b.alpha_err)
                del ffdt_b, ffdt
            except ValueError:
                continue

        alphas = np.array(alphas)
        alpha_errs = np.array(alpha_errs)
        alphas= alphas[alphas != np.array(None)]
        alpha_errs= alpha_errs[alpha_errs != np.array(None)]
        if (len(alphas) > 0) & (len(alpha_errs) > 0):
            amean = np.mean(np.array(alphas))
            aerr = np.mean(np.array(alpha_errs))
            #write results to file
            res = pd.DataFrame({'a':alphas, 'ae': alpha_errs})
            path = ('ancillary/powerlaw_statistics/data/'
                    'sample_ffd_solution_{}_{}_{}_{}.csv'.format(self.name,
                                                                 int(self.Tmin.value),
                                                                 self.mode,
                                                                 timestamp()))
            res.to_csv(path, index=False)
            self.sol_alpha_mean = amean
            self.sol_alpha_err_mean = aerr
            return