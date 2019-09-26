# The Results class lives here
import pickle

import pandas as pd
import numpy as np
import astropy.units as u


from opencluster.opencluster import OpenCluster, TEFFBINS
from opencluster.ffd import FFD
from opencluster.analysis import *
from opencluster.bayes import BayesianFlaringAnalysis

from .__init__ import logger

class Results():
    '''
    Highest level object. Collection of OpenCluster
    objects that consist of FFD objects. Methods
    include short hands to analyse everything and also
    methods to compare between OCs.

    Attributes:
    ------------
    clusters : OpenCluster objects

    units : list
        list of units treated in the sample
        typically ['s', 'erg']
    results : DataFrame
        table with all the results from postprocessing
    analysed : bool
        Did you run :func:analyse_results()?
    cutoffs_intercomp : DataFrame
        table with the thresholds common to
        all the datasets in the sample so that you
        can compare across clusters and FFDs.
    cutoffs_intercomp_teff: DataFrame
        table with the thresholds common to
        all the datasets in the sample so that you
        can compare across clusters and FFDs in each
        Teff bin.
    '''
    def __init__(self, clusters=None, units=None, results=None,
                 analysed=False, cutoffs_intercomp=None,
                 cutoffs_intercomp_teff=None):
        self.clusters = clusters
        self.units = units
        self.OCs = dict()
        self.results = results
        self.analysed = analysed
        self.cutoffs_intercomp = cutoffs_intercomp
        self.cutoffs_intercomp_teff = cutoffs_intercomp_teff
        if ((clusters is not None) & (units is not None)):
            for cluster in clusters:
                self.OCs[cluster] = dict()
                for unit in units:
                    try:
                        OC = generate_OpenCluster(cluster, unit=unit)
                        OC = OC.generate_all_FFDs()
                        self.OCs[cluster][unit] = OC
                    except FileNotFoundError:
                        pass

    def extract_FFDs_by_Teff(self, unit):
        '''
        Aggregate FFDs by Teff bins and a unit.
        '''
        aggregated_ffds = dict()
        for Teff in TEFFBINS.Teff_min:
            aggregated_ffds[Teff] = []
            for cluster in self.clusters:
                try:
                    aggregated_ffds[Teff].append(self.OCs[cluster][unit].ffds[Teff])
                except KeyError:
                    pass
        return aggregated_ffds

    def analyse_results(self, **kwargs):
        '''
        Parameters:
        ------------
        kwargs : dict
            Keyword arguments to pass to :func:add_results_to_table()
        '''
        R = deepcopy(self)
        Us = {'s':u.s, 'erg':u.erg}
        for cluster, OCs in R.OCs.items():
            for unit in R.units:
                if unit in OCs.keys():
                    OC = OCs[unit]
                    df = pd.DataFrame(columns=COLUMNS)
                    for ffd in OC.ffds.items():
                        df = add_results_to_table(ffd[1], df, Us[unit], **kwargs)

                    df.to_csv('tables/{}_{}_{}.csv'.format(cluster,unit,timestamp()), index=False)
        R.analysed = True
        R = R.add_compiled_comparison_table_to_results(save=False, n=len(R.clusters))
        return R

    def add_compiled_comparison_table_to_results(self, **kwargs):
        '''
        Parameters:
        -----------
        kwargs : dict
            keyword arguments to pass to
            :func:compile_cluster_comparison_table
        '''
        R = deepcopy(self)
        df = compile_cluster_comparison_table(**kwargs)
        R.results = df
        return R

    def define_cutoff_for_entire_sample(self):
        '''
        Find the maximum cutoff value in all bins and
        clusters to compare across both attributes.
        '''
        R = deepcopy(self)
        R.cutoffs_intercomp = R.results.groupby('unit').cutoff_MK.max()
        return R

    def define_cutoff_for_Teffbins(self):
        '''
        Find the maximum cutoff value in all bins and
        clusters to compare across both attributes.
        '''
        R = deepcopy(self)
        R.cutoffs_intercomp_teff = R.results.groupby(['unit','Tmin']).cutoff_MK.max()
        return R

    def add_FA_FR_to_analysis(self):
        '''
        Add values for FR and FA to the analysis table.
        Update results attribute.
        '''
        R = deepcopy(self)
        cols = ['unit','cluster','Tmin','FA','FR', 'FR_err', 'flarelum']
        res = pd.DataFrame(columns=cols)
        for cluster in R.clusters:
            for un in R.units:
                for Tmin, ffd in R.OCs[cluster][un].ffds.items():
                    try:
                        ffd = ffd.calculate_FR()
                        if un == 'erg':
                            ffd = ffd.calculate_FA()
                            ffd = ffd.calculate_flaring_luminosity()
                        else:
                            ffd.FA = np.nan
                            ffd.flarelum = np.nan
                        res = res.append(dict(zip(cols,[un, cluster,
                                                        Tmin,ffd.FA,
                                                        ffd.FR, ffd.FR_err,
                                                        ffd.flarelum])),
                                         ignore_index=True)
                    except (ValueError, AttributeError):
                        logger.debug('Cluster {} failed to '
                                     'produce FA and/or FR. '
                                     'Skipped.'.format(ffd.h_name))
                        continue
        data = R.results.merge(res, how='left', on=['unit','cluster','Tmin'])
        R.results = data
        return R

    def dump_to_pickle(self, add='Results'):
        '''
        Save Results object as a pickled file.

        Parameters:
        ------------
        add : 'general' or string
            specify why you saved these Results

        Return:
        --------
        path to pickled file
        '''

        obj = 'Results'
        path = 'saved_objects/{0}/{1}_{2}.p'.format(add, obj, timestamp())
        pickle.dump(self, open(path, "wb"))#wb = write + binary
        return path

    def cutoff_at_shared_threshold(self, how='all'):
        '''
        Cutoff all FFDs according to shared thresholds.
        Either Teff-bin-wise or for the entire sample.

        Parameters:
        ------------
        how : 'all' or 'teff'
            Either use R.cutoffs_intercomp or
            R.cutoffs_intercomp_teff for the cutoffs

        '''
        R = deepcopy(self)
        if how=='all':
            threshs = R.cutoffs_intercomp
        elif how=='teff':
            threshs = R.cutoffs_intercomp_teff
        for cluster in R.clusters:
            for un in R.units:
                for Tmin, ffd in R.OCs[cluster][un].ffds.items():
                    if ffd.ED is not None:#need to have some data
                        try:
                            id_ = {'all':un, 'teff':(un,Tmin)}
                            match = id_[how] in threshs.index#Teff threshold should exist
                            if match:
                                if how == 'all':
                                    cutoff = threshs[un]
                                elif how == 'teff':
                                    cutoff = threshs[un][Tmin]
                                ffd.cutoff = cutoff * u.Unit(un)
                                logger.debug('Cutoff for {} matched to {}'.format(ffd.h_name, ffd.cutoff))
                        except ValueError:
                            continue
        return R

    def fit_alpha_and_eps_bayesian(self, cluster, Teff, loglikelihood, mined=10000, deltaT=365.25, alpha_prior=2.):
        '''Determine the probability eps of a flare at >= mined
        to occur within deltaT for one Teff bin in a given cluster.
        Only works with ED so far.

        For example:

            # read in likelihood functions, prior distributions

            from opencluster.bayes import calculate_joint_posterior_distribution, uninformative_prior, gaussian_prior

            # get the analysis toolkit

            from opencluster.bayes import BayesianFlaringAnalysis

            # mix and match priors into your final posterior

            def loglikelihood(theta, *args):

                def prior(x):
                    return uninformative_prior(x, 1.25, 2.25)
                return calculate_joint_posterior_distribution(theta, *args, prior)

            BFA.fit_alpha_and_eps_bayesian("praesepe", 3000., loglikelihood)

        Parameters:
        -----------
        cluster : str
            cluster name
        Teff : float
            lower Teff bin edge
        mined : float
            default 10 ks. Lower flaring energy to predict.
        deltaT : float
            default 365.25 days. Prediction time span.
        alpha_prior : float
            default 2. prior best guess for power law exponent
        '''

        unit = "s" # maybe expand to erg later
        OCs_inTeffs = self.OCs[cluster][unit]

        # Intermediate data structures

        results = self.results
        ffd = OCs_inTeffs.ffds[Teff]
        inputs = ffd.flares[["tstart","ed_rec"]]

        # Data

        params = results[(results.cluster == cluster) &
                         (results.Tmin == Teff) &
                         (results.unit == unit)].iloc[0]
        Mprime = params.nflares
        Tprime = ffd.tot_obs_time * 365.25 # convert to days
        threshed = ffd.cutoff.value
        events = inputs.ed_rec.values

        # Sanity checks

        assert threshed == params.cutoff_MK # Check for unexplicable cutoff made after initial ML analysis

        # define prior propbability from ML estimates:

        if alpha_prior == 2.:
            beta = params.beta2 / 365.25
        else:
            beta = params.beta_MK / 365.25

        rate_prior = (beta / np.abs(alpha_prior - 1.) *
                      np.power(mined, -alpha_prior + 1.)) # evaluate cumulative FFD fit at mined
        eps_prior = 1 - np.exp(-rate_prior * deltaT) # Poisson probability

        # Do bayesian flaring analysis:

        # Init the analysis suite:
        BFA = BayesianFlaringAnalysis(mined=mined, Tprime=Tprime, deltaT=deltaT,
                                      alpha_prior=alpha_prior, eps_prior=eps_prior,
                                      threshed=threshed, Mprime=Mprime, M=Mprime,
                                      events=events, loglikelihood=loglikelihood)

        # Run MCMC to sample the posterior distribution
        BFA.sample_posterior_with_mcmc()

        # Calculate percentiles:
        BFA.calculate_percentiles()

        ffd.BFA = {"cluster":cluster,
                   "Teff":Teff,
                   "eps": BFA.eps_posterior,
                   "alpha": BFA.alpha_posterior}

        return