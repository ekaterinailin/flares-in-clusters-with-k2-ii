import glob
import os
import numpy as np
import pandas as pd
from copy import deepcopy

from astropy import units as u

from opencluster.opencluster import OpenCluster
from opencluster.utils import timestamp
from opencluster.paperplots import save_PP_plots, COLORS_LIST, MARKERS_LIST
from opencluster import PACKAGEDIR

from .__init__ import logger


COLUMNS = ['h_cluster','cluster', 'Tmin', 'Tmax', 'cutoff_MK',
           'alpha_MK', 'alpha_MK_err','beta_MK', 'R1s_MK', 'truncated_MK',
           'alpha_Bauke', 'alpha_Bauke_err','beta_Bauke', 'R1s_Bauke', 'truncated_Bauke',
           'alpha_linfit', 'alpha_linfit_err', 'beta_linfit', 'R1s_lin', 'truncated_lin',
           'beta2','beta2_err', 'R1s_2', 'truncated_2',
           'nstars', 'unit', 'nflares', 'tot_obs_time',
           'flare_frac', 'beta_Bauke_err', 'beta_MK_err',
           'sol_alpha_mean_2','sol_alpha_err_mean_2',
           'sol_alpha_mean_mk','sol_alpha_err_mean_mk',
           'sol_alpha_mean_bauke','sol_alpha_err_mean_bauke',]
CLUSTERS = pd.read_csv('{}/clusters/cluster_parameters_merged.csv'.format(PACKAGEDIR))




def generate_OpenCluster(cluster, unit=None, clusters=CLUSTERS, path=""):
    '''
    Generate an OpenCluster object using the
    parameter table stored in clusters.

    Parameters:
    -------------
    cluster : string
        cluster machine friendly name
    clusters : DataFrame
        table of cluster parameters
    '''
    if unit is None:
        logger.error('Please specify a unit to work with')
    cparams = clusters.loc[(clusters.cluster == cluster),:].iloc[0]
    h_cluster = cparams.h_cluster
    flares = pd.read_csv('{}/flares/{}_flares.csv'.format(path, cluster))
    stars = pd.read_csv('{}/luminosities/{}_luminosities.csv'.format(path, cluster))
    OC = OpenCluster(cluster=cluster, h_cluster=h_cluster,
                     age=cparams['age (Myr)'], u_age_high=cparams.u_age_high, u_age_low=cparams.u_age_low,
                     feh=cparams.FeH, u_feh=cparams.u_feh,
                     stars=stars, flares=flares, unit=unit)
    OC.merge_flares_stars()
    #Select stars and flares with Teff below a maximum value specified in "clusters"
    g = lambda x: u.Quantity(x).value
    OC.stars = OC.stars[~OC.stars.Teff_median.apply(g).isnull()]
    OC.stars = OC.stars[OC.stars.Teff_median.apply(g) < cparams.Teffmax]
    OC.flares_on_stars = OC.flares_on_stars[~OC.flares_on_stars.Teff_median.apply(g).isnull()]
    OC.flares_on_stars = OC.flares_on_stars[OC.flares_on_stars.Teff_median.apply(g) < cparams.Teffmax]
    print(OC.flares_on_stars.head(), cparams)
    OC.flares_on_stars = OC.flares_on_stars.drop_duplicates(subset=['EPIC','ed_rec','cstart','cstop'])
    print(OC.flares_on_stars)
    OC.flares = OC.flares_on_stars

    return OC



def find_cutoff_by_KS(FFD, unit, factor=1.02):
    '''
    Iteratively cut off one value at a time
    and calculate Maschberger and Kroupa power
    law exponent, then do KS test. Repeat until
    KS test confirms that the power law is a good
    idea. Finally, as rule of thumb,
    double the calculated threshold.

    Return:
    -------
    astropy Quantity
    '''
    logger.disabled = True #turn off logging
    ffd = deepcopy(FFD)
    ffd.cutoff = 0. * unit
    ffd.KS = False
    while ffd.KS == False:
        if ffd.mode == 'ED':
            assert unit == u.s
            ffd.cutoff = np.min(ffd.ED) * unit
            if len(ffd.ED) < 5:
                ffd.cutoff = np.nan * unit
                break
        elif ffd.mode == 'energy':
            assert unit == u.erg
            if len(ffd.energy) < 5:
                ffd.cutoff = np.nan * unit
                break
            ffd.cutoff = np.min(ffd.energy) * unit
            if len(ffd.energy) < 5:
                ffd.cutoff = np.nan * unit
                break
        else:
            raise AttributeError('mode attribute not properly set.')
        ffd = ffd.powerlaw_Maschberger_and_Kroupa() #compute alpha for the last time
        ffd = ffd.is_powerlaw_by_KS()
    ffd.cutoff = ffd.cutoff * factor #increase the threshold rule of thumb
    logger.disabled = False #turn on logging again
    return ffd.cutoff

def analyse_with_linfit(FFD):
    '''
    Apply full analysis with graphical method.
    '''
    try:
        ffd_lin = deepcopy(FFD)
        ffd_lin = ffd_lin.linefit()
        ffd_lin = ffd_lin.calculate_R1s()
        ffd_lin = ffd_lin.is_powerlaw_truncated()
    except ValueError:
        ffd_lin.alpha = np.nan
        ffd_lin.alpha_err = np.nan
        ffd_lin.beta = np.nan
        ffd_lin.beta_err = np.nan
        ffd_lin.alpha_mode = 'linear fit'
        ffd_lin.R1s = np.nan
        ffd_lin.truncated = False
    return ffd_lin

def analyse_with_Bauke(ffd):
    '''
    Apply full analysis to FFD fitted with Bauke
    method.
    '''
    ffd_bauke = deepcopy(ffd)
    try:
        ffd_bauke = ffd_bauke.powerlawfit()
        ffd_bauke = ffd_bauke.fit_beta_to_powerlaw()
        ffd_bauke = ffd_bauke.calculate_R1s()
        ffd_bauke = ffd_bauke.is_powerlaw_truncated()
    except ValueError:
        ffd_bauke.alpha = np.nan
        ffd_bauke.alpha_err = np.nan
        ffd_bauke.beta = np.nan
        ffd_bauke.beta_err = np.nan
        ffd_bauke.alpha_mode = 'Bauke'
        ffd_bauke.R1s = np.nan
        ffd_bauke.truncated = False
    return ffd_bauke

def analyse_with_alpha2(ffd):
    '''
    Apply full analysis after fixing alpha=2.
    '''
    try:
        ffd_2 = deepcopy(ffd)
        ffd_2.alpha = 2.0
        ffd_2.alpha_err = 0.
        ffd_2 = ffd_2.fit_beta_to_powerlaw()
        ffd_2 = ffd_2.calculate_R1s()
        ffd_2 = ffd_2.is_powerlaw_truncated()
    except ValueError:
        ffd_2.alpha_err = 0.
        ffd_2.alpha = 2.
        ffd_2.beta = np.nan
        ffd_2.beta_err = np.nan
        ffd_2.R1s = np.nan
        ffd_2.truncated = False
    ffd_2.alpha_mode = 'fix2'
    return ffd_2

def analyse_with_MK(ffd):
    '''
    Apply full analysis of FFD using Mascherger and Kroupa
    methods.
    '''
    try:
        ffd_mk = deepcopy(ffd)
        ffd_mk = ffd_mk.powerlaw_Maschberger_and_Kroupa()
        ffd_mk = ffd_mk.fit_beta_to_powerlaw()
        ffd_mk = ffd_mk.calculate_R1s()
        ffd_mk = ffd_mk.is_powerlaw_truncated()
    except (ValueError, ZeroDivisionError):
        ffd_mk.alpha = np.nan
        ffd_mk.alpha_err = np.nan
        ffd_mk.beta = np.nan
        ffd_mk.beta_err = np.nan
        ffd_mk.alpha_mode = 'MK'
        ffd_mk.R1s = np.nan
        ffd_mk.truncated = False
    return ffd_mk

def compile_cluster_comparison_table(save=True, n=4):
    '''
    Combine the results from the cluster analysis
    into one big table.

    Parameters:
    ------------
    save : bool
        If True, will save table to tables/ folder
    n : int
        number of clusters

    Return:
    -------
    DataFrame with all results
    '''
    liste = []
    for unit in ['s','erg']:
        list_of_files = glob.glob('tables/*_{}_*'.format(unit)) # * means all if need specific format then *.csv
        latest_files = sorted(list_of_files, key=os.path.getctime)
        [liste.append(file) for file in latest_files[-n:]]
    data = pd.DataFrame(columns=COLUMNS)
    for x in liste:
        data = data.append(pd.read_csv(x), ignore_index=True)
    data = data.merge(CLUSTERS.drop('cluster',axis=1),how='left',left_on='h_cluster',right_on='h_cluster')
    data.Tmin = data.Tmin.apply(u.Quantity).apply(lambda x: x.value)
    data.Tmax = data.Tmax.apply(u.Quantity).apply(lambda x: x.value)
    data.cutoff_MK = data.cutoff_MK.apply(u.Quantity).apply(lambda x: x.value)
    if save==True:
        data.to_csv('tables/cluster_comparison_{}.csv'.format(timestamp()), index=False)
    return data

def add_results_to_table(ffd, df, unit, ppplots=False, sample_solution=False):
    '''
    Analyse the ffd with a host of methods,
    add results to some shared DataFrame.

    Parameters:
    ------------
    ffd : FFD object
        FFD to analyse
    df : DataFrame
        Table to add the results to
    unit : astropy.unit
        u.erg or u.s.
    ppplots : bool
        if True will generate and save PP-plots
        for all power law fits
    sample_solution : bool
        if True will sample power law dstributions
        from the best fit parameters to probe
        the uncertainties on these results

    Return:
    -------
    Updated df
    '''
    if ffd.flares is not None:
        try:
            cutoff = find_cutoff_by_KS(ffd, unit)
            print('Cutoff FFD {} ({}-{}) at {}'.format(ffd.h_name, ffd.Tmin, ffd.Tmax, cutoff))
            if np.isnan(cutoff):
                return df
            else:
                ffd.cutoff = cutoff
                ffd = ffd.calculate_flaring_fraction()
                ffd_mk = analyse_with_MK(ffd)
                ffd_2 = analyse_with_alpha2(ffd)
                ffd_bauke = analyse_with_Bauke(ffd)
                ffd_lin = analyse_with_linfit(ffd)
                print(ffd_mk.beta_err, ffd_bauke.beta_err)
                #validate results
                assert ffd_mk.alpha_mode == 'MK'
                assert ffd_bauke.alpha_mode == 'Bauke'
                assert ffd_lin.alpha_mode == 'linear fit'
                assert ffd_2.alpha == 2.
                assert ffd_2.alpha_mode == 'fix2'
                # additional outputs
                if sample_solution == True:
                    for ffd in [ffd_bauke,ffd_mk]:
                        ffd.sample_ffd_solution(n=100, method=ffd.alpha_mode)
                if ppplots == True:
                    save_PP_plots([ffd_2,ffd_bauke,ffd_mk])

                values = [ffd_mk.h_name,ffd_mk.name,ffd.Tmin,ffd.Tmax,ffd.cutoff,
                          ffd_mk.alpha, ffd_mk.alpha_err, ffd_mk.beta, ffd_mk.R1s, ffd_mk.truncated,
                          ffd_bauke.alpha, ffd_bauke.alpha_err, ffd_bauke.beta, ffd_bauke.R1s, ffd_bauke.truncated,
                          ffd_lin.alpha, ffd_lin.alpha_err, ffd_lin.beta, ffd_lin.R1s, ffd_lin.truncated,
                          ffd_2.beta, ffd_2.beta_err, ffd_2.R1s, ffd_2.truncated,
                          ffd.nstars, str(unit), ffd.nflares, ffd.tot_obs_time,
                          ffd.flare_frac, ffd_bauke.beta_err, ffd_mk.beta_err,
                          ffd_2.sol_alpha_mean, ffd_2.sol_alpha_err_mean,
                          ffd_mk.sol_alpha_mean, ffd_mk.sol_alpha_err_mean,
                          ffd_bauke.sol_alpha_mean, ffd_bauke.sol_alpha_err_mean,]

                df = df.append(dict(zip(COLUMNS,values)), ignore_index=True)
        except ValueError:
            return df
    return df
