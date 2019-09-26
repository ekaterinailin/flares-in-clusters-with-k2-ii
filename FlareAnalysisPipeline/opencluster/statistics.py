# All the statistical tests implemented here are extracted from Maschberger and Kroupa (2009)
# with some help from Wikipedia. Read this project's README for more explanations and nice plots
# to find out why these are used here.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .utils import generate_random_power_law_distribution

def calculate_average_number_of_exceeding_values(ffd, n, **kwargs):
    '''
    Parameters:
    -----------
    ffd : FFD object

    n : int
        number of samples to average
    kwargs : dict
        Keyword arguments to pass to
        :func:calculate_number_of_exceeding_values

    Returns:
    --------
    (mean, std) : (float, float)
        average number number of exceeding values
        and standard deviation
    '''
    if (ffd.ED is None) & (ffd.flares is not None):
        data = ffd.flares.ed_rec.apply(lambda x: x.value)
    elif ffd.ED is not None:
        data = ffd.ED
    else:
        raise ValueError('No data in FFD. Check your inputs.')
    alpha = ffd.alpha
    assert alpha is not None
    assert data is not None
    exceedance_statistic = [_calculate_number_of_exceeding_values(data, alpha, **kwargs) for i in range(n)]
    exceedance_statistic = np.array(exceedance_statistic)
    return np.nanmean(exceedance_statistic), np.nanstd(exceedance_statistic)

def _calculate_number_of_exceeding_values(data, alpha, maxlim=1e8, **kwargs):
    '''
    Helper function that mimicks data similar
    to the observations (same alpha and size)
    and returns a sample from an untruncated
    distribution. The number of values that
    exceeds the maximum in the actual data is
    returned.

    Parameters:
    -----------
    data : array
        observed values
    alpha : float
        best-fit power law exponent to the data
    maxlim : float > 1.
        factor to simulate an untruncated
        version of the given power law
        distribution
    kwargs : dict
        Keyword arguments to pass to
        :func:generate_random_power_law_distribution

    Returns:
    --------
    int : number of exceeding values
    '''
    pdist = generate_random_power_law_distribution(np.min(data),
                                                   np.max(data) * maxlim,
                                                   -alpha+1,
                                                   size=data.shape[0],
                                                   **kwargs)

    if np.isnan(pdist).any():
        raise ValueError('Fake power law distribution for the'
                         ' exceedance test could not be generated.'
                         ' Check your inputs.')
    return len(np.where(pdist > np.max(data))[0])

def stabilised_KS_statistic(data, **kwargs):
    '''
    Calculate the stabilised KS statistic
    from Maschberger and Kroupa 2009, Eqn. (21)
    orginally from Michael 1983, and Kimber 1985.

    Parameters:
    --------------
    data : array
        observed values that are suspected
        to follow a power law relation
    kwargs : dict
        Keyword arguments to pass to
        :func:calculate_cumulative_powerlaw_distribution
    Return:
    --------
    float - stablised KS statistic
    '''
    sorted_data = np.sort(data)
    pp = calculate_cumulative_powerlaw_distribution(sorted_data, **kwargs)
    y = (np.arange(1, len(pp) + 1) - .5) / len(pp)
    argument = (_apply_stabilising_transformation(y)
                - _apply_stabilising_transformation(pp))
    return np.max(np.abs(argument))

def calculate_KS_acceptance_limit(n, sig_level=0.05):
    '''
    Above this limit we must reject the null-hypothesis.
    In our context, this is the hypothesis that the dis-
    tribution follows a given power law.

    Parameters:
    -----------
    n : int
        sample size
    sig_level : 0 < float < 1.
        significance level
    '''
    if ((sig_level >= 1.) | (sig_level <= 0.)):
        raise ValueError('Pass a valid significance level.')
    if n == 0:
        raise ValueError('No data to calculate KS_acceptance limit.')
    elif ((n <= 35) & (n > 0)):
        t = (pd.read_table('opencluster/static/KS_leq_35_values.csv',
                           delimiter='|', skiprows=1, header=None,
                           names=['n',.9,.95,.99])
             .set_index('n')
             .astype(float))
        return t.loc[n, 1 - sig_level]
    elif n > 35:
        return np.sqrt(-.5 * np.log((sig_level) / 2.)) / np.sqrt(n)

def _apply_stabilising_transformation(u):
    '''
    Applies the right-tail stabilising
    transformation from Kimber 1985 to
    a potentially power law distributed
    sample. Eq. 19 in Maschberger and Kroupa 2009.

    Used in :func:stabilised_KS_statistic

    Parameters:
    ------------
    u : array
        sample

    Returns:
    -----------
    array : stabilised sample
    '''
    assert (np.array(u) >= 0).all() #validate input for sqrt
    S0 = 2. / np.pi * np.arcsin(np.sqrt(u))
    return 2. * S0 * (.5 + .5 * u) -1.


def calculate_cumulative_powerlaw_distribution(data, alpha=None, truncated=True):
    '''
    Calculates the cumulative powerlaw distribution
    from the data, given the best fit power law exponent
    for y(x) ~ x^(-alpha).
    Eq. (2) in Maschberger and Kroupa 2009.

    Parameters:
    -----------
    data : array
        observed values that are suspected
        to follow a power law relation, sorted in
        ascending order
    alpha : float
        best-fit power law exponent

    Returns:
    ---------
    array : cumulative distribution
    '''
    if alpha <=1.:
        raise ValueError('This distribution function is only'
                         ' valid for alpha > 1., see also '
                         'Maschberger and Kroupa 2009.')
    data = np.sort(data)
    def expa(x, alpha):
        return np.power(x, 1. - alpha)
    if truncated == True:
        CDF = ((expa(data, alpha) - expa(np.min(data), alpha))
              / (expa(np.max(data), alpha) - expa(np.min(data), alpha)))
    elif truncated == False:
        CDF = 1. - expa(data / np.min(data), alpha)
    #fix a -0. value that occurs as the first value
    CDF[np.where(CDF==0.)[0]] = 0.
    return CDF