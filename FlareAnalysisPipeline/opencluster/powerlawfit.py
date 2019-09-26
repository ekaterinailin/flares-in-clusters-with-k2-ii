import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import zeta as Hurwitz
from scipy.optimize import fmin
from .__init__ import logger

def ML_powerlaw_estimator(EDs, a):
    '''
    Power law ML estimator from
    Maschberger and Kroupa (2009),
    formula (9).

    Parameters:
    -----------
    data : Series or np.array
        data that is suspected to follow
        a power law relation
    '''
    if np.array(a <= 1.).any():
        raise ValueError('Power law exponent must be >1.')
    n = len(EDs)
    if n == 0:
        raise ValueError('No data.')
    Y = EDs.min()
    if Y < 0:
        raise ValueError('Negative value encountered in data.')
    a = de_bias_alpha(n, a)
    Yexp = (np.power(Y,1-a))
    T = np.log(EDs).sum()
    Z = de_biased_upper_limit(EDs, a)
    Zexp = (np.power(Z,1-a))

    return n / (a - 1) + n * ((Zexp * np.log(Z) - Yexp * np.log(Y)) / (Zexp - Yexp)) - T


def de_biased_upper_limit(data, a):
    '''
    De-biases the upper limits for a
    ML power law exponent estimator.
    Uses formular (13) (and (14)) from
    Maschberger and Kroupa (2009).

    Parameters:
    -----------
    data : Series or array
        data that is suspected to follow
        a power law relation
    a : float or array of floats
        quasi de-biased ML estimator for alpha
        (de_bias_alpha before inserting here!)

    Returns:
    ---------
    Quasi de-biased upper limit.
    '''
    if len(data) == 0:
        raise ValueError('No data.')
    if (data < 0).any():
        raise ValueError('Negative values '
                         'encountered in data.')
    Xn = data.max()
    X1 = data.min()
    if Xn == X1:
        raise ValueError('Data range is zero.')
    n = len(data)
    G = (1. - a) * np.log(Xn / X1)#(14)
    base = 1.+ (np.exp(G) - 1.) / n
    exponent = 1. / (1. - a)
    return Xn * np.power(base, exponent)

def de_bias_alpha(n, alpha):
    '''
    De-biases the power law value
    according to Maschberger and Kroupa (2009),
    formula (12).

    Paramaters:
    ------------
    n : int
        Size of the data
    alpha : float or array of floats
        Power law exponent value from ML estimator

    Returns:
    -----------
    quasi de-biased ML estimator for alpha
    '''
    if np.array(np.isnan(n) | np.isnan(np.array(alpha)).any()):
        raise ValueError('de_bias_alpha: one or '
                         'both arg(s) is/are NaN')
    return (alpha - 1.) * n / (n - 2) + 1.

def lnHurwitz_Ceta(gamma,kmin,kmax):
    '''
    Calculate the modified natural log of
    the Hurwitz Ceta function from Bauke 2007.

    Parameters:
    -----------
    kmin, kmax : float or int
        min/max energy in sample that is distributed
        according to a power with:
    gamma : float
        power law exponent to a*x^(-gamma)
    '''
    if (np.asarray((gamma <= 1.)).any() | (kmin <= 0.) | (kmax <= 0.)):
        logger.warning('Hurwitz Ceta function could not be used.'
                       ' Test alpha was {}'.format(gamma))
        raise ValueError('Cannot calculate ln Hurwitz because'
                         ' either alpha or the boundaries are'
                         ' not applicable to scipy implementation'
                         ' of Hurwitz-Ceta.')
    if kmin == kmax:
        raise ValueError('Cannot calculate ln Hurwitz because'
                         ' minimum and maximum energies are the'
                         ' same! Check your sample size.')
    res = Hurwitz(gamma, kmin) - Hurwitz(gamma, kmax)
    if (np.asarray(res) > 0).all():
        return np.log(res)
    else:
        raise ValueError('Cannot calculate ln of a negative number'
                         ' or a NaN value.')

def logL_powerlaw(alpha, data):
    '''
    Calculate the term that has to be minimized
    to find the best-fit power law exponent alpha.

    Parameters:
    -----------
    alpha : Power law exponent
        needs to be > 1.
    data : array of ints
        distribution of measurements
    Return:
    --------
    -ln(likelihood) -> minimize!

    '''

    L = []
    # Formula 25 in Bauke 2007:
    if type(data) != np.ndarray:
        raise TypeError('The data are not in numpy array format. Check your input.')
    else:
        A = - data.shape[0] * lnHurwitz_Ceta(alpha, data.min(), data.max()) -1. * alpha * np.sum(np.log(data))
        return A

def powerlawfit_Bauke2007(data, minexp=1.+1e-8, maxexp=4.):
    '''
    Calculate the power law exponent using MLE
    approach from Bauke 2007.

    Parameters:
    -----------
    data : array of ints
        Distribution of measurements in no particular order.
    minexp : float >1.
        Minimum power law exponent
    maxexp : float >1.
        Maximum power law exponent

    Return:
    ----------
    Best fit value for power law exponent.
    '''
    N = 100000
    x = np.linspace(minexp, maxexp, N)
    #print('Power law exponent maximum accuracy: {}'.format((maxexp-minexp)/N))
    y = -1.*logL_powerlaw(x, data)
    yf = np.nan_to_num(y)
    yf[np.where(yf==0.)] = 1e10
    m = x[np.argmin(yf)]
    if m==minexp:
        raise ValueError('The result may be wrong because the result'
                         ' for alpha is at the boundary of considered'
                         ' exponents. Set a different range of exponents'
                         ' to check.')
    return m
