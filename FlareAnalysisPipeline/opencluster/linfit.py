import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import zeta as Hurwitz
from scipy.optimize import fmin
import astropy.units as u

def logL_lin(X,*args):
    '''
    Calculate the term that has to be minimized
    to find the best-fit line parameters in X.

    Parameters:
    -----------
    X - (b, theta=arctan(m)) intercept and slope of line fit
    *args - (x, y, sigy, sigx, sigxy)


    Return:
    --------
    -ln(likelihood) -> minimize it!

    '''

    b, theta = X
    L = []
    Lx, Ly, Lsigy, Lsigx, Lsigxy = args

    for (x, y, sigy, sigx, sigxy) in list(zip(Lx, Ly, Lsigy, Lsigx, Lsigxy)):

        c = np.cos(theta)
        s = np.sin(theta)
        B = ((-1.) * (s * x)) + (c * y) - (b * c)
        A = (sigx * (s**2)) - (2 * sigxy * s * c) + (sigy * (c**2))
        L.append(0.5 * (B**2) / A)
    if np.isnan(sum(L)):
        raise ValueError('Likelihood function got a NaN output. Check your inputs.')
        return
    else:
        return sum(L)

def linfit(data):

    '''
    Fits a linear function to a data set with errors in x and y,
    calculates errors on best fit parameters using jackknife algorithm


    Parameters:
    ------------
    data - DataFrame with x, y, sigx, sigy, rho (correlation factor)

    Return:
    (m_mean, msig) - slope with uncertainty
    (b_mean, bsig) - intercept with uncertainty
    '''

    b_, theta_ = [], []
    data = data.applymap(float)
    data['sigxy'] = data.sigy * data.sigx * data.rho
    data = data.dropna(how='any', subset=['x','y','sigy'])
    mi, ma = data.x.min(), data.x.max()
#  t = np.arange(mi, ma, (ma - mi) / 200.)
    #jackknife loop
    for id_ in data.index.values:
        d = data.drop(id_)
        b, theta = fmin(logL_lin, [60,1.5],
                        args=(d.x, d.y, d.sigy, d.sigx, d.sigxy),
                        disp=0)
#        lin = b + np.tan(theta) * t - data.norm.iloc[0]
        b_.append(b - data.norm.iloc[0])
        theta_.append(theta)

    N = data.shape[0]
    b_, theta_ = np.asarray(b_), np.asarray(theta_)
    m_ = np.tan(theta_)
    m_mean = m_.mean()
    b_mean = b_.mean()
    sigm = np.sqrt( (N-1) / N * ( (m_ - m_mean)**2 ).sum() )
    sigb = np.sqrt( (N-1) / N * ( (b_ - b_mean)**2 ).sum() )
    return (m_mean, sigm), (b_mean, sigb)

def define_linfit_arrays(df, mode):
    '''
    Calculate the cumulative distribution from flare
    EDs or energies.

    Parameters:
    --------------
    df : DataFrame
        Contains an ed_rec and a Lum_Kepler column
    mode : string -  'ED' or 'energy'
        Define the output unit
    aggregate : False or bool
        Option to bin flares with identical energies

    Return:
    -------
    a :
        energy or ED corresponding to each value in freq
    freq :
        Cumulative frequency of flares with energy >= the
        corresponding value in a. Units of 1/unit(obs_duration)
    unit : astropy unit
        unit of a
    '''
    if mode == 'ED':
        a = np.sort(df.ed_rec.values)
    elif mode == 'energy':
        energy = df.ed_rec*df.Lum_Kepler
        a = np.sort(energy.values)
    a = a[~np.isnan(a)]
    freq = (np.arange(len(a), 0, -1) )
    #aggregate flares with identical energies
    return a, freq