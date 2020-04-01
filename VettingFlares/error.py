import numpy as np
import pandas as pd

from altaipony.altai import find_iterative_median
import copy

def refine_detrended_flux_err(flcd, mask_pos_outliers_sigma=2.5, 
                              std_rolling_window_length=15, pad=25):
    """Attempt to recover a good estimate of the ligh curve noise.
    Start out from a simple standard deviation of the flux.
    Then filter out outliers above `mask_pos_outliers_sigma`.
    Apply rolling window standard deviation on the filtered array.
    Calculate a mean standard deviation from the result.
    Fill in this mean into the masked values.

    Parameters:
    -----------
    flcd : de-trended FlareLightCurve

    mask_pos_outliers_sigma : float
        sigma value above which to mask positive outliers
    std_rolling_window_length : int
        rolling window length for standard deviation calculation
    pad : int
        How many values to pad-mask around positive outliers.

    Return:
    --------
    FlareLightCurve with refined `detrended_flux_err` attribute.

    """

    # start with a first approximation to std
    flcd.detrended_flux_err[:] =  np.nanstd(flcd.detrended_flux)

    # and refine it:
    flcd = find_iterative_median(flcd)

    filtered = copy.deepcopy(flcd.detrended_flux)

    # mask strong positive outliers so that they don't add to std
    filtered[flcd.detrended_flux - flcd.it_med > mask_pos_outliers_sigma * flcd.detrended_flux_err] = np.nan

    # apply rolling window std
    flcd.detrended_flux_err[:] = pd.Series(filtered).rolling(std_rolling_window_length, min_periods=1).std()

    # set std to mean value if calculation fails to inf
    meanstd = np.nanmean(flcd.detrended_flux_err)

    # pad the excluded values not to create spikes of high error around flares
    isin = np.invert(np.isfinite(flcd.detrended_flux_err))
    x = np.where(isin)[0]
    for i in range(-pad, pad+1):
        y = x + i
        y[np.where(y > len(isin) - 1)] = len(isin) - 1
        isin[y] = True

    x = np.where(isin)[0]
    flcd.detrended_flux_err[x] = meanstd

    return flcd