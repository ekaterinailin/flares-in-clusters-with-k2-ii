import pandas as pd
import numpy as np
import pytest
import warnings
import astropy.units as u
from collections import OrderedDict
from ..ffd import FFD
from ..utils import *

def test_timestamp():
    pass

def test_mag_err_to_rel_flux_err():
    assert mag_err_to_rel_flux_err(0.) == 0.
    assert mag_err_to_rel_flux_err(1.) == pytest.approx(1.51188643)
    assert np.isnan(mag_err_to_rel_flux_err(np.nan))
    assert np.isnan(mag_err_to_rel_flux_err(np.nan))

def test_rel_flux_err_to_mag_err():
    assert rel_flux_err_to_mag_err(0.) == 0.
    assert rel_flux_err_to_mag_err(1.51188643) == pytest.approx(1.)

def test_remove_bad_2MASS_photometry():
    testdf = pd.DataFrame({'Qfl_2MASS' : ['AAA',np.nan,'DEA','DDD', 'CAB'],
                           'J_2MASS' : [12., np.nan, 23, 24, 25],
                           'H_2MASS' : [12.5, 34, 23.5, 24.5, 25.5],
                           'K_2MASS' : [12.1, np.nan, 23.1, 24.1, 25.1],})
    resultdf = pd.DataFrame({'Qfl_2MASS': {0: 'AAA', 1: np.nan, 2: 'DEA', 3: 'DDD', 4: 'CAB'},
     'J_2MASS': {0: 12.0, 1: np.nan, 2: np.nan, 3: np.nan, 4: np.nan},
     'H_2MASS': {0: 12.5, 1: np.nan, 2: np.nan, 3: np.nan, 4: 25.5},
     'K_2MASS': {0: 12.1, 1: np.nan, 2: 23.1, 3: np.nan, 4: np.nan}})
    assert resultdf.equals(remove_bad_2MASS_photometry(testdf))

def test_remove_bad_PS1_photometry():
    testdf = pd.DataFrame({'gFlags_PS1' : [8,np.nan,40],
                           'rFlags_PS1' : [0,256+64,40],
                           'iFlags_PS1' : [132,16,40],
                           'zFlags_PS1' : [16,4,40],
                           'yFlags_PS1' : [8,4,40],
                           'g_PS1' : [12., 23, np.nan],
                           'r_PS1' : [12.5, 34, 23.5,],
                           'i_PS1' : [12.1, np.nan, 25.1],
                           'z_PS1' : [12.1, np.nan, 25.1],
                           'y_PS1' : [12.1, np.nan, 25.1],
                           'e_g_PS1' : [1., .3, np.nan],
                           'e_r_PS1' : [.5, .4, .5,],
                           'e_i_PS1' : [.1, np.nan, .1],
                           'e_z_PS1' : [.1, np.nan, .1],
                           'e_y_PS1' : [.1, np.nan, .1],
                           'Qual_PS1' : [60, np.nan, 20],})

    resultdf =  pd.DataFrame({'gFlags_PS1': {0: 8.0, 1: np.nan, 2: 40.0},
                                 'rFlags_PS1': {0: 0, 1: 320, 2: 40},
                                 'iFlags_PS1': {0: 132, 1: 16, 2: 40},
                                 'zFlags_PS1': {0: 16, 1: 4, 2: 40},
                                 'yFlags_PS1': {0: 8, 1: 4, 2: 40},
                                 'g_PS1': {0: 12., 1: np.nan, 2: np.nan},
                                 'r_PS1': {0: 12.5, 1: np.nan, 2: 23.5},
                                 'i_PS1': {0: 12.1, 1: np.nan, 2: 25.1},
                                 'z_PS1': {0: 12.1, 1: np.nan, 2: 25.1},
                                 'y_PS1': {0: 12.1, 1: np.nan, 2: 25.1},
                                 'e_g_PS1': {0: 1., 1: np.nan, 2: np.nan},
                                 'e_r_PS1': {0: .5, 1: np.nan, 2: 0.5},
                                 'e_i_PS1': {0: .1, 1: np.nan, 2: 0.1},
                                 'e_z_PS1': {0: .1, 1: np.nan, 2: 0.1},
                                 'e_y_PS1': {0: .1, 1: np.nan, 2: 0.1},
                                'Qual_PS1' : {0: 60, 1: np.nan, 2:20}})
    assert resultdf.equals(remove_bad_PS1_photometry(testdf))


def test_remove_bad_Gaia_photometry():
    testdf = pd.DataFrame({'Gflux_Gaia' : [80000.,56022.,430.],
                           'BPflux_Gaia' : [np.nan,56022.,430.],
                           'RPflux_Gaia' : [80000.,56022.,430.],
                           'G_Gaia' : [9., 15., np.nan],
                           'BP_Gaia' : [9., 15., 18.],
                           'RP_Gaia' : [np.nan, 15., 18.],
                           'e_Gflux_Gaia' : [np.nan, 86., 50.],
                           'e_BPflux_Gaia' : [26., 86., 50.],
                           'e_RPflux_Gaia' : [26., 86., 50.],})
    resultdf = pd.DataFrame({'Gflux_Gaia': {0: 80000.0, 1: 56022.0, 2: 430.0},
                             'BPflux_Gaia': {0: np.nan, 1: 56022.0, 2: 430.0},
                             'RPflux_Gaia': {0: 80000.0, 1: 56022.0, 2: 430.0},
                             'G_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'BP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'RP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'e_Gflux_Gaia': {0: np.nan, 1: 86.0, 2: 50.0},
                             'e_BPflux_Gaia': {0: 26.0, 1: 86.0, 2: 50.0},
                             'e_RPflux_Gaia': {0: 26.0, 1: 86.0, 2: 50.0}})
    #hacky solution to assert equality
    assert(resultdf-remove_bad_Gaia_photometry(testdf)).sum().sum() == 0.0

def test_calculate_e_mag_for_Gaia():
    testdf = pd.DataFrame({'Gflux_Gaia': {0: 1.0, 1: 56022.0, 2: .10},
                             'BPflux_Gaia': {0: np.nan, 1: 56022.0, 2: .10},
                             'RPflux_Gaia': {0: 1.0, 1: 56022.0, 2: .10},
                             'G_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'BP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'RP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                             'e_Gflux_Gaia': {0: 99., 1: 0.0, 2: 99.9},
                             'e_BPflux_Gaia': {0: 99.0, 1: 0.0, 2: 99.9},
                             'e_RPflux_Gaia': {0: 99.0, 1: 0.0, 2: 99.9}})
    resultdf = pd.DataFrame( {'Gflux_Gaia': {0: 1.0, 1: 56022.0, 2: 0.1},
                                 'BPflux_Gaia': {0: np.nan, 1: 56022.0, 2: 0.1},
                                 'RPflux_Gaia': {0: 1.0, 1: 56022.0, 2: 0.1},
                                 'G_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                                 'BP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                                 'RP_Gaia': {0: np.nan, 1: 15.0, 2: np.nan},
                                 'e_Gflux_Gaia': {0: 99.0, 1: 0.0, 2: 99.9},
                                 'e_BPflux_Gaia': {0: 99.0, 1: 0.0, 2: 99.9},
                                 'e_RPflux_Gaia': {0: 99.0, 1: 0.0, 2: 99.9},
                                 'e_G_Gaia': {0: 5.0, 1: 0.0, 2: 7.5},
                                 'e_BP_Gaia': {0: np.nan, 1: 0.0, 2: 7.5},
                                 'e_RP_Gaia': {0: 5.0, 1: 0.0, 2: 7.5}})
    #hacky solution to assert equality
    assert(resultdf-calculate_e_mag_for_Gaia(testdf)).sum().sum() == 0.0

def test_convert_PS1_to_Sloan():
    """Serves as a test for :func:`PS_to_SDSS`, too."""
    # -----------------------------------------------------------
    # Prepare fake data:
    n = 4
    df = pd.DataFrame({"g_PS1":np.linspace(-1,2,n),
                      "e_g_PS1":np.linspace(.1,.2,n),
                      "i_PS1":np.linspace(-1,2,n),
                      "e_i_PS1":np.linspace(.1,.2,n),
                      "r_PS1":[np.nan, np.nan, np.nan, np.nan],
                      "e_r_PS1":[np.nan, np.nan, np.nan, np.nan]})

    # -----------------------------------------------------------
    # Call function:
    # Should throw two warnings because z and y bands are not in df
    with pytest.warns(UserWarning) as record:
            dff = convert_PS1_to_Sloan(df)
            assert len(record) == 2

    # Check results from PS_to_SDSS
    assert ((dff.e_i_SDSS - dff.e_i_PS1).dropna().values > 0 ).all() # uncertainties can only grow
    assert df.shape[1] == 12
    assert dff.r_SDSS.isnull().values.all()
    assert dff.e_r_SDSS.isnull().values.all()



def test_correct_BP_RP_for_extinction():
    # -------------------------------------------------------------
    # Some fake data:

    inputdf = pd.DataFrame({"BPRP_Gaia":[1,2,3.4,5],
                       "e_BP_Gaia":[.1]*4,
                       "e_RP_Gaia":[.1]*4,
                       "Ext_BPRP_Gaia":[np.nan, .3, .3, .3],
                       "eup_Ext_BPRP_Gaia":[np.nan, np.nan, .3, .3],
                       "elo_Ext_BPRP_Gaia":[np.nan, np.nan, .3, .3],})

    s = inputdf.shape

    # -------------------------------------------------------------
    # Call function:

    df = correct_BP_RP_for_extinction(inputdf)

    # -------------------------------------------------------------
    # Do some checks:

    assert df.loc[0,"BPRP_Gaia_corr"] == .7
    assert (df.e_BPRP_Gaia_corr.values == df.e_BPRP_Gaia_corr.mean()).all() # all uncertainties are the same
    assert df.shape[1]-2 == s[1]
    assert df.shape[0] == s[0]