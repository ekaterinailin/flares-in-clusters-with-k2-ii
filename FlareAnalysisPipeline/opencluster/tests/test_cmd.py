import pandas as pd
import numpy as np
import pytest

from collections import OrderedDict
from astropy import units as u
from gaia.gaia import calculate_distance_from_parallax

from ..opencluster import OpenCluster
from ..cmd import (K_2MASS_to_KBB,
                   color_2MASS_to_BB,
                   color_BB_to_Johnson,
                   Teff_Boyajian,
                   Teff_Apsis,
                   mann_formula,
                   Teff_Mann,
                   prioritize_Teff_Mann,
                   cut_Teff_Mann_on_Teff,
                   cut_Teff_Mann_on_Mks,
                   correct_for_extinction,
                   query_dustmap_w_percentiles,
                   calculate_weighted_mean)


def test_K_2MASS_to_KBB():
    K_2MASS = pd.Series([10., 10., 10])
    J_K_BB = pd.Series([0.,0.,np.nan])
    K_BB, e_K_BB = K_2MASS_to_KBB(K_2MASS, J_K_BB, e_K_2MASS=0., e_J_K_BB=0.)
    res_K_BB = pd.Series([10.039,10.039,np.nan])
    res_e_K_BB = pd.Series([0.007**2,0.007**2,np.nan])
    assert K_BB.equals(res_K_BB)
    assert e_K_BB.equals(res_e_K_BB)

def test_color_2MASS_to_BB():
    testdf = pd.DataFrame({'J_2MASS' : [10., 13., np.nan],
                       'H_2MASS' : [10., 13., np.nan],
                       'K_2MASS' : [10., 13., np.nan],
                       'g_SDSS'  : [10., 13., np.nan],
                       'r_SDSS' : [10., 13., np.nan],
                       'i_SDSS' : [10., 13., np.nan],
                       'z_SDSS' : [10., 13., np.nan],
                       'y_SDSS' : [10., 13., np.nan],
                       'e_J_2MASS' : [.10, .13, np.nan],
                       'e_H_2MASS' : [.10, .13, np.nan],
                       'e_K_2MASS' : [.10, .13, np.nan],
                       'e_g_SDSS'  : [.10, .13, np.nan],
                       'e_r_SDSS' : [.10, .13, np.nan],
                       'e_i_SDSS' : [.10, .13, np.nan],
                       'e_z_SDSS' : [.10, .13, np.nan],
                       'e_y_SDSS' : [.10, .13, np.nan],})
    resultdf = pd.DataFrame({'J_2MASS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'H_2MASS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'K_2MASS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'g_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'r_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'i_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'z_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'y_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'J_H_BB': {0: 0.049494949494949494, 1: 0.049494949494949494, 2: np.nan},
                             'J_K_BB': {0: 0.01831129196337742, 1: 0.01831129196337742, 2: np.nan},
                             'H_K_BB': {0: -0.035015447991761074, 1: -0.035015447991761074, 2: np.nan},
                             'K_BB': {0: 10.038981688708036, 1: 13.038981688708036, 2: np.nan},
                             'e_K_BB': {0: 4.9008382585334204e-05, 1: 4.9008382585334204e-05, 2: np.nan}})
    assert (resultdf-color_2MASS_to_BB(testdf)).sum().sum() == 0.

def test_color_BB_to_Johnson():
    testdf = pd.DataFrame({'J_K_BB' : [0., 1., 1.],
                           'J_H_BB' : [0., 1., np.nan],
                           'H_K_BB' : [0., 1., np.nan],
                           'e_J_K_BB' : [0.01, .1, .0],
                           'e_J_H_BB' : [0.01, .1, np.nan],
                           'e_H_K_BB' : [0.01, .1, np.nan],
                           'K_BB'  : [0., 13., 13.],
                           'e_K_BB' : [0., .2, np.nan],})
    resultdf = pd.DataFrame({'J_K_BB': {0: 0.0, 1: 1.0, 2: 1.0},
                             'J_H_BB': {0: 0.0, 1: 1.0, 2:np.nan},
                             'H_K_BB': {0: 0.0, 1: 1.0, 2:np.nan},
                             'K_BB': {0: 0.0, 1: 13.0, 2: 13.0},
                             'e_K_BB': {0: 0.0, 1: 0.2, 2:np.nan},
                             'J_H_Johnson': {0: 0.01, 1: 1.02, 2:np.nan},
                             'K_Johnson': {0: 0.0, 1: 13.0, 2: 13.0},
                             'e_K_Johnson': {0: 0.0, 1: 0.2, 2:np.nan},
                             'H_K_Johnson': {0: 0.01, 1: 0.92, 2:np.nan},
                             'J_K_Johnson': {0: 0.01, 1: 1.0, 2: 1.0},
                             'J_Johnson': {0: 0.01, 1: 14.0, 2: 14.0},
                             'H_Johnson': {0: 0.01, 1: 13.92, 2:np.nan}})
    assert (resultdf-color_BB_to_Johnson(testdf)).sum().sum() == 0.

def test_Teff_Boyajian():
    testdf = pd.DataFrame({'J_2MASS' : [10., 13., np.nan],
                           'H_2MASS' : [10., 13., np.nan],
                           'K_2MASS' : [np.nan, 13., np.nan],
                           'g_SDSS'  : [10., 13., np.nan],
                           'r_SDSS' : [10., 13., np.nan],
                           'i_SDSS' : [10., 13., np.nan],
                           'z_SDSS' : [10., 13., np.nan],
                           'y_SDSS' : [10., 13., np.nan],
                           'e_r_SDSS' : [.1,.1,.1,],
                           'e_g_SDSS' : [.1,.1,.1,],
                           'e_i_SDSS' : [.1,.1,.1,],
                           'e_y_SDSS' : [.1,.1,.1,],
                           'e_z_SDSS' : [.1,.1,.1,],
                           'e_J_2MASS' : [.1,.1,.1,],
                           'e_K_2MASS' : [.1,.1,.1,],
                           'e_H_2MASS' : [.1,.1, np.nan,],
                           'e_BP_Gaia' : [.1,.1,.1],
                           'e_RP_Gaia' : [.1,.1,.1],
                           'e_FeH' : [.01, np.nan, 0.],
                           'FeH' : [.0,.5,0.]})
    resultdf = pd.DataFrame({'J_2MASS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'H_2MASS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'K_2MASS': {0: np.nan, 1: 13.0, 2: np.nan},
                             'g_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'r_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'i_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'z_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'y_SDSS': {0: 10.0, 1: 13.0, 2: np.nan},
                             'J_H_BB': {0: 0.049494949494949494, 1: 0.049494949494949494, 2: np.nan},
                             'J_K_BB': {0: np.nan, 1: 0.01831129196337742, 2: np.nan},
                             'H_K_BB': {0: np.nan, 1: -0.035015447991761074, 2: np.nan},
                             'K_BB': {0: np.nan, 1: 13.038981688708036, 2: np.nan},
                             'e_K_BB': {0: np.nan, 1: 4.9008382585334204e-05, 2: np.nan},
                             'J_H_Johnson': {0: 0.05998989898989899, 1: 0.05998989898989899, 2: np.nan},
                             'K_Johnson': {0: np.nan, 1: 13.038981688708036, 2: np.nan},
                             'e_K_Johnson': {0: np.nan, 1: 4.9008382585334204e-05, 2: np.nan},
                             'H_K_Johnson': {0: np.nan, 1: -0.021864057672502574, 2: np.nan},
                             'J_K_Johnson': {0: np.nan, 1: 0.028128179043743644, 2: np.nan},
                             'J_Johnson': {0: np.nan, 1: 13.06710986775178, 2: np.nan},
                             'H_Johnson': {0: np.nan, 1: 13.017117631035534, 2: np.nan},
                             'Teff_g_z_1': {0: 7089.0, 1: 7089.0, 2: np.nan},
                             'Teff_g_z_2': {0: 7131.0, 1: 7131.0, 2: np.nan},
                             'Teff_g_i_1': {0: 7279.0, 1: 7279.0, 2: np.nan},
                             'Teff_g_i_2': {0: 7325.0, 1: 7325.0, 2: np.nan},
                             'Teff_g_J_1': {0: np.nan, 1: 8760.34908681167, 2: np.nan},
                             'Teff_g_J_2': {0: np.nan, 1: 8791.685229821445, 2: np.nan},
                             'Teff_g_J_3': {0: np.nan, 1: 8958.654634234164, 2: np.nan},
                             'Teff_g_J_4': {0: np.nan, 1: 8889.369291487143, 2: np.nan},
                             'Teff_g_H_1': {0: np.nan, 1: 8627.266682560836, 2: np.nan},
                             'Teff_g_H_2': {0: np.nan, 1: 8815.535760296852, 2: np.nan},
                             'Teff_g_H_3': {0: np.nan, 1: 8785.140587710297, 2: np.nan},
                             'Teff_g_H_4': {0: np.nan, 1: 8786.157998856195, 2: np.nan},
                             'Teff_g_K_1': {0: np.nan, 1: 8607.75131523672, 2: np.nan},
                             'Teff_g_K_2': {0: np.nan, 1: 8600.514380036813, 2: np.nan},
                             'Teff_g_K_3': {0: np.nan, 1: 8703.458290081515, 2: np.nan},
                             'Teff_g_K_4': {0: np.nan, 1: 8565.254767956882, 2: np.nan}})
    testres = Teff_Boyajian(testdf)

    # do a coarse consistency check:

    assert (resultdf - testres).sum().sum() == 0.

def test_Teff_Apsis():
    # create fake data
    teffs = np.linspace(3000,10000,20)
    bprp = np.linspace(.8,2.5,20)
    teffs[[3,7]] = np.nan
    
    cols = ["Teff_Apsis", "e_Teff_Apsis"]
    df = pd.DataFrame({"teff_val_Gaia":teffs,
                       "a":np.random.rand(20),
                       "BPRP_Gaia":bprp})

    # call function
    res = Teff_Apsis(df)

    # do some checks
    assert res[cols].dropna(how="all").shape[0] == 6
    assert res.columns.values.shape[0] == 5
    assert (res.Teff_Apsis.dropna().values > 4000.).all()
    assert (res.Teff_Apsis.dropna().values < 6750.).all()
    assert (res.e_Teff_Apsis.dropna().values == 175.).all()
    assert (res.a.values == df.a.values).all()

def test_mann_formula():
    p1 = pd.Series({'a':1.,'b':1.,'c':1.,'d':1.,'e':1.,'JH_FeH':1.,
                    'JH2':1.,'formula':4, 'Teff_factor': 3500.,
                    'sigma':3,'add_in_quadrature_K':4})
    p2 = pd.Series({'a':1.,'b':1.,'c':1.,'d':1.,'e':1.,'JH_FeH':1.,
                    'JH2':1.,'formula':6, 'Teff_factor': 3500.,
                    'sigma':3,'add_in_quadrature_K':4})
    p3 = pd.Series({'a':1.,'b':1.,'c':1.,'d':1.,'e':1.,'JH_FeH':1.,
                    'JH2':1.,'formula':7, 'Teff_factor': 3500.,
                    'sigma':3,'add_in_quadrature_K':4})
    color = pd.Series([1., np.nan, 0.])
    extra = pd.Series([1., 1., 0.])
    err = pd.Series([0., 0., 0.])
    extra_err = pd.Series([0., 0., 0.])
    x = [color, extra, err, extra_err]
    r1, r1err = np.array([17500.0, np.nan, 3500.0]), np.array([5., 500.024999 , 5.])
    r2, r2err = np.array([21000.0, np.nan, 3500.0]), np.array([5., 500.024999 , 5.])
    r3, r3err = np.array([24500.0, np.nan, 3500.0]), np.array([5., 500.024999 , 5.])
    assert np.allclose(mann_formula(x,p1)[0], r1, equal_nan=True)
    assert np.allclose(mann_formula(x,p2)[0], r2, equal_nan=True)
    assert np.allclose(mann_formula(x,p3)[0], r3, equal_nan=True)

    assert np.allclose(mann_formula(x,p1)[1], r1err, equal_nan=True)
    assert np.allclose(mann_formula(x,p2)[1], r2err, equal_nan=True)
    assert np.allclose(mann_formula(x,p3)[1], r3err, equal_nan=True)

def test_Teff_Mann():
    testdf = pd.DataFrame({'r_SDSS' : [14., 16., 10., 10.],
                           'z_SDSS' : [12, 14, 10., 10.],
                           'J_2MASS' : [12., 14., 10., 10.],
                           'H_2MASS' : [10., 12., 10., np.nan,],
                           'BP_Gaia' : [13., 12., 10., 10.,],
                           'RP_Gaia' : [10., 9., 10., 10.,],
                           'FeH' : [.1, np.nan, 0., 0.],
                           'e_r_SDSS' : [.1,.1,.1,.1],
                           'e_z_SDSS' : [.1,.1,.1,.1],
                           'e_J_2MASS' : [.1,.1,.1,.1],
                           'e_H_2MASS' : [.1,.1,.1, np.nan,],
                           'e_BP_Gaia' : [.1,.1,.1,.1],
                           'e_RP_Gaia' : [.1,.1,.1,.1],
                           'e_FeH' : [.01, np.nan, 0.001, 0.],
                           'BPRP_Gaia_corr' : [3.0,3.0, .0,.0],
                           'e_BPRP_Gaia_corr' : [.2,.1, .1,.1]})

    resultdf = pd.DataFrame({'r_SDSS': {0: 14.0, 1: 16.0, 2: 10.0, 3: 10.0},
                             'z_SDSS': {0: 12.0, 1: 14.0, 2: 10.0, 3: 10.0},
                             'J_2MASS': {0: 12.0, 1: 14.0, 2: 10.0, 3: 10.0},
                             'H_2MASS': {0: 10.0, 1: 12.0, 2: 10.0, 3: np.nan},
                             'BP_Gaia': {0: 13.0, 1: 12.0, 2: 10.0, 3: 10.0},
                             'RP_Gaia': {0: 10.0, 1: 9.0, 2: 10.0, 3: 10.0},
                             'FeH': {0: 0.1, 1: np.nan, 2: 0.0, 3: 0.0},
                             'e_r_SDSS' : {0:.1,1:.1,2:.1,3:.1},
                           'e_z_SDSS' : {0:.1,1:.1,2:.1,3:.1},
                           'e_J_2MASS' : {0:.1,1:.1,2:.1,3:.1},
                           'e_H_2MASS' : {0:.1,1:.1,2:.1, 3:np.nan,},
                           'e_BP_Gaia' : {0:.1,1:.1,2:.1,3:.1},
                           'e_RP_Gaia' : {0:.1,1:.1,2:.1,3:.1},
                           'e_FeH' : {0:.01, 1:np.nan,2: 0.001, 3:0.},
                             'BPRP_Gaia_corr' : {0:3.0,1:3.0, 2:.0,3:.0},
                             'e_BPRP_Gaia_corr' : {0:.2,1:.1,2:.1,3:.1},
                             'Teff_Mann_BP_RP_nan': {0: 3262.3149999999973,
                              1: 3262.3149999999973, 2: 11357.5, 3: 11357.5},
                             'Teff_Mann_r_z_nan': {0: 3359.159999999999,1: 3359.159999999999,
                              2: 5414.5,3: 5414.5},
                             'Teff_Mann_r_J_nan': {0: 4065.0119999999993, 1: 4065.0119999999993,
                              2: 8557.5, 3: 8557.5},
                             'Teff_Mann_BP_RP_FeH': {0: 3281.0644999999963, 1: np.nan, 2: 9922.5, 3: 9922.5},
                             'Teff_Mann_r_z_FeH': {0: 3376.3660000000004, 1: np.nan, 2: 5501.999999999999,
                              3: 5501.999999999999},
                             'Teff_Mann_r_J_FeH': {0: 4117.669500000002, 1: np.nan, 2: 8862.0, 3: 8862.0},
                             'Teff_Mann_BP_RP_isJH': {0: 4265.800000000004, 1: 4265.800000000004,
                              2: 11102.0, 3: np.nan},
                             'Teff_Mann_r_z_isJH': {0: 4023.0399999999995,1: 4023.0399999999995,
                              2: 4844.0,3: np.nan},
                             'Teff_Mann_r_J_isJH': {0: 5240.479999999998,1: 5240.479999999998,
                              2: 7528.499999999999,3: np.nan}})

    assert (Teff_Mann(testdf)-resultdf).sum().sum() == 0.

def test_prioritize_Teff_Mann():
    testdf = pd.DataFrame({'Teff_Mann_BP_RP_nan' : [3000., 3000., 3000.,],
                           'Teff_Mann_r_z_nan' : [3000., np.nan, 3000.,],
                           'Teff_Mann_r_J_nan' : [3000., 3000., 3000.,],
                           'Teff_Mann_BP_RP_FeH' : [3000., np.nan, np.nan,],
                           'Teff_Mann_r_z_FeH' : [3000., np.nan, 3000.,],
                           'Teff_Mann_r_J_FeH' : [3000., 3000., 3000.,],
                           'Teff_Mann_BP_RP_isJH' : [3000., np.nan, 3000.,],
                           'Teff_Mann_r_z_isJH' : [3000., 3000., 3000.,],
                           'Teff_Mann_r_J_isJH' : [3000., 3000., 3000.,],})
    resultdf = pd.DataFrame( {'Teff_Mann_BP_RP_nan': {0: np.nan, 1: 3000.0, 2: np.nan},
                              'Teff_Mann_r_z_nan': {0: np.nan, 1: np.nan, 2: np.nan},
                              'Teff_Mann_r_J_nan': {0: np.nan, 1: np.nan, 2: np.nan},
                              'Teff_Mann_BP_RP_FeH': {0: 3000.0, 1: np.nan, 2: np.nan},
                              'Teff_Mann_r_z_FeH': {0: 3000.0, 1: np.nan, 2: 3000.0},
                              'Teff_Mann_r_J_FeH': {0: 3000.0, 1: 3000.0, 2: 3000.0},
                              'Teff_Mann_BP_RP_isJH': {0: np.nan, 1: np.nan, 2: 3000.0},
                              'Teff_Mann_r_z_isJH': {0: np.nan, 1: 3000.0, 2: np.nan},
                              'Teff_Mann_r_J_isJH': {0: np.nan, 1: np.nan, 2: np.nan}})
    assert prioritize_Teff_Mann(testdf).equals(resultdf)

def test_cut_Teff_Mann_on_Teff():
    testdf = pd.DataFrame({'Teff_Mann_BP_RP_nan' : [3000., 2500., 5000.,],
                           'Teff_Mann_BP_RP_FeH' : [3000., 2500., 5000.,],
                           'Teff_Mann_BP_RP_isJH' : [3000., 2500., 5000.,],
                           'Taff_Mann_BP_RP_isJH' : [3000., 2500., 5000.,],
                           'e_Teff_Mann_BP_RP_nan' : [300., 250., 500.,],
                           'e_Teff_Mann_BP_RP_FeH' : [300., 250., 500.,],
                           'e_Teff_Mann_BP_RP_isJH' : [300., 250., 500.,],
                           'e_Taff_Mann_BP_RP_isJH' : [300., 250., 500.,]})
    resultdf = pd.DataFrame({'Teff_Mann_BP_RP_nan': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Teff_Mann_BP_RP_FeH': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Teff_Mann_BP_RP_isJH': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Taff_Mann_BP_RP_isJH': {0: 3000.0, 1: 2500.0, 2: 5000.0},
                             'e_Teff_Mann_BP_RP_nan': {0: 300.0, 1: np.nan, 2: np.nan},
                             'e_Teff_Mann_BP_RP_FeH': {0: 300.0, 1: np.nan, 2: np.nan},
                             'e_Teff_Mann_BP_RP_isJH': {0: 300.0, 1: np.nan, 2: np.nan},
                             'e_Taff_Mann_BP_RP_isJH': {0: 300.0, 1: 250.0, 2: 500.0}})
    print(cut_Teff_Mann_on_Teff(testdf))
    assert cut_Teff_Mann_on_Teff(testdf).equals(resultdf)

def test_cut_Teff_Mann_on_Mks():
    testdf = pd.DataFrame({'Teff_Mann_BP_RP_nan' : [3000., 2500., 5000.,],
                           'Teff_Mann_BP_RP_FeH' : [3000., 2500., 5000.,],
                           'Teff_Mann_BP_RP_isJH' : [3000., 2500., 5000.,],
                           'Taff_Mann_BP_RP_isJH' : [3000., 2500., 5000.,],
                           'e_Teff_Mann_BP_RP_nan' : [300., 250., 500.,],
                           'e_Teff_Mann_BP_RP_FeH' : [300., 250., 500.,],
                           'e_Teff_Mann_BP_RP_isJH' : [300., 250., 500.,],
                           'e_Taff_Mann_BP_RP_isJH' : [300., 250., 500.,],
                           'K_2MASS':[15, 15, 15],
                           'parallax_Gaia':[3,np.nan,0.5],
                           'parallax_error_Gaia':[.03,np.nan,np.nan]})
    resultdf = pd.DataFrame({'Teff_Mann_BP_RP_nan': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Teff_Mann_BP_RP_FeH': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Teff_Mann_BP_RP_isJH': {0: 3000.0, 1: np.nan, 2: np.nan},
                             'Taff_Mann_BP_RP_isJH': {0: 3000.0, 1: 2500.0, 2: 5000.0},
                             'e_Teff_Mann_BP_RP_nan' : {0: 300.0, 1: np.nan, 2: np.nan},
                               'e_Teff_Mann_BP_RP_FeH' :{0: 300.0, 1: np.nan, 2: np.nan},
                               'e_Teff_Mann_BP_RP_isJH' :{0: 300.0, 1: np.nan, 2: np.nan},
                               'e_Taff_Mann_BP_RP_isJH' : {0: 300.0, 1: 250.0, 2: 500.0},
                             'K_2MASS':{0:15, 1:15, 2:15},
                             'parallax_Gaia':{0:3,1:np.nan,2:0.5},
                             'parallax_error_Gaia':{0:.03,1:np.nan,2:np.nan}})

    assert cut_Teff_Mann_on_Mks(testdf).equals(resultdf)


def test_correct_for_extinction():

    # read in test data:
    testdf = pd.read_csv('opencluster/tests/testdf.csv')

    #----------------------------------------------------------------
    # First check if no extincion occurs

    # Init open cluster:

    OC = OpenCluster(cluster='hyades',h_cluster='Hyades',age=600,
                              age_unit=1e6*u.yr, stars=testdf, feh=0.13, u_feh=0.01)

    # Define DataFrames to compare: They should be equal.
    before = OC.stars
    after = correct_for_extinction(before, thresh=14)

    # Test the flags are set to 10 == no reliable extinction values
    assert (after.bayestar2019_flags.values == np.array([10]*10)).all()

    # Test if they are equal for different bands and uncertainties:

    cols = ['J_2MASS', 'g_PS1', 'r_PS1', 'e_J_2MASS', 'e_g_PS1', 'e_r_PS1']
    for col in cols:
        a, b = before[col], after[col]
        assert ((a == b) | (np.isnan(a) & np.isnan(b)) | np.isnan(b)).all()

    #-------------------------------------------------------------------
    # Second check if extinction correction brightens the stars and
    # increases uncertainty.

    # Init open cluster with increased distance (=more dust!)

    OC = OpenCluster(cluster='hyades',h_cluster='Hyades',age=600,
                              age_unit=1e6*u.yr, stars=testdf, feh=0.13, u_feh=0.01)
    OC.stars['parallax_Gaia'] /= 100.

    # Define DataFrames to compare:

    before = OC.stars
    after = correct_for_extinction(before, thresh=14)

    # Test the flags are set:
    assert (after.bayestar2019_flags.values == np.array([10,10,2,2,2,2,10,2,2,2])).all()

    # Test if uncertainties behave as they should:
    cols = ['e_J_2MASS', 'e_g_PS1', 'e_r_PS1']
    for col in cols:
        a, b = before[col], after[col]
        # uncertainties must grow or stay the same
        assert ((a <= b) | (np.isnan(a) & np.isnan(b)) | np.isnan(b)).all()

    # Test if band magnitudes behave as they should:
    cols = ['J_2MASS', 'g_PS1', 'r_PS1']
    for col in cols:
        a, b = before[col], after[col]
        # band magnitude must decrease or stay the same
        assert ((a >= b) | (np.isnan(a) & np.isnan(b)) | np.isnan(b)).all()


def test_query_dustmaps_w_percentiles():

    # Fetch test data:
    testdf = pd.read_csv('opencluster/tests/testdf.csv')
    OC = OpenCluster(cluster='hyades',h_cluster='Hyades',age=600,
                     age_unit=1e6*u.yr, stars=testdf, feh=0.13, u_feh=0.01)

    # extract distance
    df = OC.stars.loc[:,OC.stars.columns.str.contains("Gaia")].rename(columns = lambda x : str(x)[:-5])
    df = calculate_distance_from_parallax(df, check_GoF=False)

    # define inputs
    d = df.distance.values * u.pc
    ra = df.RAJ2000_.values * u.deg
    dec = df.DEJ2000_.values * u.deg

    # test function
    ext, ext_err, flags = query_dustmap_w_percentiles(d, ra, dec)

    # do some checks
    assert (flags == np.array([10, 10, 10, 10, 10, 10, 30, 10, 10, 10])).all()
    assert np.allclose(ext, np.array([ 0.,  0.,  0.,  0.,  0.,  0., np.nan,  0.,  0.,  0.]), equal_nan=True)
    assert np.allclose(ext_err, np.array([ 0.,  0.,  0.,  0.,  0.,  0., np.nan,  0.,  0.,  0.]), equal_nan=True)


def test_calculate_weighted_mean():
    #--------------------------------------------------------------
    # Mock some data
    n = 5
    testdf = pd.DataFrame({"Teff_a":np.linspace(3000,5000,n),
                           "e_Teff_a":[100]*n,
                           "Teff_b":np.linspace(3000-100,5000-100,n),
                           "e_Teff_b":[100]*n,
                           "Teff_c":np.linspace(3000+100,5000+100,n),
                           "e_Teff_c":[100]*n,})
    testdf.loc[0,"e_Teff_c"] = np.nan
    testdf.loc[1,["Teff_c", "e_Teff_c"]] = np.nan
    testdf.loc[2,["Teff_c", "e_Teff_c", "Teff_b", "e_Teff_b"]] = np.nan
    testdf.loc[3,["Teff_c", "e_Teff_c", "Teff_b", "e_Teff_b", "e_Teff_a"]] = np.nan

    #--------------------------------------------------------------
    # call function
    res = calculate_weighted_mean(testdf)

    #--------------------------------------------------------------
    # do some checks
    y = np.array([100/np.sqrt(2),100/np.sqrt(2),100.,100/np.sqrt(3)])
    assert res.Teff_std.dropna().values == pytest.approx(y)
    z = np.array([2950.,3450.,4000.,5000.])
    assert res.Teff_median.dropna().values == pytest.approx(z)
