import pandas as pd
import numpy as np
import pytest
from astropy import units as u
from ..analysis import *
from ..utils import generate_random_power_law_distribution
from ..ffd import FFD

from copy import deepcopy

fstars = pd.read_csv('luminosities/fake_luminosities.csv')
fflares = pd.read_csv('flares/fake_flares.csv')

clusters = CLUSTERS.append({'cluster':'fake',
                        'h_cluster':'Fake Cluster',
                        'age (Myr)': 1000.,
                        'FeH' : 0.05,
                        'dist (pc)': 200,
                        'u_age_high':200,
                        'u_age_low':100,
                        'u_feh':0.001,
                        'Teffmax':7000}, ignore_index=True)
OC = generate_OpenCluster('fake', unit='erg', clusters=clusters)

def test_analyse_with_MK():
    #test regular functionality
    n=50
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=20)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.,
                                             mode='ED')
    bffd = analyse_with_MK(testffd)

    assert bffd.alpha == pytest.approx(1.8779744872436217)
    assert bffd.alpha_err == 0.31048566187207655
    assert bffd.beta == pytest.approx(272.8914645110053)
    assert bffd.beta_err == pytest.approx(84.8561170338658)
    assert bffd.alpha_mode == 'MK'
    assert bffd.truncated == True
    assert bffd.R1s == pytest.approx(310.8193557739258)
    #test exception
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': np.full(n,np.nan),
                                           'Lum_Kepler': [1.]*n}),
                                            mode='ED')
    bffd = analyse_with_MK(testffd)
    assert np.isnan(bffd.alpha_err)
    assert np.isnan(bffd.alpha_err)
    assert np.isnan(bffd.beta)
    assert np.isnan(bffd.beta_err)
    assert np.isnan(bffd.R1s)
    assert bffd.alpha_mode == 'MK'
    assert bffd.truncated == False


def test_analyse_with_alpha2():
    #test regular functionality
    n=50
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=20)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.,
                                             mode='ED')
    bffd = analyse_with_alpha2(testffd)
    assert bffd.alpha == 2.
    assert bffd.alpha_err == 0.
    assert bffd.beta == pytest.approx(494.4630844116211)
    assert bffd.beta_err == pytest.approx(5.572031530461575)
    assert bffd.alpha_mode == 'fix2'
    assert bffd.truncated == False
    assert bffd.R1s == bffd.beta
    #test exception
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': np.full(n,np.nan),
                                           'Lum_Kepler': [1.]*n}),
                                            mode='ED')
    bffd = analyse_with_alpha2(testffd)
    assert bffd.alpha == 2.
    assert bffd.alpha_err == 0.
    assert np.isnan(bffd.beta)
    assert np.isnan(bffd.beta_err)
    assert np.isnan(bffd.R1s)
    assert bffd.alpha_mode == 'fix2'
    assert bffd.truncated == False

def test_analyse_with_Bauke():
    #test regular functionality
    n=50
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=20)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.)
    bffd = analyse_with_Bauke(testffd)
    assert bffd.alpha == pytest.approx(1.7581963894364983)
    assert bffd.alpha_err == pytest.approx(0.3408917963921351)
    assert bffd.beta == pytest.approx(145.7644081573807)
    assert bffd.beta_err == 49.808809123761506
    assert bffd.alpha_mode == 'Bauke'
    assert bffd.truncated == True
    assert bffd.R1s == pytest.approx(192.2515197753906)
    #test exception
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': np.full(n,np.nan),
                                           'Lum_Kepler': [1.]*n}))
    bffd = analyse_with_Bauke(testffd)
    assert np.isnan(bffd.alpha)
    assert np.isnan(bffd.alpha_err)
    assert np.isnan(bffd.beta)
    assert np.isnan(bffd.beta_err)
    assert np.isnan(bffd.R1s)
    assert bffd.alpha_mode == 'Bauke'
    assert bffd.truncated == False

def test_analyse_with_linfit():
    #test regular functionality
    n=20
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=12)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.)
    testffd.cutoff=0.*u.s
    lin_ffd = analyse_with_linfit(testffd)
    assert lin_ffd.alpha == pytest.approx(2.089435345331606)
    assert lin_ffd.alpha_err == pytest.approx(0.11687656149082083)
    assert lin_ffd.beta == pytest.approx(244.19963212386557)
    assert lin_ffd.beta_err == pytest.approx(0.14089739153154884)
    assert lin_ffd.mode == 'ED'
    assert lin_ffd.cutoff == 0.*u.s
    assert lin_ffd.cutoff_ED_lower == lin_ffd.cutoff
    assert lin_ffd.R1s == pytest.approx(224.1524778595606)
    assert lin_ffd.truncated == False
    assert lin_ffd.tot_obs_time == 1.
    #test exception
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': np.full(n,np.nan),
                                           'Lum_Kepler': [1.]*n}))

    testffd.cutoff=0.*u.s
    lin_ffd = analyse_with_linfit(testffd)

    assert np.isnan(lin_ffd.alpha)
    assert np.isnan(lin_ffd.alpha_err)
    assert np.isnan(lin_ffd.beta)
    assert np.isnan(lin_ffd.beta_err)
    assert lin_ffd.alpha_mode == 'linear fit'
    assert np.isnan(lin_ffd.R1s)
    assert lin_ffd.truncated == False

def test_find_cutoff_by_KS():
    n=20
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=12)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.)
    assert find_cutoff_by_KS(testffd, u.s, factor=1.02) == testffd.ED.min()*1.02 * u.s
    n=30
    data = np.empty(n)
    data[:20] = generate_random_power_law_distribution(10, 100, -1, size=20, seed=12)
    data[20:] = np.linspace(0,2,10)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}),
                                            tot_obs_time=1.)
    assert find_cutoff_by_KS(testffd, u.s, factor=1.02) == 2.04 *u.s
    #test exception
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': np.full(n,np.nan),
                                           'Lum_Kepler': [1.]*n}))
    with pytest.raises(ValueError):
        find_cutoff_by_KS(testffd, u.s)

def test_generate_OpenCluster():
    #check that the inputs are still there

    assert OC.cluster == 'fake'
    assert OC.h_cluster == 'Fake Cluster'
    assert OC.age == 1000.
    assert OC.feh == 0.05
    assert OC.u_age_low == 100
    assert OC.u_age_high == 200
    assert OC.u_feh == 0.001
    #print(OC.flares.head(), fflares.head()) This is not true anymore because we modify flares to be the same as flares_on_stars
    #assert (OC.flares.values==fflares.values).all()
    assert 'Lum_Kepler' in OC.flares_on_stars.columns.values
    #check some data points
    assert OC.flares_on_stars.EPIC[0] == fstars.EPIC[0]
   # assert OC.flares_on_stars.e_i_PS1[2] == fstars.e_i_PS1[1]
   # assert OC.flares_on_stars.i_PS1[2] == fstars.i_PS1[1]
   # assert OC.flares_on_stars.i_PS1[2] == fstars.i_PS1[1]
    assert (OC.flares_on_stars.ed_rec == fflares.ed_rec).all()

