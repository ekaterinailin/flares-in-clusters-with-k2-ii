import warnings
from copy import deepcopy

import pandas as pd
import numpy as np
import astropy.units as u
import pytest

from ..ffd import FFD
from ..utils import generate_random_power_law_distribution

def test__get_threshold_warning():
    data = generate_random_power_law_distribution(10,1e4,-1., 100, seed=2)
    ffd = FFD(ED=data)
    ffd.flares = pd.DataFrame({'Lum_Kepler' : np.linspace(1e28,1e30,100),
                               'ed_rec' : data})

    ffd.cutoff_energy_lower = 1e30 * u.erg
    assert ffd.flares.shape == (93,3)
    assert ffd.flares.en.min() > 1e30

def test_cutoff_energy():
    data = generate_random_power_law_distribution(10,1e4,-1., 100, seed=2)
    ffd = FFD(ED=data)
    ffd.flares = pd.DataFrame({'Lum_Kepler' : np.linspace(1e28,1e30,100),
                               'ed_rec' : data})
    ffd.cutoff_energy_lower = 5.*u.Joule
    with pytest.raises(AssertionError):
        ffd.cutoff_energy_lower = 1e33
    assert ffd.cutoff_energy_lower == u.Quantity(5e7*u.erg)
    assert ffd.flares.shape == (100,3)
    ffd.cutoff_energy_lower = 1e31*u.erg
    assert ffd.flares.en.min() > ffd.cutoff_energy_lower.value
    assert ffd.flares.shape == (47, 3)

def test_cutoff_ED():
    data = generate_random_power_law_distribution(10,1e4,-1., 100, seed=2)
    ffd = FFD(ED=data)
    ffd.cutoff_ED_lower = 5.*u.minute
    assert ffd.flares.shape == (4,2)
    assert ffd.cutoff_ED_lower == u.Quantity(300, 's')
    with pytest.raises(AssertionError):
        ffd.cutoff_ED_lower = 50.
    with pytest.raises(AssertionError):
        ffd.cutoff_ED_upper = 50.

def test_is_powerlaw_by_KS():
    size = 100
    campaign = size * [0]
    EPIC = size * [211119999]
    tobs = size * [8766.]
    minval, maxval = 10., 1e4
    ED = generate_random_power_law_distribution(minval, maxval, -1., size=size, seed=10)
    with pytest.raises(AttributeError):
        ffdt = FFD(ED=ED,
                   alpha=2.,
                   flares = pd.DataFrame({'ed_rec':ED,
                                          'Campaign':campaign,
                                          'EPIC':EPIC,
                                          'total_obs_time' : tobs}))
        ffdt.is_powerlaw_by_KS().KS
    ffdt = FFD(ED=ED,
               alpha=2.,
               flares = pd.DataFrame({'ed_rec':ED,
                                          'Campaign':campaign,
                                          'EPIC':EPIC,
                                          'total_obs_time' : tobs}),
               mode='ED')
    assert ffdt.is_powerlaw_by_KS().KS

def test_is_powerlaw_truncated():
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        size = 100
        campaign = size * [0]
        EPIC = size * [211119999]
        tobs = size * [8766.]
        ED = generate_random_power_law_distribution(10, 1000,-1.,size, seed=100)
        ffd = FFD(alpha=3., ED=ED, flares=pd.DataFrame({'ed_rec':ED,
                                                       'Campaign':campaign,
                                                       'EPIC':EPIC,
                                                       'total_obs_time' : tobs}))
        assert ~ffd.is_powerlaw_truncated().truncated
        # Verify some things about the warning
        assert len(w) == 1
        assert issubclass(w[-1].category, UserWarning)
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        size = 100
        campaign = size * [0]
        EPIC = size * [211119999]
        tobs = size * [8766.]
        ED = generate_random_power_law_distribution(10, 1000,-1.,200)
        ED = np.sort(ED)[:100]
        ffd = FFD(alpha=3., ED=ED, flares=pd.DataFrame({'ed_rec':ED,
                                                       'Campaign':campaign,
                                                       'EPIC':EPIC,
                                                       'total_obs_time' : tobs}))
        assert ffd.is_powerlaw_truncated().truncated
        # Verify some things
        assert len(w) == 2
        assert issubclass(w[-1].category, UserWarning)


def test_calculate_FR():
    testffd = FFD(flares=pd.DataFrame({'Campaign' : [4,4,4],
                                       'EPIC' : [222222222, 222222222, 222222222],
                                       'istart' : [20,30,40],
                                       'total_obs_time':[8766.0, 8766.0, 8766.0],
                                       'ed_rec' : [10.,10.,10.]}), tot_obs_time=1.)
    assert testffd.tot_obs_time == 1.
    assert testffd.nflares == 3
    testffd = testffd.calculate_FR()
    assert testffd.FR == 3.

def test_linefit():
    n=50
    data = generate_random_power_law_distribution(10, 100, -1, size=n, seed=10)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n * [4],
                                           'EPIC' : n * [222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':n * [8766.0],#hours
                                           'ed_rec': data,
                                           'Lum_Kepler': [1.]*n}), tot_obs_time=1.)
    testffd.mode='ED'
    testffd = testffd.linefit()

    assert testffd.cutoff_ED_lower == 0*u.s
    assert testffd.cutoff == 0*u.s
    assert testffd.mode == 'ED'
    assert testffd.alpha == pytest.approx(2.752742549358)
    assert testffd.alpha_err == pytest.approx(0.1517230969458871)
    assert testffd.beta == pytest.approx(3879.522717058634)
    assert testffd.beta_err == pytest.approx(0.18993955034827)

def test_powerlawfit():
    n = 30
    data = generate_random_power_law_distribution(10,100,-1,seed=12,size=n)*1e4
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                       'EPIC' : n*[222222222],
                                       'istart' : np.linspace(0,100,n),
                                       'total_obs_time':[0.]*(n-1) + [8766.0],
                                       'ed_rec' : data,
                                       'Lum_Kepler' : [1.] * n}), tot_obs_time=1.)

    testffd = testffd.powerlawfit()
    assert testffd.alpha == pytest.approx(1.897839790,rel=1e-4)
    assert testffd.alpha_err == pytest.approx(0.3812958738116, rel=1e-3)

    rr = generate_random_power_law_distribution(10,100,-1,seed=12,size=n)
    testffd.flares['ed_rec'] = np.rint(rr*1e4)
    testffd.flares.loc[2,'ed_rec'] = np.nan
    testffd = testffd.powerlawfit()
    testffd.flares.loc[2,'ed_rec'] = -5.
    testffd = testffd.powerlawfit()
    assert (testffd.flares.ed_rec >= 0.).all

def test_powerlaw_Maschberger_and_Kroupa():
    n = 30000
    #test1 alpha=2
    data = generate_random_power_law_distribution(10,100,-1,seed=12,size=n)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                           'EPIC' : n*[222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':[0.]*(n-1) + [8766.0],
                                           'ed_rec' : data,
                                           'Lum_Kepler' : [1] * n}), tot_obs_time=1.)
    testffd = testffd.powerlaw_Maschberger_and_Kroupa( alims=[1.01,3.])
    assert testffd.alpha == pytest.approx(2.,rel=5e-3)

    #test2 alpha=3
    data = generate_random_power_law_distribution(10,100,-2,seed=12,size=n)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                           'EPIC' : n*[222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':[0.]*(n-1) + [8766.0],
                                           'ed_rec' : data,
                                           'Lum_Kepler' : [1] * n}), tot_obs_time=1.)
    testffd = testffd.powerlaw_Maschberger_and_Kroupa( alims=[1.01,3.])
    assert testffd.alpha == pytest.approx(3.,rel=5e-3)

def test_calculate_R1s():
    testffd = FFD(alpha=2., alpha_err=0., beta=1., beta_err=0.)
    testffd_res = testffd.calculate_R1s()
    assert testffd_res.R1s == 1.
    assert testffd_res.R1s_err == 0.
    testffd_res = testffd.calculate_R1s(sqr=True)
    assert testffd_res.R1s == 1.
    assert testffd_res.R1s_err == 0.
    testffd = FFD(alpha=2., alpha_err=1., beta=1., beta_err=0.)
    testffd_res = testffd.calculate_R1s()
    assert testffd_res.R1s == 1.
    assert testffd_res.R1s_err == 1.
    testffd_res = testffd.calculate_R1s(sqr=True)
    assert testffd_res.R1s == 1.
    assert testffd_res.R1s_err == 1.
    testffd = FFD(alpha=2., alpha_err=0., beta=10., beta_err=1.)
    testffd_res = testffd.calculate_R1s()
    assert testffd_res.R1s == 10.
    assert testffd_res.R1s_err == 1.
    testffd_res = testffd.calculate_R1s(sqr=True)
    assert testffd_res.R1s == 10.
    assert testffd_res.R1s_err == 1.

def test_calculate_flaring_fraction():
    n = 40
    data = generate_random_power_law_distribution(10,100,-2,seed=12,size=n)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                           'EPIC' : n*[222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':[0.]*(n-1) + [8766.0],
                                           'ed_rec' : data,
                                           'Lum_Kepler' : [1] * n}),nstars=1)
    testffd = testffd.calculate_flaring_fraction()
    assert testffd.flare_frac == 1.

def test_fit_beta_to_powerlaw():
    n = 40
    data = generate_random_power_law_distribution(10,100,-2,seed=12,size=n)
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                           'EPIC' : n*[222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':[0.]*(n-1) + [8766.0],
                                           'ed_rec' : data,
                                           'Lum_Kepler' : [1] * n}),
                  nstars=1.,
                  tot_obs_time=1.,
                  alpha=3., alpha_err=.1)
    edffd = deepcopy(testffd)
    edffd.mode = 'ED'
    edffd = edffd.fit_beta_to_powerlaw()
    assert edffd.beta == pytest.approx(2.*8128.093667984009)
    assert edffd.beta_err == pytest.approx(1632.305300084492)
    assert edffd.alpha == 3.
    assert edffd.alpha_err == 0.1
    enffd = deepcopy(testffd)
    enffd.mode = 'energy'
    enffd = enffd.fit_beta_to_powerlaw()
    assert enffd.beta == pytest.approx(2.*8128.0936619203785)
    assert enffd.beta_err == pytest.approx(1632.305300084492)
    assert enffd.alpha == 3.
    assert enffd.alpha_err == 0.1
    assert enffd.beta == pytest.approx(edffd.beta, rel=1e-3)
    assert enffd.beta_err == pytest.approx(edffd.beta_err, rel=1e-2)

def test_clean_FFD():
    n = 40
    data = generate_random_power_law_distribution(10,100,-2,seed=12,size=n)
    data[3:5] = -4.
    data[6:8] = np.nan
    testffd = FFD(flares=pd.DataFrame({'Campaign' : n*[4],
                                           'EPIC' : n*[222222222],
                                           'istart' : np.linspace(0,100,n),
                                           'total_obs_time':[0.]*(n-1) + [8766.0],
                                           'ed_rec' : data,
                                           'Lum_Kepler' : [1] * n}),
                  nstars=1.,
                  tot_obs_time=1.,
                  alpha=3., alpha_err=.1)
    assert len(testffd.ED) == 40-4