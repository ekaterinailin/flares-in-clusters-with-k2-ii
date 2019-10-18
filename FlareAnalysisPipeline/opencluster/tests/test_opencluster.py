import pandas as pd
import numpy as np
import pytest
from astropy import units as u
from ..opencluster import OpenCluster
from ..lum import read_Kepler_response
from ..analysis import *
from .. import PACKAGEDIR
from copy import deepcopy

import specmatchemp.library
Kp = read_Kepler_response()
lib = specmatchemp.library.read_hdf(wavlim=[Kp.wav.min(),Kp.wav.max()])

testdf = pd.read_csv('opencluster/tests/testdf.csv')
TestCluster = OpenCluster(cluster='hyades',h_cluster='Hyades',age=600,
                          age_unit=1e6*u.yr, stars=testdf, feh=0.13, u_feh=0.01,
                          extinction=False)

analysis_cluster = deepcopy(TestCluster)
analysis_cluster.determine_Teff()
analysis_cluster.find_Rstar(lib)
h = analysis_cluster.find_SED(lib)

analysis_cluster.calculate_luminosities()

analysis_cluster.flares = pd.DataFrame({'EPIC':[210942999,210942999,210563410],
                                        'Campaign':[4,4,4,],
                                        'istart' : [20,30,20],
                                        'ed_rec':[40,400,50],
                                        'ampl_rec':[.1,.5,.08,],
                                        'total_obs_time' : [500,500,500]})




def test_init():
    testcluster = deepcopy(analysis_cluster)
    assert testcluster.age_unit == 1e6*u.yr


def test_determine_Teff():
    hyades = deepcopy(TestCluster)
    before = hyades.stars.columns.values
    h = hyades.determine_Teff()
    after = hyades.stars.columns.values
    assert (set(np.setdiff1d(after, before)) ==
            set(['H_Johnson', 'H_K_BB', 'H_K_Johnson', 'J_H_BB', 'J_H_Johnson', 'J_Johnson',
                 'J_K_BB', 'J_K_Johnson', 'K_BB', 'K_Johnson', 'Teff_Mann_BP_RP_FeH',
                 'Teff_Mann_BP_RP_isJH', 'Teff_Mann_BP_RP_nan', 'Teff_Mann_r_J_FeH',
                 'Teff_Mann_r_J_isJH', 'Teff_Mann_r_J_nan', 'Teff_Mann_r_z_FeH',
                 'Teff_Mann_r_z_isJH', 'Teff_Mann_r_z_nan', 'Teff_g_H_Boy', 'Teff_g_J_Boy',
                 'Teff_g_K_Boy', 'Teff_g_i_Boy', 'Teff_g_z_Boy', 'Teff_median',
                 'Teff_std', 'e_H_Johnson', 'e_H_K_BB', 'e_H_K_Johnson', 'e_J_H_BB',
                 'e_J_H_Johnson', 'e_J_Johnson', 'e_J_K_BB', 'e_J_K_Johnson', 'e_K_BB',
                 'e_K_Johnson', 'e_Teff_Mann_BP_RP_FeH', 'e_Teff_Mann_BP_RP_isJH',
                 'e_Teff_Mann_BP_RP_nan', 'e_Teff_Mann_r_J_FeH', 'e_Teff_Mann_r_J_isJH',
                 'e_Teff_Mann_r_J_nan', 'e_Teff_Mann_r_z_FeH', 'e_Teff_Mann_r_z_isJH',
                 'e_Teff_Mann_r_z_nan', 'e_Teff_g_H_Boy', 'e_Teff_g_J_Boy', 'e_Teff_g_K_Boy',
                 'e_Teff_g_i_Boy', 'e_Teff_g_z_Boy', 'e_Teff_Apsis', 'Teff_Apsis',
                 'Teff_g_r_Boy', 'e_Teff_g_r_Boy']))

    assert np.array_equal(h.columns.values, ['EPIC', 'Teff_median', 'Teff_std'])
    assert hyades.stars[hyades.stars.Teff_median.isnull()].shape[0] ==  0
    print(h.Teff_median)
    assert (h.Teff_median.dropna() < 7000.).all()
    assert (h.Teff_median.dropna() > 3000.).all()

def test_find_Rstar():
    hyades = deepcopy(TestCluster)
    with pytest.raises(AttributeError) as e:
        hyades.find_Rstar(lib)
        assert str(e.value) == 'Teff_median and/or FeH not assigned in to any star.\nRun determine_Teff first.'
    hyades.determine_Teff()
    before = hyades.stars.columns.values
    hyades.find_Rstar(lib, mode="specmatchemp")
    after = hyades.stars.columns.values
    after = hyades.stars.columns.values
    assert np.array_equal(np.setdiff1d(after, before),
                         ['Rstar', 'Rstar_consistent', 'Rstar_double_check', 'e_Rstar', 'e_Rstar_double_check'])
    assert hyades.stars[hyades.stars.Rstar.isnull()].shape[0] ==  0
    assert hyades.stars[hyades.stars.e_Rstar.isnull()].shape[0] ==  0
    
    assert (hyades.stars.Rstar_consistent.dropna().values == np.array([True])).all()

    # Now test the mann mode
    hyades = deepcopy(TestCluster)
    hyades.determine_Teff()
    before = hyades.stars.columns.values
    hyades.find_Rstar(lib, mode="mann")
    after = hyades.stars.columns.values
    print(hyades.stars.Teff_median)
    assert np.array_equal(np.setdiff1d(after, before),
                         ['Rstar', 'Rstar_consistent', 'Rstar_double_check', 'e_Rstar', 'e_Rstar_double_check'])
    assert hyades.stars[hyades.stars.Rstar.isnull()].shape[0] ==  4
    assert hyades.stars[hyades.stars.e_Rstar.isnull()].shape[0] ==  4
    print(hyades.stars.Rstar_consistent)
    assert (hyades.stars.Rstar_consistent.dropna().values == np.array([True])).all()


def test_find_SED():
    hyades = deepcopy(TestCluster)
    with pytest.raises(AttributeError) as e:
        hyades.find_SED(lib)
    assert str(e.value) == 'Teff_median and/or FeH not assigned to any star.\nRun determine_Teff first.'
    hyades.determine_Teff()
    before = hyades.stars.columns.values
    hyades.find_SED(lib)
    after = hyades.stars.columns.values
    after = hyades.stars.columns.values
    assert np.array_equal(np.setdiff1d(after, before),
                          ['SED', 'SED_SJS', 'SED_SME',
                           'SpT', 'e_SED', 'e_SED_SJS',
                           'e_SED_SME', 'w_SED', 'w_SED_SJS',
                           'w_SED_SME'])
    assert hyades.stars[hyades.stars.SED_SME.isnull()].shape[0] ==  0
    assert hyades.stars[hyades.stars.e_SED_SME.isnull()].shape[0] ==  0
    assert hyades.stars[hyades.stars.w_SED_SME.isnull()].shape[0] ==  0
    assert (hyades.stars.w_SED_SME[0].value > 1000.).all()
    #test that masked wavelengths are replaced with 1. flux and 0. uncertainty
    assert np.where(~(hyades.stars.SED_SME[0] >= 0.))[0].shape[0] == 0
    assert np.where(~(hyades.stars.e_SED_SME[0] >= 0.))[0].shape[0] == 0

def test_calculate_luminosities():
    hyades = deepcopy(TestCluster)
    assert hyades.stars.shape[0] == 10
    hyades.determine_Teff()
    hyades.find_Rstar(lib)
    h= hyades.find_SED(lib)
    hyades.calculate_luminosities()
    #check if columns are there
    print(hyades.stars[["Rstar", "Teff_median"]])
    #print(hyades.stars.columns.values, np.isin(hyades.stars.columns.values, 'Lum_SED'))
    assert np.isin(hyades.stars.columns.values, 'Lum_SED').any()
    assert np.isin(hyades.stars.columns.values, 'Lum_SB').any()
    assert np.isin(hyades.stars.columns.values, 'Lum_Kepler').any()
    assert np.isin(hyades.stars.columns.values, 'planck_curve').any()
    assert np.isin(hyades.stars.columns.values, 'e_Lum_SED').any()
    assert np.isin(hyades.stars.columns.values, 'e_Lum_SB').any()
    assert np.isin(hyades.stars.columns.values, 'e_Lum_Kepler').any()
    assert np.isin(hyades.stars.columns.values, 'e_planck_curve').any()
    wonan = hyades.stars[['Lum_SED','Lum_SB','Lum_Kepler']].dropna(how='any')
    wonan = wonan.dropna(how="any")
    #two targets have no values yet
    assert wonan.shape[0] == 10

    #loose check if values have a reasonable ratio
    assert (wonan.Lum_Kepler < wonan.Lum_SED).all()
    assert (wonan.Lum_SED < wonan.Lum_SB).all()
    #stricter check if values have a reasonable ratio
    assert (wonan.Lum_SED/wonan.Lum_SB > 0.1).all()
    assert (wonan.Lum_Kepler/wonan.Lum_SED > 0.35).all()
    assert (wonan.Lum_Kepler/wonan.Lum_SB > 0.1).all()

def test_merge_flares_stars():

    analysis_cluster.merge_flares_stars()
    assert (analysis_cluster.stars.columns.values.shape[0] +
            analysis_cluster.flares.columns.values.shape[0] ==
            analysis_cluster.flares_on_stars.shape[1]+2)
    assert analysis_cluster.flares_on_stars.shape[0] == analysis_cluster.flares.shape[0]
    assert analysis_cluster.flares_on_stars.EPIC.equals(pd.Series([210942999,210942999,210563410]))

def test_generate_FFD():
    analysis_cluster.merge_flares_stars()
    a = deepcopy(analysis_cluster)

    print(a.stars[["Teff_median", "EPIC"]], "TEFFEpIC\n\n")
    #everything with 'erg'
    a.unit = 'erg'
    a.generate_FFD(Tlims=(2500*u.Kelvin,3330*u.Kelvin))
    ffd = a.ffds[2500]
    print(ffd)
    assert ffd.Tmin == 2500*u.Kelvin
    assert ffd.Tmax == 3330*u.Kelvin
    assert ffd.flares.shape[0] == 3
    assert ffd.name == 'hyades'
    assert ffd.h_name == 'Hyades'
    assert ffd.ED.shape[0] == 3
    print(ffd.ED)
    assert (ffd.ED==np.array([40, 400, 50])).all()
    assert ffd.tot_obs_time == 2000./24./365.25
    assert ffd.agemin == ffd.agemax
    assert ffd.agemin == 600.
    assert ffd.alpha is None
    assert ffd.beta is None
    assert ffd.alpha_err is None
    assert ffd.beta_err is None
    assert ffd.multiplicity is None
    assert ffd.cutoff_energy_lower is None
    assert ffd.cutoff_energy_upper is None
    assert ffd.cutoff_ED_lower is None
    assert ffd.cutoff_ED_upper is None
    assert ffd.energy_err is None
    assert ffd.ED_err is None
    assert ffd.mode is 'energy'
    assert ffd.nstars == 4
    cols = ['Lum_SB','Lum_Kepler','Lum_SED','Teff_median','EPIC','Rstar','Campaign','istart']
    df = analysis_cluster.flares_on_stars[cols]
    df = df.iloc[[0,1],:]
    compareto = ffd.flares[cols]
    compareto.Lum_Kepler = compareto.Lum_Kepler


    print(df.iloc[0,:],compareto.iloc[0,:])
    print(compareto -df)
    assert (compareto - df).sum().sum() == 0.
    #everything with 's'
    a = deepcopy(analysis_cluster)
    a.unit = 's'
    a.generate_FFD(Tlims=(2500*u.Kelvin,3330*u.Kelvin))
    ffd = a.ffds[2500]
    assert ffd.Tmin == 2500*u.Kelvin
    assert ffd.Tmax == 3330*u.Kelvin
    assert ffd.flares.shape[0] == 3
    assert ffd.name == 'hyades'
    assert ffd.h_name == 'Hyades'
    assert ffd.ED.shape[0] == 3
    assert (ffd.ED==np.array([40, 400, 50])).all()
    assert ffd.tot_obs_time == 2000./24./365.25 # this is likely to change but should be a multiple of 500, like 2000
    assert ffd.agemin == ffd.agemax
    assert ffd.agemin == 600.
    assert ffd.alpha is None
    assert ffd.beta is None
    assert ffd.alpha_err is None
    assert ffd.beta_err is None
    assert ffd.multiplicity is None
    assert ffd.cutoff_energy_lower is None
    assert ffd.cutoff_energy_upper is None
    assert ffd.cutoff_ED_lower is None
    assert ffd.cutoff_ED_upper is None
    assert ffd.energy_err is None
    assert ffd.ED_err is None
    assert ffd.mode is 'ED'
    assert ffd.nstars == 4
    cols = ['Lum_SB','Lum_Kepler','Lum_SED','Teff_median','EPIC','Rstar','Campaign','istart']
    df = analysis_cluster.flares_on_stars[cols]
    df = df.iloc[[0,1],:]
    compareto = ffd.flares[cols]
    assert (compareto - df).sum().sum() == 0.

def test_generate_all_FFDs():
    clusters = CLUSTERS.append({'cluster':'fake',
                        'h_cluster':'Fake Cluster',
                        'age (Myr)': 1000.,
                        'FeH' : 0.05,
                        'dist (pc)': 200,
                        'u_age_high':200,
                        'u_age_low':100,
                        'u_feh':0.001,
                        'Teffmax':7000}, ignore_index=True)
    #generate fake OC with observations:
    OC = generate_OpenCluster('fake', clusters=clusters, unit='erg', path=PACKAGEDIR)
    #add this new step of the analysis and add a mode:
    OC = OC.generate_all_FFDs()
    #test results
    assert len(set(OC.ffds)) == 8
    assert (OC.ffds[3500].flares.Teff_median == OC.flares_on_stars.Teff_median).all()
    assert OC.ffds[3250].ED is None
    assert OC.ffds[5000].ED is None
    with pytest.raises(KeyError) as e:
        OC.ffds[5560]
    ffd = OC.ffds[3500]
    assert ffd.name == 'fake'
    assert ffd.h_name == 'Fake Cluster'
    assert ffd.age == 1000.
    assert ffd.agemin == 900
    assert ffd.agemax == 1200.
    assert ffd.cutoff is None
    assert ffd.mode is 'energy'
    assert ffd.ED == pytest.approx(np.array([17.71664063, 10.26588968, 22.17849684, 17.69557813]))

    #Now everything with ED:
    #generate fake OC with observations:
    OC = generate_OpenCluster('fake', clusters=clusters, unit='s', path=PACKAGEDIR)
    OC.unit = 's'
    OC = OC.generate_all_FFDs()
    #test results
    assert len(set(OC.ffds)) == 8
    assert (OC.ffds[3500].flares.Teff_median == OC.flares_on_stars.Teff_median).all()
    assert OC.ffds[3250].ED is None
    assert OC.ffds[5000].ED is None
    with pytest.raises(KeyError) as e:
        OC.ffds[5560]
    ffd = OC.ffds[3500]
    assert ffd.name == 'fake'
    assert ffd.h_name == 'Fake Cluster'
    assert ffd.age == 1000.
    assert ffd.agemin == 900
    assert ffd.agemax == 1200.
    assert ffd.cutoff is None
    assert ffd.mode is 'ED'
    assert ffd.ED == pytest.approx(np.array([17.71664063, 10.26588968, 22.17849684, 17.69557813]))
