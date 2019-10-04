import pandas as pd
import numpy as np
import pytest
from astropy import units as u

from ..lum import (read_Kepler_response,
                   planck_full_blackbody,
                   planck,
                   interpolate_spectrum,
                   Kepler_mean_wavelength,
                   projected_luminosity_Stefan_Boltzmann,
                   projected_luminosity_SED,
                   projected_luminosity_Kepler,
                   SME_find_library_match,
                   assign_SpT_to_Teff,
                   remove_outlier_Teffs)

def test_remove_outlier_Teffs():
    teffs = [5000.] * 10
    df = pd.DataFrame({'todrop' :['g-J outlier, BP-RP outlier',
                                  np.nan, 'BP-RP outlier', np.nan,
                                  np.nan, np.nan, np.nan, np.nan,
                                  np.nan, 'r-i outlier'],
                       'Teff_Mann_r_i': teffs,
                       'Teff_BP_RP_3': teffs,
                       'Teff_median': teffs,})
    remove_outlier_Teffs(df)
    assert np.isnan(df.Teff_BP_RP_3[0])
    assert np.isnan(df.Teff_BP_RP_3[2])
    assert np.isnan(df.Teff_BP_RP_3[2])
    assert np.isnan(df.Teff_Mann_r_i[9])
    assert (df.Teff_median.values == 5000.).all()

def test_assign_SpT_to_Teff():
    df = pd.DataFrame({'Teff_std' : [500,300,200,250,200],
                       'Teff_median' : [1540,2750,2750,3000,4600],
                       'EPIC' : [211119999, 211119996, 211119995, 211119997, 211119998]})
    df = assign_SpT_to_Teff(df)
    assert (df.SpT.values == np.array(["L5","M6","M6","M5",""])).all()

def test_SME_find_library_match():
    testcut = pd.DataFrame({'Teff' : [3400., 3600.,3800.,],
                            'Terr' : [200, 200, 200],
                            'feh' : [.1,.2,.3]})

    T, feh, Terr = 3400., .6, 200
    assert np.isnan(SME_find_library_match(testcut, T, Terr, feh)).all()
    T, feh, Terr = 3500., .15, 200
    
    assert SME_find_library_match(testcut, T, Terr, feh)[0] == 0
    assert SME_find_library_match(testcut, T, Terr, feh)[1].equals(pd.Series({'Teff': 3400.0, 'feh': 0.1}))
    T, feh, Terr = np.nan, .15, Terr
    assert np.isnan(SME_find_library_match(testcut, T, Terr, feh)).all()

def test_read_Kepler_response():
    Kp = read_Kepler_response()
    assert Kp.wav.min()==3480.
    assert Kp.wav.max()==9700.

def test_planck_full_blackbody():
    wav, pl = planck_full_blackbody(0.)
    assert pl.sum() == 0.
    wav, pl = planck_full_blackbody(np.nan)
    assert np.isnan(pl).all()
    wav, pl = planck_full_blackbody(1e4)
    assert wav[np.argmax(pl)].to('angstrom').value == pytest.approx(2.8977729e-7*u.m.to('angstrom'), abs=2)#Wien's law

def test_planck():
    w = np.arange(300,10000)*u.nm
    p = planck(w, 1e4)
    assert w[np.argmax(p)].value == 300. #max at <300nm
    assert np.isnan(planck(w, np.nan)).all()
    assert planck(w, 0.).sum() == 0.
    assert w[np.argmax(planck(w, 5777))].value == pytest.approx(501.605, abs=2)#Wien's law

def test_interpolate_spectrum():
    target_wav = np.arange(2,10)
    spectrum_wav = np.linspace(0,20,41)
    spectrum = spectrum_wav**2
    target_spectrum = interpolate_spectrum(target_wav, spectrum_wav, spectrum,)
    assert target_spectrum.shape[0] == target_wav.shape[0]
    assert np.min(target_spectrum) == 4.
    assert np.max(target_spectrum) == 81.

def test_Kepler_mean_wavelength():
    Kp = pd.DataFrame({'wav' : np.arange(3000,3030),
                       'resp' : np.linspace(0,1,30)})
    res = Kepler_mean_wavelength(Kp)
    assert min(res.to('angstrom')) == 3000.*u.angstrom
    assert res.unit == u.cm
    assert res[-1].value == pytest.approx(3.0285e-5)
    assert res[0].value == pytest.approx(3.0e-5)

def test_projected_luminosity_Stefan_Boltzmann():
    # Check main functionality
    assert 4 * projected_luminosity_Stefan_Boltzmann(5777., R=1.) == pytest.approx(3.84e33, abs=1e31)
    assert projected_luminosity_Stefan_Boltzmann(0., R=1.)== 0.
    assert projected_luminosity_Stefan_Boltzmann(5777., R=0.)== 0.
    assert projected_luminosity_Stefan_Boltzmann(0., R=0.)== 0.
    assert np.isnan(projected_luminosity_Stefan_Boltzmann(np.nan, R=0.))
    assert np.isnan(projected_luminosity_Stefan_Boltzmann(np.nan, R=np.nan))
    assert np.isnan(projected_luminosity_Stefan_Boltzmann(0., R=np.nan))

    #Check the uncertainties
    assert ((projected_luminosity_Stefan_Boltzmann(1., R=1.,
                                                   e_R=.1, e_T=0.1, deriv=True) /
             projected_luminosity_Stefan_Boltzmann(1., R=1.)) ==
            0.2 *np.sqrt(5.))

def test_projected_luminosity_SED():
    SED = np.array([1,2,4,8,16,5,1,0.5,0.25,0.1,])#*u.erg/u.s/u.cm**3
    spectrum = np.array([1., .9, .6, .1, .6, .9, 1., 1., 1., 1.,])
    R = 1
    e_R = .2
    e_SED = 0.01*SED
    wavelengths = np.arange(1,11)*100*u.angstrom
    assert projected_luminosity_SED(np.zeros(10), spectrum, wavelengths, R)== 0.
    assert projected_luminosity_SED(SED, spectrum, wavelengths, 0.)== 0.
    assert projected_luminosity_SED(SED, np.zeros(10), wavelengths, R)== 0.
    assert np.isnan(projected_luminosity_SED(SED*np.nan, np.zeros(10), wavelengths, R))

    # Test some back-of-the-envelope calculations
    assert ((projected_luminosity_SED(SED, spectrum, wavelengths, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED) /
             projected_luminosity_SED(SED, spectrum, wavelengths, R,
                                        deriv=False, e_R=e_R, e_planck=e_SED))==
            pytest.approx(.4+1e-4, rel=1e-3))
    assert ((projected_luminosity_SED(SED, spectrum, wavelengths, R,
                                        deriv=True, e_R=0, e_planck=e_SED) /
             projected_luminosity_SED(SED, spectrum, wavelengths, R,
                                        deriv=False, e_R=0, e_planck=e_SED))==
            pytest.approx(3.5e-3, rel=.2))


def test_projected_luminosity_Kepler():
    SED = np.array([1,2,4,8,16,5,1,0.5,0.25,0.1,])#*u.erg/u.s/u.cm**3
    e_SED = SED * .01
    spectrum = np.array([1., .9, .6, .1, .6, .9, 1., 1., 1., 1.,])
    spectrum2 = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,])
    R = 1
    e_R = .2
    wavelengths = np.arange(1,11) * 100 * u.angstrom
    resp = pd.Series([1., .9, .6, .1, .6, .9, 1., 1., 1., 1.,])
    resp2 = pd.Series([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,])
    # Test the luminosities
    assert projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R) < projected_luminosity_Kepler(SED, spectrum, wavelengths, resp2, R)
    assert projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R) < projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R)
    # Now test the uncertainties
    assert (projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED) <
            projected_luminosity_Kepler(SED, spectrum, wavelengths, resp2, R,
                                          deriv=True, e_R=e_R, e_planck=e_SED))
    assert (projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED) <
            projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED))
    ## Check a back of the envelope calculation for the given values:
    assert ((projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED) /
             projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R,
                                        deriv=False, e_R=e_R, e_planck=e_SED)) ==
            pytest.approx(.4+1e-4, rel=1e-3))
    assert ((projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=True, e_R=e_R, e_planck=e_SED) /
             projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=False, e_R=e_R, e_planck=e_SED)) ==
            pytest.approx(.4+1e-4, rel=1e-3))
    assert ((projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=True, e_R=0, e_planck=e_SED) /
             projected_luminosity_Kepler(SED, spectrum, wavelengths, resp, R,
                                        deriv=False, e_R=0, e_planck=e_SED)) ==
            pytest.approx(3.5e-3, rel=.2))
    assert ((projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R,
                                        deriv=True, e_R=0, e_planck=e_SED) /
             projected_luminosity_Kepler(SED, spectrum2, wavelengths, resp, R,
                                        deriv=False, e_R=0, e_planck=e_SED)) ==
            pytest.approx(3.6e-3, rel=.2))

#intergration test on solar values
def test_solar_luminosity():
    R = 1.
    T = 5777.
    spectrum = np.ones(623)
    Kp = read_Kepler_response()
    Kepler_wav = Kp.wav.values*u.angstrom.to('cm')*u.cm
    Planck_curve = planck(Kepler_wav, T)

    lSED_Sun = 4 * projected_luminosity_SED(Planck_curve, spectrum, Kepler_wav, R)
    lKepler_Sun = 4 * projected_luminosity_Kepler(Planck_curve, spectrum, Kepler_wav, Kp.resp, R)
    lSB_Sun = 4 * projected_luminosity_Stefan_Boltzmann(T, R)

    BB_wav, BB_curve = planck_full_blackbody(T)
    spectrum2 = np.ones(10000)
    lBB_Sun = 4 * projected_luminosity_SED(BB_curve, spectrum2, BB_wav, R)
    assert lSB_Sun-lBB_Sun < 1e31
    assert lSB_Sun > lSED_Sun
    assert lSED_Sun == pytest.approx(2.43614943e33, rel=1e-2)
    assert lSB_Sun == pytest.approx(3.8483599264e33, rel=1e-2)
    assert lKepler_Sun == pytest.approx(1.0909422313e33, rel=1e-2)
