import pandas as pd
import numpy as np
import astropy.units as u
import pytest

from ..linfit import (linfit,
                      define_linfit_arrays,
                      logL_lin)
from ..ffd import FFD

def test_define_linfit_arrays():
    #typical input
    df = pd.DataFrame({'ed_rec' : np.array([1,1,1,1,2,2,4]),
                      'Lum_Kepler' : np.ones_like(7)})
    df.ed_rec = df.ed_rec
    df.Lum_Kepler = df.Lum_Kepler
    obs_duration = 1.
    a, freq= define_linfit_arrays(df, 'ED')
    assert (a == df.ed_rec.values).all()
    print(freq)
    assert (freq == np.arange(7,0,-1)).all()
    a, freq = define_linfit_arrays(df, 'energy')
    assert (a == df.ed_rec.values).all()
    assert (freq == np.arange(7,0,-1)).all()
    #input with a nan
    df = pd.DataFrame({'ed_rec' : np.array([1,1,1,np.nan,2,2,4]),
                      'Lum_Kepler' : np.ones_like(7)})
    df.ed_rec = df.ed_rec
    df.Lum_Kepler = df.Lum_Kepler
    a, freq = define_linfit_arrays(df, 'ED')
    assert (a == df.ed_rec.dropna().values).all()
    assert (freq == np.arange(6,0,-1)).all()

    #input with a nan in Lum_Kepler
    df.loc[2,'Lum_Kepler'] = np.nan
    a, freq = define_linfit_arrays(df, 'energy')
    assert (freq == np.arange(5,0,-1)).all()
    assert (a == np.array([1., 1., 2., 2., 4.])).all()

def test_linfit():
    y = [8.,6.,4.,2.]
    obs_time = 10.
    data = pd.DataFrame({'x':[1.,2.,3.,4.,],
                                'y':y,
                                'sigy':np.sqrt(y),
                                'sigx':0,
                                'rho': 0,
                                'norm':np.log10(obs_time)})
    m, b = linfit(data)
    assert m[0] == pytest.approx(-2.)
    assert m[1] == pytest.approx(2.3e-5,rel=.1)
    assert b[0] == pytest.approx(9.)
    assert b[1] == pytest.approx(4.7e-5,rel=.1)

    data.loc[0,'x'] = np.nan
    m, b = linfit(data)
    assert m[0] == pytest.approx(-2., rel=1e-3)
    assert m[1] == pytest.approx(4.8e-6,rel=.1)
    assert b[0] == pytest.approx(9.)
    assert b[1] == pytest.approx(2.3e-5,rel=.1)

def test_logL_lin():
    args = ([1,2,3], [4,5,6],[.1,.1,.1], [0,0,0],  [0,0,0])
    Ls = np.array([logL_lin([i,1], *args) for i in range(10)])
    assert np.argmin(Ls) == 2
    with pytest.raises(ValueError):
        logL_lin([np.nan,1], *args) #my function to be tested
    with pytest.raises(ValueError):
        logL_lin([6,np.nan], *args)
    args = ([1,2,np.nan], [4,5,6],[.1,.1,.1], [0,0,0],  [0,0,0])
    with pytest.raises(ValueError):
        logL_lin([6,1.], *args)