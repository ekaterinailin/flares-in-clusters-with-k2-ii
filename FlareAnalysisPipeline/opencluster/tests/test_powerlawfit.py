import pandas as pd
import numpy as np
import astropy.units as u
import pytest

from ..powerlawfit import (lnHurwitz_Ceta,
                           logL_powerlaw,
                           powerlawfit_Bauke2007,
                           de_biased_upper_limit,
                           de_bias_alpha,
                           ML_powerlaw_estimator,)
from ..utils import generate_random_power_law_distribution

def test_ML_powerlaw_estimator():
    dataformat = [np.array, pd.Series]
    for df in dataformat:
        EDs = df([1,1,1,1,2,2,4])
        a = 1.
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(EDs, a)
        a = 1.5
        assert ML_powerlaw_estimator(EDs, a) == pytest.approx(1.6190804181576444)
        EDs = df([])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(EDs, a)
        EDs = df([-3,2,2,2])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(EDs, a)
        EDs = df([1,1,1,1])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(EDs, a)

def test_de_biased_upper_limit():
    dataformat = [np.array, pd.Series]
    for df in dataformat:
        data = df([])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 2.)
        data = df([1,10, 100])
        assert de_biased_upper_limit(data, 1000000.) == pytest.approx(100., rel=1e-4)
        data = df([1,1,1])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 3.)
        data = df([-1,1,1])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 3.)

def test_de_bias_alpha():
    assert de_bias_alpha(200,1) == 1.
    with pytest.raises(ValueError):
        de_bias_alpha(np.nan,2)
    with pytest.raises(ValueError):
        de_bias_alpha(30,np.nan)
    with pytest.raises(ValueError):
        de_bias_alpha(np.nan,np.nan)
    with pytest.raises(ZeroDivisionError):
        de_bias_alpha(2,2)

def test_lnHurwitz_Ceta():
    with pytest.raises(ValueError):
        lnHurwitz_Ceta(4,1,1)
    with pytest.raises(ValueError):
        lnHurwitz_Ceta(1.,1,2)
    with pytest.raises(ValueError):
        lnHurwitz_Ceta(2, -1, 9)
    with pytest.raises(ValueError):
        lnHurwitz_Ceta(2, -1, 1)
    with pytest.raises(ValueError):
        lnHurwitz_Ceta(2, 2, 1)
    assert lnHurwitz_Ceta(2, 1, 2) == 0. #checked with wolfram alpha

def test_logL_powerlaw():
    rr = generate_random_power_law_distribution(10,100,-1,seed=90,size=50)
    data = np.rint(rr*1e10).astype(int)
    with pytest.raises(ValueError):
        logL_powerlaw(np.nan, data)

    data[4] = -3
    with pytest.raises(ValueError):
        logL_powerlaw(3, data)

    data[4] = 0
    with pytest.raises(ValueError):
        logL_powerlaw(3, data)

    with pytest.raises(TypeError):
        logL_powerlaw(3,3)

def test_powerlawfit_Bauke2007():
    rr = generate_random_power_law_distribution(10, 100, -1, seed=90, size=500)
    data = np.rint(rr * 1e10).astype(int)
    m = powerlawfit_Bauke2007(data, minexp=1.+1e-8, maxexp=3.)
    assert m == pytest.approx(2.03,rel=.01)
    with pytest.raises(ValueError):
        powerlawfit_Bauke2007(data, minexp=3, maxexp=3.)