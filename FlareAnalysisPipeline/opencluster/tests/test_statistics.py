import pandas as pd
import numpy as np
import astropy.units as u
import pytest

from opencluster.ffd import FFD
from ..statistics import (calculate_cumulative_powerlaw_distribution,
                          _apply_stabilising_transformation,
                          calculate_KS_acceptance_limit,
                          stabilised_KS_statistic,
                          _calculate_number_of_exceeding_values,
                          calculate_average_number_of_exceeding_values)
from ..utils import generate_random_power_law_distribution

def test_calculate_average_number_of_exceeding_values():
    data = np.linspace(10,1e4, 300)
    ffd = FFD(ED=data, alpha=2.)
    mean, std = calculate_average_number_of_exceeding_values(ffd, 1000, seed=2311)
    assert mean == 3.
    assert std == 0.
    mean, std = calculate_average_number_of_exceeding_values(ffd, 1000, seed=10)
    assert mean == 0.
    assert std == 0.
    with pytest.raises(AssertionError):
        ffd = FFD(ED=data)
        mean, std = calculate_average_number_of_exceeding_values(ffd, 1000, seed=2311)
    with pytest.raises(ValueError):
        ffd = FFD()
        calculate_average_number_of_exceeding_values(ffd, 1000, seed=2311)

def test__calculate_number_of_exceeding_values():
    data = np.linspace(10,1e4, 300)
    assert _calculate_number_of_exceeding_values(data, 2., seed=10) == 0
    assert _calculate_number_of_exceeding_values(data, 2., seed=2311) == 3
    with pytest.raises(ValueError):
        _calculate_number_of_exceeding_values(np.arange(3), 2., seed=2311)

def test_stabilised_KS_statistic():
    sizes = [1e2,1e3,1e4]
    minval, maxval = 10, 1e4
    datas = [generate_random_power_law_distribution(minval, maxval, -1., size=int(size), seed=10) for size in sizes]
    KSlist = [stabilised_KS_statistic(data, alpha=2.) for data in datas]
    assert KSlist[0] > KSlist[1]
    assert KSlist[1] > KSlist[2]

def test_calculate_KS_acceptance_limit():
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=0.05)
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=-0.05)
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=1.05)
    assert (calculate_KS_acceptance_limit(100, sig_level=0.05)
            == pytest.approx(0.13581, rel=1e-4))
    assert (calculate_KS_acceptance_limit(100, sig_level=0.01)
            == pytest.approx(0.16276, rel=1e-4))
    assert (calculate_KS_acceptance_limit(100, sig_level=0.01)
            > calculate_KS_acceptance_limit(1000, sig_level=0.01))

def test__apply_stabilising_transformation():
    _apply_stabilising_transformation(1.) == 1.
    _apply_stabilising_transformation(0.) == -1.
    with pytest.raises(AssertionError):
        _apply_stabilising_transformation(-.1)

def test_calculate_cumulative_powerlaw_distribution():
    sizes = [1e2,1e3,1e4]
    minval, maxval = 10, 1e4
    datas = [generate_random_power_law_distribution(minval, maxval, -1., size=int(size), seed=105) for size in sizes]
    CDFlist = np.array([calculate_cumulative_powerlaw_distribution(data, alpha=2.) for data in datas])
    for i in range(3):
        assert CDFlist[i][0] == 0.
        assert CDFlist[i][-1] == 1.
        assert (CDFlist[i]>=0.).all()
        assert (CDFlist[i]<=1.).all()
        assert (np.diff(CDFlist[i]) > 0.).all()

    with pytest.raises(ValueError):
        calculate_cumulative_powerlaw_distribution(datas[0], alpha=1.)