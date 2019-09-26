import numpy as np
import pandas as pd
import pytest

from ..rotation import (rossby_wright,
                        prepare_rotation_activity_table)

def test_rossby_wright():

    # check too low V-K
    with pytest.raises(ValueError):
        rossby_wright(10, 1.)

    # check expected results from generic input
    assert np.isnan(rossby_wright(np.nan, 10))
    assert rossby_wright(0, 10) == 0.
    assert rossby_wright(1, 10.) == pytest.approx(1.1735108709918103)

def test_prepare_rotation_activity_table():

    # Start with basic generic tables:
    cluster = "generic"
    EPIC = [211119999, 211119998, 211119988, 211119888]
    Prot = [.8, 3, 15]
    Campaign = [4, 5, 16, 18]
    VK = [3.8, 4.7, 6.1]
    ed_rec = [1000, 2000, 400, 600, 67, 1400]
    df = pd.DataFrame({"Campaign": Campaign,
                       "EPIC": EPIC,})
    dff = pd.DataFrame({"EPIC": EPIC + [211119999, 211119999],
                        "ed_rec": ed_rec})
    rot = pd.DataFrame({"Prot": Prot,
                        "EPIC": EPIC[:-1],
                        "V-K" : VK})

    #----------------------------------------------------------

    # run the func on the most complete version
    df2 = prepare_rotation_activity_table(cluster, df, dff, rot=rot, path=False,
                                        EDthreshold=300)

    # do some checks
    assert df2.shape == (3,9)
    print(df2.nflares)
    assert df2.loc['0', "nflares"] == 2 # sum all flares
    assert df2.loc['0', "ed_rec_sum"] == 2400 # sum all flares
    assert np.array_equal(df2.hasflare.values, np.array([1,1,1]))
    assert np.array_equal(df2.EPIC.values, np.array([211119999,211119998,211119988]))

    #----------------------------------------------------------

    # now truncate the flare sample such that
    # there is a star with rotation but no flares

    dff = pd.DataFrame({"EPIC": EPIC[1:],
                    "ed_rec": [310, 500, 1000]})

    # run the func on the modified version
    df2 = prepare_rotation_activity_table(cluster, df, dff, rot=rot, path=False,
                                          EDthreshold=300)

    # do some checks
    assert np.array_equal(df2.hasflare.values, np.array([0,1,1]))
    assert df2.loc['0', "nflares"] == 0 # sum all flares
    assert df2.loc['0', "ed_rec_sum"] == 0 # sum all flares
    assert df2.shape == (3,9)
