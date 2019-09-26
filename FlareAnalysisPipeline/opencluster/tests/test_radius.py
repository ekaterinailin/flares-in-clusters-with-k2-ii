import pytest
import pandas as pd
import numpy as np

from ..radius import (radius_specmatchemp,
                      radius_mann,
                      calculate_double_check_radii_mann)
import specmatchemp.library

def test_calculate_double_check_radii_mann():
    # --------------------------------------------------------------------
    # Produce fake data:

    testdf = pd.DataFrame({"K_2MASS":[6.,27.,15.,16.728],
                           "e_K_2MASS":[.06,.07,0.01,.094],
                          "parallax_Gaia":[1.3,2.3,3.3,3.3],
                          "parallax_error_Gaia":[.001,.001,.001,.001,],
                          "FeH":[np.nan,.1,np.nan,.1],
                          "e_FeH":[.03,np.nan,.03,.03],}
                          )
    # --------------------------------------------------------------------
    # Call the function

    res = calculate_double_check_radii_mann(testdf)

    # --------------------------------------------------------------------
    # Do some checks

    # Check that limits on MKs are applied properly
    assert np.isnan(res.loc[0,"Rstar_double_check"])
    assert np.isnan(res.loc[0,"e_Rstar_double_check"])
    assert np.isnan(res.loc[1,"Rstar_double_check"])
    assert np.isnan(res.loc[1,"e_Rstar_double_check"])

    # Check that missing metallicity does not let the function fail to produce a result:
    assert ~np.isnan(res.loc[2,"Rstar_double_check"])
    assert ~np.isnan(res.loc[2,"e_Rstar_double_check"])

    # Check against value from famous Mamajek table:
    assert res.loc[3,"Rstar_double_check"] == pytest.approx(0.127, rel=0.08)

def test_radius_mann():
    df = pd.DataFrame({'Teff_std' : [500,300,200,250,200],
                       'Teff_median' : [2550,2750,2750,3000,3600],
                       'EPIC' : [211119999, 211119996, 211119995, 211119997, 211119998],
                       'FeH': [.01, .01, np.nan, .01, .01],
                       'e_FeH': [.001, .001, np.nan, np.nan, .001]})
    res = radius_mann(df)
    assert res.Rstar.fillna(0).iloc[0] == 0.
    print(res.Teff_median.iloc[3])
    # Checked values with Pecaut and Mamajek, should be correct within uncertainties:
    assert res.Rstar.iloc[3] == pytest.approx(0.149, abs=res.e_Rstar.iloc[3])
    assert res.Rstar.iloc[1] == pytest.approx(0.128, abs=res.e_Rstar.iloc[1])
    assert res.Rstar.iloc[2] == pytest.approx(0.128, abs=res.e_Rstar.iloc[2])

def test_radius_specmatchemp():
    lib = specmatchemp.library.read_hdf(wavlim=[4000,9000])
    df = pd.DataFrame({'Teff_std' : [500,0,200,250,200],
                   'Teff_median' : [2550,2750,2750,3000,4600],
                   'EPIC' : [211119999, 211119996, 211119995, 211119997, 211119998],
                   'FeH': [.01, .01, np.nan, .01, .01],
                   'e_FeH': [.001, .0, np.nan, np.nan, .001]})
    df = radius_specmatchemp(df, lib)
    assert "Rstar" in df.columns
    assert "e_Rstar" in df.columns
    print(df.e_Rstar)
    assert (~np.isnan(df.e_Rstar.values[1:])).all()
    assert np.isnan(df.e_Rstar.values[0])
    assert df.e_Rstar.iloc[1] == 0.093 * df.Rstar.iloc[1]
    assert df.Rstar.iloc[-1] == pytest.approx(0.726, abs=df.e_Rstar.iloc[-1])