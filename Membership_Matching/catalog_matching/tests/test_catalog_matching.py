import os
import pytest

import pandas as pd
import numpy as np

from collections import OrderedDict

from ..catalog_matching import (crossmatch,
                                select_min_dist,
                                post_k2_clean,
                                find_campaigns,
                                match_k2_epic_campaign,
                                extract_extensions,
                                assign_PMem_mean,
                                merge_k2_info_and_protocol,
                                crossmatch_multiple_catalogs,
                                pick_members_and_produce_k2_search_input,
                                )



ra = 'RAJ2000'
dec = 'DEJ2000'

def test_pick_members_and_produce_k2_search_input():
    #---------------------------------------------------------
    # Produce fake data

    cross = pd.DataFrame(dict(zip(["RAJ2000_1", "DEJ2000_1",
                                   "somecolumn", "PMem_1",
                                   "PMem_2", "PMem_3"],
                                   [[20, 20, 20],
                                   [20, 20, 20],
                                   ["rolf", "hagen", "busch"],
                                   [.1, .8, .9],
                                   [.1, .8, .9],
                                   [.9, np.nan, .9]],
                                 )))
    sname, name = "test", "Test"
    coords = "1"
    series = cross.loc[2,:] #this row should be preserved in result
    outfile = ('catalog_matching/matched_catalogs/'
              'membership_matches/radec/{}_radec.csv'
              .format(sname))

    #---------------------------------------------------------
    # Call function

    res = pick_members_and_produce_k2_search_input(cross, sname,
                                                   name, coords=coords)
    df = pd.read_csv(outfile, header=None)

    #---------------------------------------------------------
    # Check if the RA/Dec file is correct:

    assert df.loc[0,0] == 20
    assert df.loc[0,1] == 20
    assert df.shape[0] == 1
    assert df.shape[1] == 2

    # Remove output file

    os.remove(outfile)

    # Check if the DataFrame was processed correctly

    assert res.shape[0] == 1 # only one member is left
    assert (res.loc[2,series.index] == series).all() # input columns are preserved


def test_crossmatch_multiple_catalogs():
    #-----------------------------------------------------------
    # Create a fake data set

    diff = 1.49/3600 # 1.5 arcsec distance
    c1 = pd.DataFrame(dict(zip(["RAJ2000_1","DEJ2000_1","PMem_1"],
                               [[21,20,19],[10,10,10],[.9,.8,.7]])))
    c2 = pd.DataFrame(dict(zip(["RAJ2000_2","DEJ2000_2","PMem_2","binary_2"],
                               [[21,20+diff,19+3*diff],[10,10,10],
                                [.75,.85,.3],[.1,.02,.11]])))
    c3 = pd.DataFrame(dict(zip(["RAJ2000_3","DEJ2000_3","PMem_3"],
                               [[np.nan,20-diff,19],[10,10,10],[.9,.9,.9]])))

    d = {"1":c1, "2":c2, "3":c3}
    renamed_catalogs = OrderedDict(sorted(d.items(), key=lambda t: t[0])) # order the dicts, not necessary for performance but helpful for testing
    name = "Test"
    sname = "test"

    #-----------------------------------------------------------
    # Call the function

    res = crossmatch_multiple_catalogs(renamed_catalogs, name, sname,
                                       arcsec=3., plot=True, union=True,
                                       bijective=True)
    #-----------------------------------------------------------
    # Do some checks

    # Check that the table size is correct
    assert res.shape[0] == 5
    assert res.shape[1] == 16

    # Check that relevant columns are created with the right names/values

    assert  "DEJ2000_1_2_3" in res.columns.tolist()
    assert set(c1.columns.values).issubset(set(res.columns.values))
    assert set(c2.columns.values).issubset(set(res.columns.values))
    assert set(c3.columns.values).issubset(set(res.columns.values))

    # Check that the distance calculation was done correctly

    assert res.loc[1, "dist_1_2_3"] == pytest.approx(2.235, rel=.1)
    assert res.loc[2, "dist_1_2_3"] == 0.

    # Check that NaNs stay NaNs:

    assert np.isnan(res.loc[4, "RAJ2000_3"])

    # Check individual values and NaNs:

    assert res.loc[2, "RAJ2000_1_2_3"] == 19.
    assert (res.DEJ2000_1_2_3.values == 10.).all()
    assert res.dist_1_2.dropna().shape[0] == 2
    assert res.dist_1_2_3.dropna().shape[0] == 2


def test_merge_k2_info_and_protocol():

    # Fake data

    folder = "catalog_matching/tests/exfiles/"
    sname = "tiny"
    df = pd.read_csv('catalog_matching/tests/exfiles/select_min_dist_union_k2.csv')
    racols = list(filter(lambda k: ra in k, df.columns.values))
    deccols = list(filter(lambda k: dec in k, df.columns.values))
    df["{}_mean".format(ra)] = df[racols].mean(axis=1)
    df["{}_mean".format(dec)] = df[deccols].mean(axis=1)


    # Call function

    merge_k2_info_and_protocol(df, "tiny", "mean", folder=folder)

    data = pd.read_csv("{}members.txt".format(folder))

    assert data.shape[1] == 6
    assert data.shape[0] == 1
    assert data.nlikelymembers[0] == 7
    assert data.LCs[0] == 5
    assert data.cluster[0] == "tiny"

    # Remove the test file, as it would get longer and longer over time.

    os.remove("{}members.txt".format(folder))

def test_assign_PMem_mean():

    # Test a working example:

    df = pd.DataFrame(index=np.arange(10),
                      columns=["RAJ2000_1","DEJ2000_1","PMem_1",
                               "RAJ2000","DEJ2000_holla","PMem_holla",
                               "RAJ2000_2","DEJ2000_BRm2","PMem_dd",])
    df["PMem_dd"] = np.arange(0,1.0,.1)
    df["PMem_holla"] = np.arange(0,1.0,.1)
    df["PMem_1"] = [.1] * 5 + [np.nan] * 5

    df = assign_PMem_mean(df)
    assert (df.PMem_mean.tolist()
            == pytest.approx([.1/3, .1, .5/3, .7/3, .3, .5, .6, .7, .8, .9]))

def test_extract_extensions():

    #---------------------------
    # Create fake data

    df = pd.DataFrame(columns=["RAJ2000_1","DEJ2000_1","PMem_1",
                               "RAJ2000","DEJ2000_holla","PMem_holla",
                               "RAJ2000_2","DEJ2000_BRm2","PMem_dd",])

    # Test output as extensions
    assert extract_extensions(df, prefix="RAJ2000") == ["1","2"]
    assert extract_extensions(df, prefix="DEJ2000") == ["1","holla","BRm2"]
    assert extract_extensions(df, prefix="PMem") == ["1","holla","dd"]

    # Test output as column names
    assert extract_extensions(df, prefix="PMem", retcolnames=True) == ['PMem_1', 'PMem_holla', 'PMem_dd']
    assert extract_extensions(df, prefix="RAJ2000", retcolnames=True) == ["RAJ2000_1","RAJ2000_2"]
    assert extract_extensions(df, prefix="DEJ2000", retcolnames=True) == ["DEJ2000_1","DEJ2000_holla","DEJ2000_BRm2"]

    # Must pass a string-convertible prefix
    with pytest.raises(TypeError):
        extract_extensions(df, prefix=None)
    with pytest.raises(TypeError):
        extract_extensions(df, prefix=1.44)

    assert set(extract_extensions(df, prefix="", retcolnames=True)) == set(df.columns.tolist()) - {"RAJ2000"}

def test_match_k2_epic_campaign():
    # set up fake data
    testdf = pd.DataFrame({"RAJ2000":[57.13268195367, 132.8329500, 59.],
                           "DEJ2000":[24.03288651412, 11.7834400, -25.],
                           "smh":["blo","blo","blo"]})

    # run function
    resdf = match_k2_epic_campaign(testdf)

    # test some values from the results
    assert resdf.Campaign[0] == 4.
    assert resdf.Campaign[1] == 5.
    assert resdf.Campaign[2] == 16.
    assert resdf.Campaign[3] == 18.
    assert np.isnan(resdf.Campaign[4])
    assert (resdf.smh.values == "blo").all()
    assert resdf.loc[4,:].DEJ2000 == -25.
    assert resdf.loc[4,:].RAJ2000 == 59.
    assert resdf.loc[0,:].EPIC == '211066477'
    assert resdf.loc[1,:].EPIC == '211409345'
    assert resdf.loc[2,:].EPIC == '211409345'
    assert resdf.loc[3,:].EPIC == '211409345'

def test_find_campaigns():
    # success modes are tested in match_k2_epic_campaign

    # test failure modes
    for t1 in [find_campaigns("blo ddd"), find_campaigns("67777 9888"),
               find_campaigns("0 nan")]:
        assert np.isnan(t1[0][0])
        assert np.isnan(t1[1])

def test_crossmatch():
    name = 'Tiny Cluster'
    sname ='tiny'
    a = pd.read_csv('catalog_matching/tests/exfiles/tiny_a.csv')
    b = pd.read_csv('catalog_matching/tests/exfiles/tiny_b.csv')
    d = pd.read_csv('catalog_matching/tests/exfiles/tiny_d.csv')#modified b
    keys = ['a','b']

    df = crossmatch(a, b, keys, name, sname, arcsec=5., plot=False)
    t = pd.read_csv('catalog_matching/tests/exfiles/crossmatch_result.csv')
    assert t.equals(df)

    df = crossmatch(a, b, keys, name, sname, arcsec=5., plot=False, union=True)
    t = pd.read_csv('catalog_matching/tests/exfiles/crossmatch_result_union.csv')
    assert t.equals(df)

    df = crossmatch(b, a, keys[::-1], name, sname, arcsec=5., plot=False, union=True)
    t = pd.read_csv('catalog_matching/tests/exfiles/crossmatch_result_union_reverse.csv')
    assert t.equals(df)

    df = crossmatch(a, d, ['a','d'], name, sname, arcsec=5., plot=False, union=True, bijective=True)
    t = pd.read_csv('catalog_matching/tests/exfiles/select_min_dist_union.csv')
    assert df.equals(t)

def test_select_min_dist():
    name = 'Tiny Cluster'
    sname ='tiny'
    a = pd.read_csv('catalog_matching/tests/exfiles/tiny_a.csv')
    b = pd.read_csv('catalog_matching/tests/exfiles/tiny_b.csv')
    c = pd.read_csv('catalog_matching/tests/exfiles/tiny_c.csv')#modified a
    d = pd.read_csv('catalog_matching/tests/exfiles/tiny_d.csv')#modified b

    df = crossmatch(c,b,['c','b'], name, sname, arcsec=5)
    df = select_min_dist(df,['c','b'])
    t = pd.read_csv('catalog_matching/tests/exfiles/select_min_dist.csv')
    assert df.equals(t)

    df = crossmatch(a,d,['a','d'], name, sname, arcsec=5, union=True)
    df = select_min_dist(df,['a','d'])
    t = pd.read_csv('catalog_matching/tests/exfiles/select_min_dist_union.csv')
    assert df.equals(t)


def test_post_k2_clean():
    """Test data cleaning and merging procedure
    on a tiny fake data set.
    """
    # Read in data file and process it:

    C = 'tiny'
    df = pd.read_csv('catalog_matching/tests/exfiles/select_min_dist_union_k2.csv')

    # This is created in the actual procedure as an intermediate step
    # This file is matched with the K2 catalog before further processing

    racols = list(filter(lambda k: ra in k, df.columns.values))
    deccols = list(filter(lambda k: dec in k, df.columns.values))
    df["{}_mean".format(ra)] = df[racols].mean(axis=1)
    df["{}_mean".format(dec)] = df[deccols].mean(axis=1)
    df[["{}_mean".format(ra),"{}_mean".format(dec)]].to_csv('catalog_matching/tests/exfiles/{}_radec.csv'.format(C),
                         index=False, header=False)

    # Read in data from K2 search:
    folder = 'catalog_matching/tests/exfiles/'
    k2 = pd.read_csv('{}{}_k2_search.txt'.format(folder,C), skiprows=[1], header=[0])

    # Call function
    k2mem = post_k2_clean(df, k2, "mean")
    k2mem = k2mem.dropna(how="all", axis=1).dropna(how="all", axis=0).T

    # ---------------------------------------------------------------------------
    # Validate some results:


    assert df.iloc[0].Pmem_a == k2mem['0'].Pmem_a
    for i in [0,1,2]:
        assert df.iloc[i].RAJ2000_a == pytest.approx(k2mem[str(i)].RAJ2000_mean, rel=1e-7)
        assert df.iloc[i].DEJ2000_a == pytest.approx(k2mem[str(i)].DEJ2000_mean, rel=1e-7)
    for i in [3,4,5,6]:
        assert k2mem[str(i)].RAJ2000_d == pytest.approx(k2mem[str(i)].RAJ2000_mean, rel=1e-7)
        assert k2mem[str(i)].DEJ2000_d == pytest.approx(k2mem[str(i)].DEJ2000_mean, rel=1e-7)
    assert df.shape[0]==k2mem.shape[1]

    assert k2mem.shape[1] == 7
