# Functions and classes related to the
# rotation-flaring-activity relation analysis live here.
#
#

import numpy as np
import pandas as pd

def rossby_wright(R, VK):
    '''Calculate the Rossby number
    based on Wright et al. (2011)
    Eq. (10) for V-K > 3.5'''
    if VK <= 3.5:
        raise ValueError
    else:
        return R / np.exp(-2.16 + 1.5 * VK - .13 * VK**2)


def prepare_rotation_activity_table(cluster, df, dff, rot=None, path=True,
                                    EDthreshold=300):
    '''Script that reads flare data, stellar parameters,
    and rotation information for a cluster and merges the tables.

    Parameters:
    ------------
    cluster : str
        cluster name
    df : DataFrame
        stellar parameters
    dff : DataFrame
        flare table
    rot : DataFrame or None
        rotation information
    path : bool
        if true will read rotation information  for rot
        from path
    EDthreshold : 300
        minimum ED to consider from flare table

    Return:
    -------
    DataFrame with rotation information with
    stellar parameters, and flaring information added.
    '''
    cl = {"praesepe": {"author":"rebull", "year":17, "color":"V-Ks", "Pname":"PPer", },
          "pleiades": {"author":"rebull", "year":16, "color":"(V-K)0", "Pname":"Prot"},
          "hyades": {"author":"douglas", "year":16, "color":"Mass", "Pname":"K2Per"},
          "generic" : {"author":"nemo", "year":0, "color":"V-K", "Pname":"Prot"}}

    obst = {5:1723, 13:1821.5, 16:1822.5, 18:1185.5, 4:1643.5}


    # init a cluster analysis
    CL = cl[cluster]
    author = CL["author"]
    year = CL["year"]
    color = CL["color"]
    rotation = CL["Pname"]


    if path == True:
        rot = pd.read_csv("ancillary/rotation/{}_rotation_{}{}.tsv".format(cluster, author, year),
                          delimiter="\t", skiprows=69, na_values=['       '])

    # attribute total observation time to stars and campaigns
    df["total_obs_time"] = df.Campaign.apply(lambda x: obst[x])

    #each LC has one observation time assigned:
    totalobstime = df[["EPIC","Campaign","total_obs_time"]].drop_duplicates()

    #now sum all observation times across campaigns for a given EPIC ID
    totalobstime = df.groupby("EPIC").total_obs_time.sum().reset_index()

    #rename and merge columns
    totalobstime = totalobstime.rename(index=str, columns={"total_obs_time":"tot_obs_time"})
    df = df.merge(totalobstime, how="left", on=["EPIC"])

    ## Cut on detection threshold or another minimum ED and calculate total energy released on each target

    # only use flare energies above a certain threshold
    dff = dff[dff.ed_rec > EDthreshold]

    # sum the flare energy in every target
    totalflare_energ = dff.groupby("EPIC").ed_rec.sum().reset_index()

    # rename and merge columns
    totalflare_energ = totalflare_energ.rename(index=str, columns={"ed_rec":"ed_rec_sum"})
    df = df.merge(totalflare_energ, how="left", on=["EPIC"])

    ## Assign a number of flares above a certain threshold observed for each target

    # count the flares per target from the flares table
    dff.ed_rec = dff.ed_rec.fillna(0)
    nflares = dff.groupby("EPIC").ed_rec.count().reset_index()
    # rename and merge columns
    nflares = nflares.rename(index=str, columns={"ed_rec":"nflares"})
    df = df.merge(nflares, how="left", on=["EPIC"])

    ## Add rotation information


    # add rotation periods after all

    # but first:
    # the hyades table needs more love

    if cluster == "hyades":
        rot["Mass"] = rot["Mass"].str.strip().replace('',np.nan)
        rot["K2Per"] = rot["K2Per"].str.strip().replace('',np.nan)
        rot["Mass"] = rot["Mass"].astype(float)
        rot["K2Per"] = rot["K2Per"].astype(float)

    # the praesepe table needs more love

    if cluster == "praesepe":
        rot["V-Ks"] = rot["V-Ks"].str.strip().replace('',np.nan)
        rot["V-Ks"] = rot["V-Ks"].astype(float)
        rot["PPer"] = rot["PPer"].astype(float)

    df2 = df.merge(rot, how="right", on=["EPIC"])
    df2["hasflare"] = 0
    df2.loc[df2.nflares > 0,"hasflare"] = 1

    # finally, the praesepe table need more love
    if cluster == "praesepe":
        df2 = df2[~df2.PPer.isnull()]
        df2 = df2[~df2["V-Ks"].isnull()]
        subset = ["EPIC", "V-Ks", "PPer", "nflares", "hasflare","tot_obs_time", "ed_rec_sum","Teff_median"]
        df2 = df2.drop_duplicates(subset=subset)

    # finally, the hyades table need more love
    if cluster == "hyades":
        df2 = df2[~df2.Mass.isnull()]
        df2 = df2[~df2.K2Per.isnull()]

    df2.nflares = df2.nflares.fillna(0)
    df2.ed_rec_sum = df2.ed_rec_sum.fillna(0)

    df2 = df2.rename(index=str, columns={color:"V-K", rotation:"Prot"})
    return df2