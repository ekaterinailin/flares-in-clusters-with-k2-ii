# Catalog matching functions are here.

from time import asctime as time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
from astroML.crossmatch import crossmatch_angular
from lightkurve import search_targetpixelfile

def pick_members_and_produce_k2_search_input(cross, sname, name,
                                             coords="GaiaC",
                                             PMem_mean_threshold=0.8):
    """Routine to pick the sufficiently certain members
    from a membership table and write out a RA/Dec table.

    Parameters:
    ------------
    df : DataFrame
        membership table with PMem_XXX columns where
        XXX is a catalog suffix
    sname, name : str, str
        file prefix and full cluster name
    coords : str
        RA/Dec value column extension, default Gaia DR2
    PMem_mean_threshold: float
        Threshold probability for likely membership
    """
    # Pick members

    print("There are {} candidate members in {}."
          .format(cross.shape[0], name))
    cross = assign_PMem_mean(cross)
    cross = cross[cross.PMem_mean > PMem_mean_threshold]
    print("There are {} members in {} with"
          " membership probability above {}."
          .format(cross.shape[0], name, PMem_mean_threshold))

    # Prepare and write out RA, Dec file to match with K2

    radeccols =["RAJ2000_{}".format(coords),"DEJ2000_{}".format(coords)]
    cross[radeccols].to_csv('catalog_matching/matched_catalogs/'
                            'membership_matches/radec/{}_radec.csv'
                            .format(sname),
                            index=False, header=False)
    return cross


def crossmatch_multiple_catalogs(l, name, sname, ra='RAJ2000',
                                 dec='DEJ2000', **kwargs):
    """Cross-match multiple catalogs, calculate
    relevant distances and preserve columns.

    Parameters:
    ------------
    l : dict
        Keys are catalog names (i.e. author names),
        values are the catalogs (DataFrames)
    name, sname : str, str
        cluster name and file prefix string
    kwargs : dict
        Keyword arguments to pass to func:crossmatch

    Return:
    DataFrame with all catalogs from input cross-matched.
    """

    keys = list(l.keys())

    # Initialize iteration with the first catalog on the list:

    cinterim = l[keys[0]]
    keyinterim = keys[0]

    # Iterate over the rest:

    for i, c in enumerate(keys[1:]):

        # Cross-match current table with the next one on the list:

        cinterim = crossmatch(cinterim, l[c], [keyinterim, c], name, sname, **kwargs)

        # Re-define columns

        racols = list(filter(lambda k: ra in k, cinterim.columns.values))
        deccols = list(filter(lambda k: dec in k, cinterim.columns.values))
        keyinterim = '_'.join(keys[:i+2])
        newra = '{}_{}'.format(ra, keyinterim)
        newdec = '{}_{}'.format(dec, keyinterim)
        cinterim[newra] = cinterim[racols].mean(axis=1)
        cinterim[newdec] = cinterim[deccols].mean(axis=1)

    return cinterim

def extract_extensions(df, prefix='DEJ2000', retcolnames=False):
    '''
    Return the extensions to existing position columns with a given
    prefix. Passing an empty string as prefix will return all suffixes
    to all column names that have one.

    Paramaters:
    --------------
    df : DataFrame

    prefix : 'DEJ2000' or 'RAJ2000' or 'PMem'
        Data extracted from multiple membership
        catalogs.

    Returns:
    --------
    List of extension strings.
    '''
    cols = df.columns.values.astype(str)
    id1 = np.flatnonzero(np.core.defchararray.find(cols,prefix)!=-1)
    id2 = np.where(np.array([i.count('_') for i in cols[id1]]) == 1)
    if retcolnames==False:
        return [x.split('_')[-1] for x in cols[id1][id2]]
    elif retcolnames==True:
        return cols[id1][id2].tolist()

def assign_PMem_mean(df, prefix="PMem"):
    cols = extract_extensions(df, prefix=prefix, retcolnames=True)
    df['PMem_mean'] = df[cols].mean(axis=1)
    return df

def merge_k2_info_and_protocol(cross, sname, col,
                               folder=('catalog_matching/matched_'
                                       'catalogs/k2_matches_to_union/')):
    """Read in the results from K2 catalog query.
    Merge with membership table.
    Write out membership numbers to file.

    Parameters:
    -----------
    cross : DataFrame
        membership table
    col : str
        suffix for RA and Dec columns
        to pass to post_k2_clean
    sname : string
        filename prefix for cluster
    """
    k2 = pd.read_csv('{}{}_k2_search.txt'.format(folder, sname), skiprows=[1], header=[0])

    k2mem = post_k2_clean(cross, k2, col)
    k2mem.to_csv('{}{}_k2_mem.csv'.format(folder,sname), index=False)
    nums = [time(),
            cross.shape[0],
            k2mem[~k2mem.EPIC.isnull()].shape[0],
            k2mem[k2mem["Target Type"]=="SC"].shape[0],
            k2mem[k2mem["Target Type"]=="LC"].shape[0],
            sname,]

    with open("{}members.txt".format(folder), "a") as f:
        f.write(",".join(["date","nlikelymembers","LCs","SLCs","LLCs","cluster"]))
        f.write("\n")
        f.write(",".join([str(x) for x in nums]))
        f.write("\n")

def match_k2_epic_campaign(df, ra="RAJ2000", dec="DEJ2000"):
    '''Merge a DataFrame with corresponding EPIC ID and
    K2 campaign in which the target was observed, producing
    at least a TargetPixelFile.

    Parameters:
    ------------
    df: DataFrame
        table with RA and Dec
    ra : str
        column name for RA
    dec : str
        column name for Dec

    Return:
    --------
    DataFrame with extra columns "EPIC" and "Campaign"
    '''
    radecs = df[[ra,dec]].apply(lambda x: "{} {}".format(x[ra],x[dec]),axis=1)
    campaigns = radecs.apply(find_campaigns)
    d1 = pd.DataFrame(campaigns.tolist(), index=campaigns.index)
    d2 = d1.join(df)
    d2 = d2.rename(index=str, columns={1:"EPIC"})
    s = d2.apply(lambda x: pd.Series(x[0]),axis=1).stack().reset_index(level=1, drop=True).astype(int)
    s.name = "Campaign"
    d2 = d2.drop(0, axis=1).join(s)
    d2 = d2.reset_index().drop("index",axis=1)

    return d2

def find_campaigns(s):
    '''Find campaigns and EPIC IDs
    for a given RA and Dec pair
    using lightkurve query function.

    Parameters:
    ------------
    s : str
        "RA Dec" in deg

    Return:
    -------
    list, int : list of campaigns, EPIC ID
    '''
    res = search_targetpixelfile(s, radius=2)
    pdf = res.table.to_pandas()

    if pdf.shape[0] > 0:
        pdf['Campaign'] = pdf.obs_id.str[-5:-3].astype(int)
        pdf['EPIC'] = pdf.target_name.str[4:]
        try:
            assert (pdf.EPIC.values == pdf.EPIC.iloc[0]).all()
        except:
            print("Something went wrong here\n")
            print(pdf.EPIC.values)
            return [np.nan], np.nan

        return pdf.Campaign.tolist(), pdf.EPIC.iloc[0]
    else:
        return [np.nan], np.nan


def crossmatch(a, b, keys, name, sname, arcsec=5., bijective =False,
               plot=False, union=False):

    '''
    Match targets from catalog c2 into reference catalog c1
    in a cone with a maximum search radius.

    Parameters:
    ------------
    a : DataFrame
    b : DataFrame
    keys : list with length 2
        Contains extensions into RA, Dec, and PMem in a and b,
        respectively.
    arcsec : 5. or float
        Maximum search radius in arcsec.
    plot : False or bool
        Export a plot of matching distances.
    bijective : True or bool
        Force 1-on-1 matching

    Return:
    --------
    DataFrame of length c1 with c2 matched if possible. Otherwise
    the row is filled with the last entry in c2 and the distance is
    filled with inf.
    '''
    c1 = copy(a).astype(float)
    c2 = copy(b).astype(float)

    # get imaging data
    imX = np.empty((c1.shape[0], 2), dtype=np.float64)
    imX[:, 0] = c1['RAJ2000_{}'.format(keys[0])].values
    imX[:, 1] = c1['DEJ2000_{}'.format(keys[0])].values

    # get standard stars
    stX = np.empty((c2.shape[0], 2), dtype=np.float64)
    stX[:, 0] = c2['RAJ2000_{}'.format(keys[1])].values
    stX[:, 1] = c2['DEJ2000_{}'.format(keys[1])].values

    # crossmatch catalogs
    max_radius = arcsec / 3600
    dist, ind = crossmatch_angular(imX, stX, max_radius)
    match = ~np.isinf(dist)#indices into imX

    #matching distances
    dist_match = dist[match]
    dist_match *= 3600

    #merge data frames
    ind = list(filter(lambda x: x != c2.shape[0], ind))#throw out fill values
    for colname in c2.columns.values:
        c1[colname] = np.nan
        c1.iloc[match, c1.columns.get_loc(colname)] = c2.iloc[ind, c2.columns.get_loc(colname)].tolist()

    #add distances to results
    distcol = 'dist_{}_{}'.format(keys[0], keys[1])
    c1[distcol] = np.nan
    c1.iloc[match, c1.columns.get_loc(distcol)] = dist_match

    #Force 1-on-1 matching if requested
    if bijective == True:
        c1 = select_min_dist(c1, keys)

    #append the rest of c2 if required
    if union == True:
        allid = list(range(c2.shape[0]))
        ids, counts = np.unique(allid+ind, return_counts=True)
        unmatched_ids = ids[np.where(counts==1)]
        c1 = c1.append(c2.iloc[unmatched_ids,:],ignore_index=True, sort=False)

    #plotting optional
    if plot == True:
        plot_radius_of_match(dist_match,name, sname,
                             keys, arcsec, match, imX, stX)
    return c1

def select_min_dist(dataframe, keys):
    '''
    Makes the match unique. keys[1] targets may fit keys[0] targets
    multiple times. Filter only those that are closest fits.

    Parameters:
    ------------
    dataframe : DataFrame
        Contains the crossmatched dataframe as returned by
        crossmatch().
    keys : list of strings
        The two extensions to column names are stored in the
        list of length 2, first the reference catalog, second
        the image catalog.

    Return:
    -------
    DataFrame with non-optimal fits filtered out.
    '''
    df = copy(dataframe)
    nonnan = df.groupby(['RAJ2000_{}'.format(keys[1]),
                         'DEJ2000_{}'.format(keys[1])])['dist_{}_{}'.format(keys[0], keys[1])].idxmin().dropna()
    nonnan = nonnan.astype(int).values
    nonnan = nonnan[~np.isnan(nonnan)]
    nan = df[df['RAJ2000_{}'.format(keys[1])].isna()].index.values
    distnan = df[~df['RAJ2000_{}'.format(keys[1])].isnull() & df['dist_{}_{}'.format(keys[0], keys[1])].isnull()].index.values
    uniques = list(np.concatenate([nan,nonnan,distnan]))
    allid = list(range(df.shape[0]))
    ids, counts = np.unique(allid+uniques, return_counts=True)
    unmatched_ids = ids[np.where(counts==1)]
    unmatched_cols = list(filter(lambda k: '_{}'.format(keys[1]) in k, df.columns.values))
    df.loc[unmatched_ids,unmatched_cols] = np.nan
    return df

def plot_radius_of_match(dist_match,name, sname, keys, arcsec, match, imX, stX):
    '''Plot the output from crossmatch()'''
    ax = plt.axes()
    ax.hist(dist_match,
         histtype='stepfilled', ec='k', fc='#AAAAAA')
    ax.set_xlabel('radius of match (arcsec)')
    ax.set_ylabel('N(r, r+dr)')
    ax.text(0.95, 0.95,
            "Total unique objects: %i\nNumber with match: %i" % (imX.shape[0] + stX.shape[0] - np.sum(match),
                                                          np.sum(match)),
            ha='right', va='top', transform=ax.transAxes)

    ax.set_title('{}: {} matches to {}'.format(name, keys[1], keys[0]))
    plt.savefig('plots/matched_catalogs/{}_{}_{}_maxradius{}arcsec.png'.format(sname, keys[0], keys[1],arcsec), dpi=300)
    plt.close()
    return

def get_catalogs(tc, pathlist):
    '''
    Read in all catalogs into a dictionary with keys
    that are the extensions to the respective columns.

    Parameters:
    ------------
    tc : DataFrame
        Characteristica of catalog tables
    pathlist : list
        Paths to the catalogs whose characteristics are
        noted in tc.
    '''
    cats = {}
    for i, row in tc.iterrows():
        for p in pathlist:
            s = '{}/{}.{}'.format(row.folder, row['name'], row.ext)
            if p.find(s) != -1:
                if row.nskip != 0:
                    r = int(row.nskip)
                else:
                    r = None
                df = pd.read_table(p,skiprows=r, usecols=row.usecols, delimiter=row['sep'])
                cats[row.folder[9:]] = df
    return cats

def rename_columns(df, ext, cols):
    '''
    Rename columns in dataframe with useful extensions.

    Parameters:
    ------------
    df : DataFrame
        Catalog without extensions.
    ext : string
        Extension to use for this catalog.
    cols : list of strings
        column name bases to replace old ones.

    Return:
    ---------
    Catalog with useful extensions for later matching.
    '''
    new_colnames = ['{}_{}'.format(x,ext) for x in cols]
    return df.rename(index=str, columns=dict(zip(df.columns.values,new_colnames)))

def post_k2_clean(df, k2, col):
    '''Clean up the k2 matched catalog and merge back
    in to the cluster membership catalog using the crossmatch
    function.

    Parameters:
    ------------
    df : DataFrame
        matched cluster membership catalog
    k2 : DataFrame
        K2 data that the input of RA and Dec from
        df returns.
    col : str
        suffix for RA and Dec columns in df

    Return:
    -----------
    Dataframe with k2 and df merged on RA and Dec.

    '''
    ra = 'RAJ2000'
    dec = 'DEJ2000'
    rak2 = "RAJ2000_K2"
    dek2 = "DEJ2000_K2"
    k2[rak2] = k2["RA (J2000)"]
    k2[dek2] = k2["Dec (J2000)"]

    k2mem = crossmatch(df, k2[[rak2,dek2,"Campaign","K2 ID"]], [col,"K2"], col, col, arcsec=3, bijective =False,
               plot=False, union=False)

    #rename columns
   # k2mem.columns = k2mem.columns.str.replace("_"+col,"")
    k2.columns = k2.columns.str.replace("_K2","")

    #merge input and crossmatch output
    k2mem = k2mem.merge(k2, how="outer", on=["Campaign", "K2 ID"])
    #rename to final form
    k2mem = k2mem.rename(index=str, columns={'K2 ID':'EPIC',
                                             'RA (J2000)': 'RAJ2000_K2',
                                             'Dec (J2000)': 'DEJ2000_K2',
                                             "Ang Sep (')" : 'dist_mean_K2_arcmin',
                                             'RAJ2000' : 'RAJ2000_mean',
                                             'DEJ2000' : 'DEJ2000_mean',})
    df = df.drop([ra+"_"+col,dec+"_"+col], axis=1)

    # Remove duplicates in rows, then in columns
    k2mem = k2mem.drop_duplicates().loc[:,~k2mem.columns.duplicated()]

    keepcolumns = ['EPIC', 'Dataset Name', 'Campaign', 'Data Release',
       'RAJ2000_K2', 'DEJ2000_K2', 'Target Type',
       'RA PM', 'RA PM Err', 'Dec PM', 'Dec PM Err', 'Plx', 'Plx Err',
       'U Mag', 'U Mag Err', 'B Mag', 'B Mag Err', 'V Mag', 'V Mag Err',
       'G Mag', 'G Mag Err', 'R Mag', 'R Mag Err', 'I Mag', 'I Mag Err',
       'Z Mag', 'Z Mag Err', 'J Mag', 'J Mag Err', 'H Mag', 'H Mag Err',
       'K Mag', 'K Mag Err', 'KEP Mag', 'Kep Flag', 'Hip ID', 'Tyc ID',
       'SDSS ID', 'UCAC ID', '2MASS ID', '2MASS Flag',
       'crowding', 'contamination', 'flux fraction', 'cdpp3', 'cdpp6',
       'cdpp12', 'Module', 'Output', 'Channel', 'Nearest Neighbor',
       'Nomad ID', 'dist_mean_K2_arcmin',  'RAJ2000_mean',  'DEJ2000_mean'] + list(df.columns.values)
    k2mem = k2mem[keepcolumns]
    return k2mem

#NOT USED A LOT, may deprecate soon, tried using in Reino table treatment

def match_gaia(name, k2, mem_prob=.8, mem='CG'):
    '''
    Match EPIC IDs with long cadence LCs to membership
    lists in Gaia DR2 obtained from Cantat-Gaudin (2018).
    Saves the matched catalog into the original folder.

    Parameters
    -----------
    name :  str
        cluster name in file name
    k2 : DataFrame
        pandas dataframe with EPIC and DR2 IDs
    mem_prob: .8 or float
        minimum membership probability according to
        Cantat-Gaudin (2018)
    mem : 'CG' or 'GaiaC'
        membership source

    '''
    if mem=='CG':
        sr = list(range(76))+[77, 78]
        source = 'Cantat-Gaudin (2018)'
    elif mem=='GaiaC':
        sr = list(range(52))+[53, 54]
        source = 'Gaia Collaboration (2018)'
    cg_2018 = pd.read_csv('{}_clusters/{}.tsv'.format(mem,name),
                      skiprows=sr, delimiter='\t')
    if mem=='CG':
        cg_prob = cg_2018[cg_2018.PMemb >= mem_prob]
    else:
        cg_prob = cg_2018
    cg_prob_match = cg_prob.merge(k2, left_on='Source',right_on='source_id', how='left')
    print('{} has {} long cadence LCs in the list of'
          ' {} members from {} with '
          'membership probability >= {}'.format(name, cg_prob_match.shape, cg_prob.shape, source, mem_prob))
    cg_prob_match.to_csv('{}_clusters/{}_k2_matched_prob_{}.csv'.format(mem,name, mem_prob))
    return