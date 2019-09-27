import warnings

import pandas as pd
import numpy as np

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Button
from bokeh.models.callbacks import CustomJS
from bokeh.models.layouts import Column
from bokeh.io import output_file, show

from astropy import units as u
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.coordinates import SkyCoord, Distance

from dustmaps.bayestar import BayestarWebQuery

from gaia.gaia import calculate_distance_from_parallax

from opencluster import PACKAGEDIR

def interactive_CMD(specs,cid1='g_SDSS',cid2='i_SDSS'):
    '''
    Simplistic tool to create an interactive
    bokeh plot where outliers can be marked and saved in
    '/home/ekaterina/Documents/appaloosa/stars_shortlist/share/temp'
    '''
    # Create some random data and put it into a ColumnDataSource
    s = specs.set_index('EPIC')
    s = s[[cid1, cid2, 'e_'+cid1,'e_'+cid2]].dropna(how='any')
    x = list(s[cid1]-s[cid2])
    y = list(s[cid2])
    size = list(np.sqrt(s['e_'+cid1]**2+s['e_'+cid2]**2)*100.)
    z = list(s.index.values)
    source_data = ColumnDataSource(data=dict(x=x, y=y,desc=z))

    # Create a button that saves the coordinates of selected data points to a file
    savebutton = Button(label="Save", button_type="success")
    savebutton.callback = CustomJS(args=dict(source_data=source_data), code="""
            var inds = source_data.selected['1d'].indices;
            var data = source_data.data;
            var out = "";
            for (i = 0; i < inds.length; i++) {
                out += data['desc'][inds[i]] + " ";
            }
            var file = new Blob([out], {type: 'text/plain'});
            var elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(file);
            elem.download = 'selected-data.txt';
            document.body.appendChild(elem);
            elem.click();
            document.body.removeChild(elem);
            """)

    # Plot the data and save the html file
    p = figure(plot_width=800, plot_height=400,
               #y_range=(20,7),
               tools="lasso_select, reset, hover",)
    p.circle(x='x', y='y',  source=source_data, fill_alpha=0.8)#add size='desc' to visualize uncertainties on mag
    p.xaxis.axis_label = '{}-{}'.format(cid1,cid2)
    p.yaxis.axis_label = cid1
    plot = Column(p, savebutton)
    output_file("test.html")
    show(plot)
    return


def correct_for_extinction(df, dustmap='bayestar2019', thresh=10):
    '''Correct for dust extinction using 3D dustmaps.

    Parameters:
    --------------
    df : DataFrame
        stellar attributes including Gaia, 2MASS, PS1
    dustmap : str
        Which dustmap to use. Check dustmap doc for
        more dustmaps. Default: Bayestar2019

    '''
    dft = df.loc[:,df.columns.str.contains("Gaia")].rename(columns = lambda x : str(x)[:-5])
    dft = calculate_distance_from_parallax(dft, check_GoF=False)

    # Fill in the median cluster distance for targets without Gaia parallax:
    dft.loc[dft['distance'].isnull(),"distance"] = dft['distance'].median()

    d = np.abs(dft["distance"].values) * u.pc
    ra = dft.ra.values * u.deg
    dec = dft.dec.values * u.deg
    # validate distance, RA, and Dec
    assert len(d) == dft.shape[0]
    assert len(ra) == dft.shape[0]
    assert len(dec) == dft.shape[0]

    ext, ext_err, flags = query_dustmap_w_percentiles(d, ra, dec, dustmap=dustmap)

    # Add results from dustmap query to stars:

    df['{}_ext'.format(dustmap)] = ext
    df['{}_ext_err'.format(dustmap)] = ext_err
    df['{}_flags'.format(dustmap)] = flags
    #print(list(zip(df['{}_flags'.format(dustmap)].tolist(),d)))
    # Apply flags
    df.loc[df['{}_flags'.format(dustmap)] > thresh, #i.e. star too close or too far, or algorithm did not converge
           ['{}_ext_err'.format(dustmap),
            '{}_ext'.format(dustmap)] ] = np.nan

    # Convert to ext values for 2MASS and PanSTARRS

    # This table is Table 1 in http://argonaut.skymaps.info/usage

    conversion = dict(zip(['g_PS1', 'r_PS1', 'i_PS1', 'z_PS1',
                           'y_PS1', 'J_2MASS', 'H_2MASS', 'K_2MASS'],
                          [3.518, 2.617, 1.971, 1.549,
                           1.263, 0.7927, 0.4690, 0.3026]))
    for key in conversion:
        # This nans everything that does not have extinction values given:
        # m_band_intrinsic = m_band_observed - extinction_value * conversion_factor
        df[key] -= (conversion[key] * df['{}_ext'.format(dustmap)])
        # Propagate uncertainty in quadrature on extinction value to photometry:
        df['e_{}'.format(key)] = np.sqrt(df['e_{}'.format(key)]**2 +
                                         (conversion[key] * df['{}_ext_err'.format(dustmap)])**2)
    return df


def query_dustmap_w_percentiles(d, ra, dec, dustmap='bayestar2019'):
    '''Query 3D dustmaps. Cite Green (2018) if
    you use them. This func queries the Bayestar2019
    dust map remotely in default mode.
    (The web interface takes the same arguments
    as the local interface.)

    Parameters:
    -----------
    d : 1d array with astropy units
        distance along the line of sight
    ra : 1d array with astropy units
        RA in ICRS
    dec: 1d array with astropy units
        Declination ICRS

    Return:
    ---------
    ext : 1-d array of floats
        Extinction value from given dustmap
    ext_err : 1-d array of floats
        uncertainty on ext
    flags : 1-d array of ints
        0 : no extinction given
        1 : uncertainties are symmetric below 0.05
        2 : uncertainties are symmetric below 0.10
    '''
    with warnings.catch_warnings():

        # silence warnings from astropy
        warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)

        # Query dustmaps frm Bayestar2019:
        coords = SkyCoord(ra, dec, distance=d, frame='icrs')
        q = BayestarWebQuery(version=dustmap)
        E, quality = q(coords, mode='percentile', pct=[36,50,84], return_flags=True)

    # Make output numpy-ish
    E = np.array(E)
    Efl = np.nan_to_num(E)
    #print(quality)
    # Flag outputs that either have asymmetric uncertainties, or are NaNs;
    flags = (np.abs(2 * Efl[:,1]-Efl[:,0]-Efl[:,2]) < 0.1).astype(int) * 1 #acceptable uncertainty
    flags = (np.abs(2 * Efl[:,1]-Efl[:,0]-Efl[:,2]) < 0.05).astype(int) * 2 #low uncertainty
    flags[np.isnan(E[:,1])] += 4 #flags 1 and 2 failed to compute anything
    flags[~quality["reliable_dist"]] += 8 #too close or too far target
    flags[~quality["converged"]] += 16 #Algorithm did not converge
    # define extinction value and uncertainty as mean of percentiles:
    ext = E[:,1]
    ext_err = (E[:,2]-E[:,0]) / 2.
    return ext, ext_err, flags


def K_2MASS_to_KBB(K_2MASS, J_K_BB, e_K_2MASS=0., e_J_K_BB=0.):
    '''
    Kwargs are not yet tested!
    Parameters:
    -----------
    K_2_MASS: Series of floats
        K_s 2MASS magnitude
    J_K_BB : Series of floats
        J-K Bessel and Brett
    e_K_2MASS : float
        uncertainty on K_s 2MASS magnitude
    e_JKBB : float
        uncertainty on J-K Bessel and Brett color

    Return:
    --------
    K_BB : K in Bessel and Brett system
    e_K_BB : uncertainty on K_BB
    '''
    e_K_BB = e_K_2MASS**2 + 0.007**2 + J_K_BB**2 * 0.005**2 + 0.001**2 * e_J_K_BB**2
    K_BB = K_2MASS + 0.039 - 0.001 * J_K_BB
    return K_BB, e_K_BB

def color_2MASS_to_BB(df):
    '''
    Convert 2MASS J-K, H-K, and J-H to Bessel and Brett 1988
    J-K, H-K, J-H, and K using relations from Carpenter2001.
    '''
    def e_A_B(row, A, B, e_A, e_B):
        '''Calculate uncertainty on BB color.'''
        e_AB = np.sqrt(e_A**2 + e_B**2)
        return np.sqrt(e_AB**2 +
                       ((A-B) * row.a1_sigma / row.a1)**2 +
                       row.a0_sigma**2) / row.a1

    conv = pd.read_csv('{}/static/BB_to_TWOMASS.csv'.format(PACKAGEDIR))

    for index, row in conv.iterrows():
        b1, b2 = list(row.color)
        dfb1, dfb2 = b1 + '_2MASS', b2 + '_2MASS'
        e_dfb1, e_dfb2 = 'e_' + b1 + '_2MASS', 'e_' + b2 + '_2MASS'
        f = lambda x: (x - row.a0) / row.a1
        df['{}_{}_BB'.format(b1,b2)] = (df[dfb1] - df[dfb2]).apply(f)
        df['e_{}_{}_BB'.format(b1,b2)] = e_A_B(row, df[dfb1], df[dfb2], df[e_dfb1], df[e_dfb2])
    df['K_BB'], df['e_K_BB'] = K_2MASS_to_KBB(df.K_2MASS, df.J_K_BB)
    return df

def color_BB_to_Johnson(df):
    '''
    Convert Bessel/Brett J-K, H-K, and J-H to Johnson
    J-K, H-K, J-H, and K using relations from
    Bessel and Brett 1988.
    '''
    conv = pd.read_csv('{}/static/BB_to_Johnson.csv'.format(PACKAGEDIR))
    for index, row in conv.iterrows():

        b1, b2 = list(row.color) # "JH" -> ["J", "H"]

        bbcol = '{}_{}_BB'.format(b1,b2)
        e_bbcol = 'e_{}_{}_BB'.format(b1,b2)

        f = lambda x: row.a1 * x + row.a0
        df['{}_{}_Johnson'.format(b1,b2)] = df[bbcol].apply(f)

        e_f = lambda e_x: row.a1 * e_x
        df['e_{}_{}_Johnson'.format(b1,b2)] = df[e_bbcol].apply(e_f)

        df['K_Johnson'] = df['K_BB']
        df['e_K_Johnson'] = df['e_K_BB']

    df['J_Johnson'] = df.J_K_Johnson+df.K_Johnson
    df['H_Johnson'] = df.H_K_Johnson+df.K_Johnson

    df['e_J_Johnson'] = np.sqrt(df.e_J_K_Johnson**2 + df.e_K_Johnson**2)
    df['e_H_Johnson'] = np.sqrt(df.e_H_K_Johnson**2 + df.e_K_Johnson**2)

    return df

def color_2MASS_to_Johnson(df):
    '''
    Wrap up the full conversion needed to use color-temperature
    relations in Boyajian+2013.

    Individual functions are tested.
    '''
    df = color_2MASS_to_BB(df)
    df = color_BB_to_Johnson(df)
    return df

def Teff_Boyajian(df):
    '''
    Apply color-temperature relations from Boyajian et al. (2013).
    Drop values if the fall out of range.
    '''
    df = color_2MASS_to_Johnson(df)
    colorcols = dict(zip(list('JHKgrizyBV'),
                     ['J_Johnson','H_Johnson','K_Johnson',
                      'g_SDSS','r_SDSS','i_SDSS',
                      'z_SDSS','y_SDSS','B_K2','V_K2']))
    bm = pd.read_csv('{}/static/boyajian13_optimized_for_Teff_calculation.csv'.format(PACKAGEDIR), skiprows=15)

    for index, row in bm.iterrows():

        b1, b2 = list(row.color)
        dfb1, dfb2 = colorcols[b1], colorcols[b2]
        x = df[dfb1] - df[dfb2]
        feh = df['FeH']
        x_err = np.sqrt(df["e_" + dfb1]**2 + df["e_" + dfb2]**2)

        #apply range for each CTR:
        x[((x < row.mag_min) | (x > row.mag_max))] = np.nan

        # apply [Fe/H] range
        # Boyajian's values are centered on solar with a spread about (-0.25; 0.25)
        x[((feh < -.25) | (feh > .25))] = np.nan

        teff = row.a0 + row.a1 * x + row.a2 * x**2 + row.a3 * x**3
        Teffdes = 'Teff_{}_{}_Boy'.format(b1, b2)
        df[Teffdes] = teff
        df['e_' + Teffdes] = np.sqrt((row.a1 + row.a2 * 2 * x + row.a3 * 3 * x**2)**2 * x_err**2 +
                                                            (row.sigma_percent / 100. * teff)**2)
        print("Number of {} ".format(Teffdes), df[Teffdes].count())
    return df


def Teff_Apsis(df, av_uncertainty=175., trange=(4099,6750)):
    """Use Apsis (Andrae et al. 2018) Teff
    values using the median uncertainty
    in a region with reliable estimates
    according to the paper.

    Parameters:
    ------------
    df : DataFrame

    Return:
    ------
    DataFrame

    """
    # use Teffs from Apsis catalog
    cols = ["Teff_Apsis", "e_Teff_Apsis"]
    df[cols[0]] = df["teff_val_Gaia"]
    df[cols[1]] = av_uncertainty

    # remove all value that fall outside of reliabel regions
    condition = (df[cols[0]] > trange[0]) & (df[cols[0]] < trange[1]) & (df.BPRP_Gaia < 2.)
    df.loc[~condition, cols] = np.nan

    return df


def mann_formula(x, p):
    '''
    Parameters:
    ------------
    x : 4-tuple of arrays of floats
        0: colors,
        1: J-H or Fe/H optional,
        2: uncertainty on colors,
        3: uncertainty on J-H or Fe/H
    p : Series
        coefficients from Mann+2015 Table 2
    '''
    sig_phot_squared = (p.b + 2 * x[0] * p.c + 3 * x[0]**2 * p.d + 4 * x[0]**3 * p.e)**2 * x[2]**2 * p.Teff_factor**2
    if p.formula == 4:
        teff = p.a + x[0]*p.b + x[0]**2 * p.c + x[0]**3 * p.d + x[0]**4 * p.e

    elif p.formula == 6:

        teff = p.a + x[0]*p.b + x[0]**2 * p.c + x[0]**3 * p.d + x[0]**4 * p.e + x[1]*p.JH_FeH
        sig_phot_squared += p.JH_FeH**2 * x[3]**2 * p.Teff_factor**2
    elif p.formula == 7:
        teff = p.a + x[0]*p.b + x[0]**2 * p.c + x[0]**3 * p.d + x[0]**4 * p.e + x[1]*p.JH_FeH + x[1]**2 * p.JH2
        sig_phot_squared += (p.JH_FeH + p.JH2 * 2)**2 * x[3]**2 * p.Teff_factor**2

    # if no uncertainties on colors are given, use 500 K to blow up uncertainties
    sig_phot_squared[np.isnan(sig_phot_squared)] = 250000
    sig_teff = np.sqrt(sig_phot_squared + p.add_in_quadrature_K**2 + p.sigma**2)
    return teff * p.Teff_factor, sig_teff

def Teff_Mann(df):
    '''
    Apply color-temperature relations from the Erratum to Mann+2015.
    Gaia uncertainties on extinction NOT included.
    '''
    mm = pd.read_csv('{}/static/mann15_optimized_for_Teff_calculation.csv'.format(PACKAGEDIR))
    colorcols = dict(zip(['r','z','J','BP','RP', 'E_BP_RP'],
                     ['r_SDSS','z_SDSS', 'J_2MASS',
                      'BP_Gaia','RP_Gaia','e_bp_min_rp_val_Gaia']))
                       # the last one is for extinction
    for index, row in mm.iterrows():

        #if Gaia apply corrected BP-RP
        if row.c1=="BP":
            color = df["BPRP_Gaia_corr"]
            err = df["e_BPRP_Gaia_corr"]

        else:
            color = df[colorcols[row.c1]] - df[colorcols[row.c2]]
            err = np.sqrt(df["e_" + colorcols[row.c1]]**2 + df["e_" + colorcols[row.c2]]**2)

        if row.extra == 'isJH':
            extra = df.J_2MASS - df.H_2MASS
            extra_err = np.sqrt(df.e_J_2MASS**2 + df.e_H_2MASS**2)
        elif row.extra == 'FeH':
            extra = df.FeH
            extra_err = df.e_FeH
        else:
            extra = pd.Series()
            extra_err = pd.Series()
        X = (color, extra, err, extra_err)
        Teffdes = 'Teff_Mann_{}_{}_{}'.format(row.c1,row.c2,row.extra)
        df[Teffdes] = mann_formula(X, row)[0]
        df['e_'+ Teffdes] = mann_formula(X, row)[1]
        df.loc[df['e_'+Teffdes].isnull(), Teffdes] = np.nan
        df.loc[df[Teffdes].isnull(), 'e_'+Teffdes] = np.nan
        print("Number of {} ".format(Teffdes), df[Teffdes].count())
    return df

def prioritize_Teff_Mann(df):
    '''
    As stated in the paper, we prefer Teff from
    colors+[Fe/H]
    to colors+(J-H)
    to colors without any extra information.
    The less informative values are dropped from the table.
    '''
    colors = ['r_z', 'BP_RP','r_J']
    for col in colors:
        rz = ((df.columns.str.contains('Teff_Mann')) & (df.columns.str.contains(col)))
        feh = (df.columns.str.contains('FeH'))
        nan = (df.columns.str.contains('nan'))
        isjh = df.columns.str.contains('isJH')
        #if Teff is available from FeH, drop the other two, because Teff depend on metallicity
        ids = ~df.loc[:,feh & rz].isnull().iloc[:,0].values
        df.loc[ids, (rz & (isjh | nan))] = np.nan
        #if Teff is not available from FeH, but from J-H, drop the pure color-derived, because J-H is a proxy to FeH
        ids = ~df.loc[:,isjh & rz].isnull().iloc[:,0].values
        df.loc[ids, (rz & nan)] = np.nan
    return df

def cut_Teff_Mann_on_Teff(df):
    '''
    Apply the Teff cuts that are given in Mann+2015.
    '''
    Tmin = 2700.
    Tmax = 4100.
    teffmann = df.loc[:,df.columns.str.startswith('Teff_Mann')].columns.values
    e_teffmann = df.loc[:,df.columns.str.startswith('e_Teff_Mann')].columns.values
    for col, e_col in zip(teffmann, e_teffmann):
        #print("misfits Teff: ", df.loc[((df.loc[:, col]>Tmax) | (df.loc[:, col]<Tmin)), col].shape[0])
        df.loc[((df.loc[:, col]>Tmax) | (df.loc[:, col]<Tmin)), [col, e_col]] = np.nan

    return df


def cut_Teff_Mann_on_Mks(df):
    '''
    Apply the Mks cuts that are given in Mann+2015.
    '''
    # Find distances:
    dft = df.loc[:,df.columns.str.contains("Gaia")].rename(columns = lambda x : str(x)[:-5])
    dft = calculate_distance_from_parallax(dft, check_GoF=False)

    #Convert to distance modulus:
    dft["distmod"] = Distance(np.abs(dft.distance.values), u.pc).distmod

    #Calculate Mks
    dft["Mks"] = df.K_2MASS - dft.distmod

    # Do the cuts:
    teffmann = df.loc[:,df.columns.str.startswith('Teff_Mann')].columns.values
    e_teffmann = df.loc[:,df.columns.str.startswith('e_Teff_Mann')].columns.values
    for col, e_col in zip(teffmann, e_teffmann):
        #print("misfits mks: ",
        #      df.loc[(dft["Mks"] > 9.8), col].shape[0],
        #      df.loc[(dft["Mks"] < 4.6),col].shape[0],
        #      df.loc[dft.Mks.isnull(), col].shape[0])
        df.loc[((dft["Mks"] > 9.8) | (dft["Mks"] < 4.6) | dft.Mks.isnull()), [col, e_col]] = np.nan
#

    return df

def calculate_weighted_mean(df, prefix="Teff_", colval = ("Teff_median", "Teff_std")):
    """Calculate the weighted Teff mean and uncertainty.

    Parameters:
    -----------
    df : DataFrame
        table with all Teff and e_Teff columns that are
        relevant

    Return:
    --------
    DataFrame with Teff_median and Teff_std
    """
    teff = df.loc[:,df.columns.str.startswith(prefix)].columns.values

    for i,t in enumerate(teff):
        df["__xp" + str(i)] = df[t] / df["e_"+t]**2 # x * p
        df["__p" + str(i)] = 1./df["e_"+t]**2 # p

    xp = df.loc[:, df.columns.str.startswith('__xp')].columns.values
    p = df.loc[:, df.columns.str.startswith('__p')].columns.values

    df[colval[0]] = df[xp].sum(axis=1) / df[p].sum(axis=1)
    df[colval[1]] = 1./np.sqrt(df[p].sum(axis=1))

    df = df.drop(xp, axis=1)
    df = df.drop(p, axis=1)
    df = df.replace([np.inf, -np.inf], np.nan)

    return df
