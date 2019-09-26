# The formatting of latex tables for Ilin+2019(2nd) lives here.
import pandas as pd
import astropy.units as u
from opencluster.layout import ValRepr
from num2tex import num2tex
from copy import deepcopy
from opencluster.analysis import *

from .__init__ import logger

def format_upper_lower_err(df, value, up, low):
    '''
    Format three series into one regex expression
    with upper and lower error bounds

    Parameters:
    --------------
    df : DataFrame
        cluster_compilation table
    value : str
        column name with expected value
    up : str
        column name with upper bound
    low : str
        column name with lower bound
    Return:
    -------

    '''

    def string_format_upper_lower_err(series, v, up, low):
        s = deepcopy(series)
        val = str(num2tex(s[v]))
        a = '{' + str(num2tex(s[low])) + '}'
        b = '{' + str(num2tex(s[up])) + '}'
        return r'${}\pm _{}^{}$'.format(val, a, b)

    ndf = deepcopy(df)
    key = deepcopy(value)
    ndf[key] = ndf.apply(lambda x: string_format_upper_lower_err(x, value, up, low),
                       axis=1)
    #drop column that you do not need anymore
    ndf = ndf.drop([up, low],axis=1)
    return ndf

def latexify_results(df, unit):
    '''
    Script that takes the cluster_compilation table
    and makes a beautiful latex table out of it.

    Parameters:
    ------------
    df : DataFrame
        cluster_compilation table
    unit : 's' or 'erg'
        Two different tables will be generated depending
        on this parameter

    Return:
    -------
    formatted DataFrame
    '''
    VR = ValRepr(unit)
    dfselect = deepcopy(df)
    dfselect = df[VR.selection]
    dfselect = dfselect[dfselect.unit == unit]

    for l in ['truncated_2','truncated_Bauke','truncated_MK','Tmax','Tmin','age (Myr)','u_age_high','u_age_low']:
        dfselect.loc[~dfselect[l].isnull(),l] = dfselect.loc[~dfselect[l].isnull(),l].astype(int)
    if unit=='s':
        dfselect.beta2 = dfselect.apply(lambda x: r'${:.2f}\pm{:.2f}$'.format(num2tex(x.beta2, precision=2),x.beta2_err),axis=1)
        dfselect.beta_Bauke = dfselect.apply(lambda x: r'${:.2f}\pm{:.2f}$'.format(num2tex(x.beta_Bauke, precision=2),x.beta_Bauke_err),axis=1)
        dfselect.beta_MK= dfselect.apply(lambda x: r'${:.2f}\pm{:.2f}$'.format(num2tex(x.beta_MK, precision=2),x.beta_MK_err),axis=1)
    elif unit == 'erg':
        #dfselect.FA = dfselect.FA.apply(lambda x: r'${}$'.format(num2tex(x, precision=precision)))
        dfselect.beta2 = dfselect.apply(lambda x: r'${:.2e}$'.format(num2tex(x.beta2)) + r'$\pm{:.2e}$'.format(num2tex(x.beta2_err)),axis=1)
        dfselect.beta_Bauke = dfselect.apply(lambda x: r'${:.2e}$'.format(num2tex(x.beta_Bauke)) + r'$\pm{:.2e}$'.format(num2tex(x.beta_Bauke_err)),axis=1)
        dfselect.beta_MK = dfselect.apply(lambda x: r'${:.2e}$'.format(num2tex(x.beta_MK)) + r'$\pm{:.2e}$'.format(num2tex(x.beta_MK_err)),axis=1)

    dfselect.alpha_Bauke = dfselect.apply(lambda x: r'${:.2f}$'.format(num2tex(x.alpha_Bauke, precision=2))
                                          + r'$\pm{:.2f}$'.format(num2tex(x.alpha_Bauke_err, precision=2)),
                                          axis=1)
    dfselect.alpha_MK= dfselect.apply(lambda x: r'${:.2f}\pm{:.2f}$'.format(num2tex(x.alpha_MK, precision=2),x.alpha_MK_err),axis=1)

    dfselect.FeH = dfselect.apply(lambda x: r'${:.2f}\pm{:.2f}$'.format(x.FeH,x.u_feh),axis=1)
    dfselect.cutoff_MK = dfselect.cutoff_MK.apply(lambda x: u.Quantity(x).value)
    dfselect.cutoff_MK = dfselect.cutoff_MK.apply(lambda x: r'${}$'.format(num2tex(x, precision=2)))

    dfselect = format_upper_lower_err(dfselect, 'age (Myr)', 'u_age_high', 'u_age_low')
    dfselect = dfselect.drop(['beta2_err','u_feh','unit', 'beta_Bauke_err','beta_MK_err','alpha_Bauke_err','alpha_MK_err'],axis=1)
    dfselect = dfselect.fillna('-')
    dfselect = dfselect.replace('$\mathrm{NaN}$','-')
    dfselect = dfselect.replace('$\mathrm{NaN}\pmnan$','-')
    dfselect = dfselect.replace('$\mathrm{NaN}\pm\mathrm{NaN}$','-')
    dfselect = dfselect.replace('$\mathrm{NaN}$$\pm\mathrm{NaN}$','-')

    dfselectrenamed = dfselect.rename(index=str, columns=VR.valrepr)
    return dfselectrenamed

def write_results_to_texfile(df,unit,t=False):
    '''
    Takes a formatted DataFrame and writes it
    to .tex

    Parameters:
    ------------
    df : DataFrame
        cluster_compilation table
    unit : 's' or 'erg'
        Two different tables will be generated depending
        on this parameter
    t : bool
        if True will swap columns and rows
    '''
    with open('tables/latex_allresults_{}.tex'.format(unit), 'w') as f:
        dfselectrenamed = latexify_results(df, unit)
        if t == True:
            dfselectrenamed = dfselectrenamed.set_index('cluster')
            dfselectrenamed = dfselectrenamed.T
        st = dfselectrenamed.to_latex(index=False, escape=False)
        st = st.replace('times','cdot')
        f.write(st)
    return