# The plot prodcution and formatting for Ilin+2019(2nd) lives here
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from copy import deepcopy

from .opencluster import OpenCluster
from .ffd import FFD
from .utils import timestamp
from .layout import ValRepr
MARKERS_LIST = ['p','<','>', 'X','*','o','s','d']
COLORS_LIST=['k','darkred','red','orangered','orange','yellowgreen','green','blue']
TEFFS_LIST = [2500,3000,3250,3500,3750,4000,5000,6000]
MARKERS_DICT = dict(zip(TEFFS_LIST,MARKERS_LIST))
COLORS_DICT = dict(zip(TEFFS_LIST,COLORS_LIST))
CB_COLOR_CYCLE = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

import matplotlib as mpl

def import_matplotlib_style(label_size = 16):


    font = {'family' : 'monospace',
            'weight' : 'normal',
            'size'   : label_size}

    mpl.rc('font', **font)

    params = {'legend.fontsize': label_size,
              'figure.figsize': (7, 5),
              'axes.labelsize': label_size,
              'axes.titlesize':label_size,
              'xtick.labelsize':label_size,
              'ytick.labelsize':label_size}
    mpl.rcParams.update(params)

def plot_age_mass_AM(val, valerr, data, unit, ax=None,
                     fmt=MARKERS_DICT, s=10, fontsize=13,
                     starcolors=COLORS_DICT):
    '''
    Plot activity measure (AM) from data against
    age and effective temperature. Save file to folder.

    Parameters:
    ------------
    val : str
        column name in data that contains AM
    valerr : str
        column name in data that contains AM error
    data : DataFrame
        full table containing at least the
        columns 'age (Myr)', 'Tmin', 'Tmax',
        and val
    unit : 's' or 'erg'
        either work on ED or flare energy
    ax : None or Axis
        Panel to plot to
    fmt : list or str
        marker format as in matplotlib
    s : float
        marker size
    fontsize : int

    starcolors : list or str
        plot color(s)

    Return:
    --------
    modified ax
    '''
    VR = ValRepr(unit) #create the right layout settings
    data = data[~data[val].isnull()]
    data = data[data.unit == unit]
    data = data[data.Tmax-data.Tmin <1000] #avoid helper bins

    n = data.groupby(['Tmin','Tmax']).count().shape[0]
    for it in [fmt, starcolors]:
        if isinstance(it,str):
            it = [it]*n
        elif isinstance(it,list):
            it = it[:n]
    for label, group in data.groupby(['Tmin','Tmax']):
        group = group.sort_values(by='age (Myr)')
        ax.errorbar(x=group['age (Myr)'], y=group[val],
                    xerr=(group['u_age_low'],group['u_age_high']),
                    yerr=group[valerr],
                    color=starcolors[int(group.Tmin.iloc[0])],fmt=fmt[int(group.Tmin.iloc[0])],markersize=s,
                    label='{}-{} K'.format(int(group.Tmin.iloc[0]),int(group.Tmax.iloc[0]+1)))
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(VR.valrepr[val], fontsize=fontsize)
    ax.set_xlabel('age (Myr)', fontsize=fontsize)
    ax.legend(loc=1,fontsize=fontsize)
    return ax

def save_PP_plots(ffds):
    '''
    Create and save PP plots for all fitted
    power laws.
    '''
    #validate the ffds:
    map_alphamode = {'fix2' : r'$\alpha=2$', 'Bauke' : 'Bauke (2007)', 'MK':'Maschberger and Kroupa (2009)'}
    assert len(list(set([ffd.name for ffd in ffds]))) == 1
    assert len(list(set([ffd.Tmin for ffd in ffds]))) == 1
    assert len(list(set([ffd.mode for ffd in ffds]))) == 1
    fig, axes = plt.subplots(nrows=1, ncols=len(ffds), figsize=(len(ffds)*4,3.5), sharey=True, sharex=True)
    for ffd, ax in zip(ffds, axes):
        ffd.plot_percentile_percentile(ax,sig_level=0.05,c='r')
        ax.set_title(r'{}, $tr={}$'.format( map_alphamode[ffd.alpha_mode],int(ffd.truncated)))
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/PP/{}_{}_{}_{}.png'.format(ffd.name, ffd.Tmin, ffd.mode, timestamp()))
    return

def plot_all_FFDs(OC, ax, mode='ED'):
    '''
    Plot all FFDs in one plot
    '''
    l = list(OC.ffds.items())
    l.sort()
    for key, ffd in l:
        ff = deepcopy(ffd)
        if ((mode=='energy') and (ffd.flares.Lum_Kepler.values.shape[0]==0)):
            #no energies will be calculated
            continue
        ff.mode = mode
        ff.plot_FFD(ax, setplotlims=False, c=COLORS_DICT[ff.Tmin], marker=MARKERS_DICT[ff.Tmin], s=100,
                    label='{:.0f}-{:.0f}'.format(key,ff.Tmax))
    plt.legend()
    plt.title(OC.h_cluster, fontsize=13)
    return

def save_all_FFD_plots_in_ED_and_energy(OC):
    fig, ax = plt.subplots(figsize=(7,6))
    plot_all_FFDs(OC, ax, mode='ED')
    plt.savefig('plots/FFDs/all_ED_FFDs_{}_{}.png'.format(OC.cluster, timestamp()), dpi=300)
    fig, ax = plt.subplots(figsize=(7,6))
    plot_all_FFDs(OC, ax, mode='energy')
    plt.savefig('plots/FFDs/all_energy_FFDs_{}_{}.png'.format(OC.cluster, timestamp()), dpi=300)
    return

def save_PP_plot(ffd):
    fig, ax = plt.subplots(figsize=(7,6))
    ffd.plot_percentile_percentile(ax)
    plt.title('PP plot: {}\n {} cutoff at {}\n Fitting method: {}'.format(ffd.h_name,ffd.mode,ffd.cutoff,ffd.alpha_mode))
    plt.savefig('plots/PP/PPplot_{}_{}_{}_{}.png'.format(ffd.name, ffd.mode, ffd.alpha_mode.strip(), timestamp()), dpi=300)
    return

def plot_interconnected_subplots(df, xs, ys, err=False,
                                 xerr=None, yerr=None,
                                 color="k", alpha=.4, width=3, height=4):
    """Plot multiple subplots in a row.

    Parameters:
    ----------
    df : DataFrame
        has columns xs and ys
    xs : list of strings
        column names for abscissa
    ys : list of string
        column names for ordinate
    color : str
        lines and scatter dot color
    alpha : float
        alpha of lines
    width : float
        width of indiv. subplot
    height : float
        height of each subplot
    """
    n = len(ys)
    assert n == len(xs)

    fig, axes = plt.subplots(nrows=1, ncols=n, sharey=True,
                             figsize=(n*width, height))
    for l, row in df.iterrows():

        ddt = pd.DataFrame({"x":row[xs].values,"y":row[ys].values})
        for m in range(n-1,0,-1):
            j = 0
            vals = np.array([ddt.loc[m, "x"], ddt.loc[m, "y"],
                             ddt.loc[m-1, "x"], ddt.loc[m-1, "y"]])

            while (np.isnan(vals) == True).any() & (j < m-1):
                j += 1
                vals = np.array([ddt.loc[m, "x"], ddt.loc[m, "y"],
                                 ddt.loc[m-1-j, "x"], ddt.loc[m-1-j, "y"]])

            if (np.isnan(vals) == False).all():
                con = ConnectionPatch(xyA=vals[:2], xyB=vals[2:],
                                      coordsA="data", coordsB="data",
                                      axesA=axes[m], axesB=axes[m-1-j],
                                      color=color, alpha=alpha, zorder=0)
                axes[m].add_artist(con)

    if err==False:

        for i,col,ind,ax in zip(range(1, n+1), ys, xs, axes):
            ax.scatter(df[ind], df[col], c=color, zorder=20)

    elif err==True:

        assert n == len(xerr)
        assert n == len(yerr)
        for i,col, ind, yr, xr, ax in zip(range(1, n+1), ys, xs, yerr, xerr,axes):
            ax.errorbar(x=df[ind], y=df[col], xerr=df[xr], yerr=df[yr],
                       c=color, fmt="o", zorder=20)

    return fig, axes