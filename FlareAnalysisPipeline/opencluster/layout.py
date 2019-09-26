SELECTION = {'s':['h_cluster','Tmin', 'Tmax', 'nstars', 'nflares',
                  'age (Myr)', 'FeH','u_age_high', 'u_age_low', 'u_feh',
                  'cutoff_MK', #'R1s_2', 'R1s_Bauke', 'R1s_MK',#FA,
                  'alpha_Bauke', 'alpha_MK', 'alpha_Bauke_err', 'alpha_MK_err',
                  'beta2', 'beta2_err',
                  'beta_Bauke', 'beta_MK', 'beta_Bauke_err', 'beta_MK_err',
                  'truncated_2', 'truncated_Bauke', 'truncated_MK',
                  'unit',],
             'erg':['h_cluster','Tmin', 'Tmax', 'nstars', 'nflares',
                  'age (Myr)', 'FeH','u_age_high', 'u_age_low', 'u_feh',
                  'cutoff_MK', #'FA', #'R35_2', 'R35_Bauke', 'R35_MK',
                  'alpha_Bauke', 'alpha_MK',  'alpha_Bauke_err', 'alpha_MK_err',
                  'beta2', 'beta2_err',
                  'beta_Bauke', 'beta_MK','beta_Bauke_err', 'beta_MK_err',
                  'truncated_2', 'truncated_Bauke', 'truncated_MK',
                  'unit',]}

UNITMAP = {'s':r'$ED$',
           'erg':r'$E_\mathrm{Kp,flare,}$'}

TITLEMAP = {'s':r'Calculated from $ED$ distribution',
            'erg':r'Calculated from $E_\mathrm{flare}$ distribution'}

class ValRepr(object):
    '''
    Human-friendly variable representations.
    LaTeX enhanced.

    Attributes:
    ------------
    unit : 's' or 'erg'
        define if you are working with flare energies or EDs
    unitmap : dict
        Map units to observables in form of readable strings
    valrepr : dict
        Map column values from results to human-friendly LaTeX
        enhanced strings.
    '''
    def __init__(self, unit='s'):
        self.unit = unit
        self.valrepr = {'R1s_2':r'$R1s_2\;[\mathrm{yr}^{-1}]$',
                       'R1s_Bauke':r'$R1s_\mathrm{B}\;[\mathrm{yr}^{-1}]$',
                       'R1s_MK':r'$R1s_\mathrm{MK}\;[\mathrm{yr}^{-1}]$',
                       'Tmax' : r'$T_\mathrm{max} [\mathrm{K}]$',
                       'Tmin' : r'$T_\mathrm{min} [\mathrm{K}]$',
                       'alpha_Bauke' : r'$\alpha_\mathrm{B}$',
                       'alpha_MK' : r'$\alpha_\mathrm{MK}$',
                       'beta2' : r'$\beta_2\;[\mathrm{yr}^{-1}]$',
                       'beta_Bauke' : r'$\beta_\mathrm{B}\;[\mathrm{yr}^{-1}]$',
                       'beta_MK' : r'$\beta_\mathrm{MK}\;[\mathrm{yr}^{-1}]$',
                       'cutoff_MK' : UNITMAP[unit] + r'$_\mathrm{min}\;[$'+unit+r'$]$',
                       'nstars'  : r'$n_*$',
                       'nflares' : r'$n_\mathrm{flares}$',
                       'truncated_2' : r'$tr_2$',
                       'truncated_Bauke' : r'$tr_\mathrm{B}$',
                       'truncated_MK' : r'$tr_\mathrm{MK}$',
                       'age (Myr)' : r'age [Myr]',
                       'FeH' : r'[Fe/H]',
                       'h_cluster' : 'cluster',
                       'FA' : r'$FA$',
                       'FR' : r'$FR\;[\mathrm{yr}^{-1}]$',
                       'tot_obs_time' : r'$t_\mathrm{obs}\;[\mathrm{yr}]$',
                       'flare_frac' : r'$n_\mathrm{*,flaring} / n_*$',
                       'dist (pc)' : r'distance [pc]',
                       'flarelum' : r'$L_\mathrm{Kp,flare}\;[\mathrm{erg\;s}^{-1}]$',
                       'm67_freq_s':r'cumulative flare frequency$\;[\mathrm{yr}^{-1}]$',
                       'm67_freq_erg':r'cumulative flare frequency$\;[\mathrm{yr}^{-1}]$'}
        self.selection = SELECTION[unit]
        self.title = TITLEMAP[unit]
        self.FFDxlabel = UNITMAP[unit] + r'$\;[$'+unit+r'$]$'