import os
from   pathlib import Path

import numpy                as np
import matplotlib           as mpl
import matplotlib.pyplot    as plt
import seaborn              as sns
import pandas               as pd

from scipy.stats    import kruskal, wilcoxon, mannwhitneyu, linregress
from scipy.optimize import curve_fit
from lmfit import Model, Parameters
from statsmodels.stats.multitest import multipletests


sns.set_context('paper')
mpl.rcParams['font.size'] = 16
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['lines.linewidth'] = 2

# make a colour map viridis
viridis = mpl.colormaps["viridis"]
flare   = mpl.colormaps["rocket"]
crest   = mpl.colormaps["mako"]
magma   = mpl.colormaps["magma"]
rocket_r    = mpl.colormaps['rocket_r']

color_freq = {1:magma(0.05), 5:magma(0.1), 10:magma(0.2), 20:magma(.4), 30:magma(.5), 40:magma(.6), 50:magma(.7), 100:magma(.9)}
color_squares = {1:viridis(0.2), 5:viridis(.4), 7:viridis(.6), 15:viridis(.8), 20:viridis(1.0)}
Fs = 2e4

freq_sweep_pulses = np.arange(9)

def sdnfunc(expected, gamma):
    return gamma * expected / (gamma + expected)

def nosdn(expected, m):
    return m * expected


# ── Inlined helpers (from utils, pattern_index, plot_tools) ───────────────────

def expand_list_column(df, col, prefix):
    expanded = pd.DataFrame(df[col].tolist(), index=df.index)
    expanded.columns = [f'{prefix}{i}' for i in range(expanded.shape[1])]
    return pd.concat([df, expanded], axis=1)


# Mapping from pattern ID → list of constituent 1-sq spot IDs.
_patternID = {
    1:[101], 2:[105], 3:[109], 4:[113], 5:[117],
    6:[147], 7:[151], 8:[155], 9:[159], 10:[163],
    11:[197], 12:[201], 13:[205], 14:[209], 15:[213],
    16:[243], 17:[247], 18:[251], 19:[255], 20:[259],
    21:[293], 22:[297], 23:[301], 24:[305], 25:[309],
    26:[339], 27:[343], 28:[347], 29:[351], 30:[355],
    31:[389], 32:[393], 33:[397], 34:[401], 35:[405],
    36:[435], 37:[439], 38:[443], 39:[447], 40:[451],
    41:[485], 42:[489], 43:[493], 44:[497], 45:[501],
    46:[209,247,259,301,393], 47:[205,251,297,389,447],
    48:[197,255,347,401,439], 49:[201,293,351,355,443],
    50:[251,305,343,397,451],
    51:[105,109,113,117,155,159,243,309,343,351,355,405,443,451,485],
    52:[101,109,117,147,155,197,305,309,339,343,351,401,451,485,497],
    53:[151,163,197,201,209,213,259,301,339,347,393,401,435,439,489],
    54:[113,159,205,209,243,251,255,301,347,355,393,405,439,443,447],
    55:[105,151,163,201,213,247,259,293,297,389,397,435,489,493,501],
    56:[101,113,147,159,209,243,485], 57:[101,113,159,205,209,213,497],
    58:[147,163,209,243,255,443,485], 59:[109,201,247,301,309,355,501],
    60:[117,151,259,347,351,389,439], 61:[163,197,201,247,301,447,501],
    62:[117,151,259,355,393,405,439], 63:[105,117,339,347,351,389,493],
    64:[147,163,197,255,443,447,501], 65:[109,201,309,355,393,405,439],
    66:[101,205,213,397,401,451,497], 67:[155,251,293,297,305,489,493],
    68:[105,293,297,339,389,489,493], 69:[155,251,305,343,435,451,489],
    70:[305,343,397,401,435,451,497],
    71:[101,113,147,159,163,197,205,209,213,243,255,443,447,485,501],
    72:[101,113,147,159,205,209,213,243,255,343,397,401,451,485,497],
    73:[109,147,163,197,201,209,243,247,255,301,309,443,447,485,501],
    74:[109,117,151,201,247,259,301,309,347,351,355,389,393,405,439],
    75:[101,113,155,159,205,213,251,305,343,397,401,435,451,489,497],
    76:[105,117,155,251,259,293,297,305,339,347,351,389,435,489,493],
    77:[109,151,163,197,201,247,301,309,355,393,405,439,443,447,501],
    78:[105,117,151,259,293,297,339,347,351,355,389,393,405,439,493],
    79:[105,155,251,293,297,305,339,343,397,401,435,451,489,493,497],
    80:[101,147,197,243,293,339,435], 81:[105,151,201,247,297,343,439],
    82:[109,155,205,251,301,347,443], 83:[113,159,209,255,305,351,447],
    84:[117,163,213,259,309,355,451], 85:[101,151,197,293,389,435,485],
    86:[105,155,201,297,393,439,489], 87:[109,159,205,301,397,443,493],
    88:[113,163,209,305,401,447,497], 89:[117,147,213,309,405,451,501],
    90:[151,247,293,343,389,439,485], 91:[155,251,297,347,393,443,489],
    92:[147,243,309,339,405,435,501], 93:[159,255,301,351,397,447,493],
    94:[163,259,305,355,401,451,497],
    95:[101,105,147,151,197,201,243,247,293,339,343,389,435,439,485],
    96:[105,109,155,201,205,247,251,297,301,343,347,393,439,443,489],
    97:[101,117,147,151,197,213,243,293,309,339,389,405,435,485,501],
    98:[105,151,155,197,201,247,251,293,297,343,389,393,439,485,489],
    99:[109,113,159,205,209,255,301,305,347,351,397,401,443,447,493],
    100:[101,117,147,163,213,243,259,309,339,355,405,435,451,497,501],
    101:[109,155,159,205,251,255,297,301,347,351,393,397,443,489,493],
    102:[113,117,163,209,213,259,305,309,355,401,405,447,451,497,501],
    103:[113,159,163,209,255,259,305,351,355,397,401,447,451,493,497],
    104:[105,163,255,347,401], 105:[109,201,259,351,447],
    106:[113,205,251,297,355], 107:[155,209,301,393,443],
    108:[159,251,305,397,451], 109:[203,257,353,395,438],
    110:[200,230,298,306,342], 111:[150,201,253,271,307],
    112:[105,163,255,347,401],
    999:[101,105,109,113,117,147,151,155,159,163,197,201,205,209,213,
         243,247,251,255,259,293,297,301,305,309,339,343,347,351,355,
         389,393,397,401,405,435,439,443,447,451,485,489,493,497,501],
}


def _get_patternID(sq_set):
    for k, v in _patternID.items():
        if v == sq_set:
            return int(k)


def get_patternIDlist_for_nSq_pattern(patternIDnSq):
    return [_get_patternID([spot]) for spot in _patternID[patternIDnSq]]


def ax_to_partial_dist_heatmap_ax(
        pivotdf, numdf, fig, ax,
        barw=0.03, pad=0.01, shrink=0.8, palette='viridis',
        annotate=False, show_marginals=True, cbar_label=''):
    bbox = ax.get_position()
    x0, y0 = bbox.x0, bbox.y0
    w,  h  = bbox.width, bbox.height
    ax.remove()

    if show_marginals:
        ax  = fig.add_axes([x0, y0, shrink*w, shrink*h])
        axx = fig.add_axes([x0, y0+shrink*h+pad, shrink*w, barw], aspect='auto')
        axy = fig.add_axes([x0+shrink*w+pad, y0, barw, shrink*h], aspect='auto')
        axc = fig.add_axes([x0+shrink*w+barw+2*pad, y0, barw, shrink*h], aspect='auto')
    else:
        main_w = w - barw - 2*pad
        ax  = fig.add_axes([x0, y0, main_w, h])
        axx = axy = None
        axc = fig.add_axes([x0+main_w+pad, y0, barw, h], aspect='auto')

    maxlim = np.round(np.max(pivotdf.values), 2)
    minlim = np.round(np.min(pivotdf.values), 2)
    ax.imshow(pivotdf, cmap=palette, vmin=minlim, vmax=maxlim, aspect='auto')

    if show_marginals:
        pw = pivotdf.mean(axis=0).values.reshape(1, -1)
        fw = pivotdf.mean(axis=1).values.reshape(-1, 1)
        axx.imshow(pw, cmap=palette, vmin=minlim, vmax=maxlim)
        axy.imshow(fw, cmap=palette, vmin=minlim, vmax=maxlim, origin='lower')

    if annotate:
        pw_n = numdf.sum(axis=0).values.reshape(1, -1)
        fw_n = numdf.sum(axis=1).values.reshape(-1, 1)
        for i in range(fw_n.shape[0]):
            for j in range(pw_n.shape[1]):
                ax.text(j, i, f'{pivotdf.values[i,j]:.2f}',
                        ha='center', va='center', color='white', fontsize=12)
                ax.text(j-0.2, i-0.2, f'{numdf.values[i,j]:.0f}',
                        ha='center', va='center', color='yellow', fontsize=10)
                if show_marginals:
                    if i == 0:
                        axx.text(j, 0, f'{pw[0,j]:.2f}',
                                 ha='center', va='center', color='white', fontsize=12)
                    axy.text(0, i, f'{fw[i,0]:.2f}',
                             ha='center', va='center', color='white', fontsize=12)

    ax.set_xticks(np.arange(9), labels=np.arange(9))
    ax.set_ylim([-0.5, 3.5])
    ax.set_yticks([0, 1, 2, 3], labels=[20, 30, 40, 50])
    ax.set_xlabel('Pulse Index')
    ax.set_ylabel('Frequency (Hz)')
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(ax.get_images()[0], cax=axc)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=14)

    if show_marginals:
        for _ax in [axx, axy, axc]:
            for spine in _ax.spines.values():
                spine.set_visible(False)
        axx.get_xaxis().set_visible(False)
        axx.get_yaxis().set_visible(False)
        axy.get_xaxis().set_visible(False)
        axy.get_yaxis().set_visible(False)
        axx.set_aspect('auto')
        axy.set_aspect('auto')
    else:
        for spine in axc.spines.values():
            spine.set_visible(False)

    return ax, axx, axy, axc, cbar


# load datapaths

# Paths
paper_figure_export_location = Path("./figure3_output")
paper_figure_export_location.mkdir(parents=True, exist_ok=True)
data_path_FS = Path("../../DATA_Nov2025")

### Load CC FreqSweep data
CC_FS_shortdf_withkernelfit_datapath = data_path_FS / "all_cells_FreqSweep_CC_kernelfit_response_measurements.h5"
cc_FS_shortdf = pd.read_hdf(CC_FS_shortdf_withkernelfit_datapath, key='data')
print(cc_FS_shortdf.shape)

cc_FS_shortdf_slice = cc_FS_shortdf[
    (cc_FS_shortdf['location'] == 'CA1') &
    (cc_FS_shortdf['numSq'].isin([1,5,7,15])) &
    (cc_FS_shortdf['stimFreq'].isin([20,30,40,50])) &
    (cc_FS_shortdf['condition'] == 'Control') &
    (cc_FS_shortdf['ch0_response']==1) &
    (cc_FS_shortdf['IR'] >50) & (cc_FS_shortdf['IR'] < 400) &
    (cc_FS_shortdf['tau'] < 40) &
    (cc_FS_shortdf['spike_in_baseline_period'] == 0) &
    (cc_FS_shortdf['ac_noise_power_in_ch0'] < 40)
]
print(cc_FS_shortdf.shape, '--screened-->', cc_FS_shortdf_slice.shape)
cc_FS_shortdf_slice = cc_FS_shortdf_slice.copy()
cc_FS_shortdf_slice['patternList'] = cc_FS_shortdf_slice['patternList'].astype('int32')
patternIDs = np.sort(cc_FS_shortdf_slice[cc_FS_shortdf_slice['numSq'] != 1]['patternList'].unique())
cc_FS_shortdf_slice = expand_list_column(cc_FS_shortdf_slice, 'peaks_field_norm', 'pfn_')
print(f"CC FS cells: {cc_FS_shortdf_slice['cellID'].nunique()}  sweeps: {cc_FS_shortdf_slice['trialID'].nunique()}")

### Load LTM CC long data (provides 7-sq and 15-sq patterns)
cc_LTM_longdf = pd.read_hdf(Path("../../2022/VC_DATA/all_cells_LTMRand_CC_long.h5"), key='data')
cc_LTM_slice = cc_LTM_longdf[
    (cc_LTM_longdf['location'] == 'CA1') &
    (cc_LTM_longdf['numSq'].isin([7, 15])) &
    (cc_LTM_longdf['stimFreq'].isin([20, 30, 40, 50])) &
    (cc_LTM_longdf['condition'] == 'Control') &
    (cc_LTM_longdf['IR'] > 50) & (cc_LTM_longdf['IR'] < 400) &
    (cc_LTM_longdf['tau'] < 40)
].copy()
cc_LTM_slice = expand_list_column(cc_LTM_slice, 'peaks_cell', 'PSC_')
cc_LTM_slice = expand_list_column(cc_LTM_slice, 'peaks_cell_norm', 'normPSC_')
cc_LTM_slice['PSC_8'] = np.nan  # LTM has 8 pulses; pad to match FS column set
cc_LTM_slice['numChannels'] = 2  # no field data in LTM
cc_LTM_slice['patternList'] = cc_LTM_slice['patternList'].astype('int32')
print(f"CC LTM cells: {cc_LTM_slice['cellID'].nunique()}  sweeps: {cc_LTM_slice.shape[0]}")

### Load VC FreqSweep data
VC_FS_shortdf_withkernelfit_datapath = data_path_FS / "all_cells_FreqSweep_VC_kernelfit_response_measurements.h5"
vc_FS_shortdf = pd.read_hdf(VC_FS_shortdf_withkernelfit_datapath, key='data')
print(vc_FS_shortdf.shape)

vc_FS_shortdf_slice = vc_FS_shortdf[
    (vc_FS_shortdf['location'] == 'CA1') &
    (vc_FS_shortdf['numSq'].isin([1,5,15])) &
    (vc_FS_shortdf['stimFreq'].isin([20,30,40,50])) &
    (vc_FS_shortdf['condition'] == 'Control') &
    (vc_FS_shortdf['ch0_response']==1) &
    (vc_FS_shortdf['IR'] >40) & (vc_FS_shortdf['IR'] < 400) &
    (vc_FS_shortdf['tau'] < 40) &
    (vc_FS_shortdf['ac_noise_power_in_ch0'] < 40) &
    (vc_FS_shortdf['valley_0'].notnull())
]
print(vc_FS_shortdf.shape, '--screened-->', vc_FS_shortdf_slice.shape)
print(f"VC cells: {vc_FS_shortdf_slice['cellID'].nunique()}  sweeps: {vc_FS_shortdf_slice['trialID'].nunique()}")

def calculate_expected_response(celldf, pulse_index, freq, patternID,):
    """
    Calculate the expected response of a pattern based on the response to individual spots in the pattern
    """
    # constants
    Fs      = 2e4
    cellID  = celldf['cellID'].iloc[0]
    
    # checks
    field_data=True if celldf['numChannels'].iloc[0] == 4 else False

    # check if the given cell has 1sq data
    if not 1 in celldf['numSq'].unique():
        # print('No 1Sq data for this cell', celldf['numSq'].unique())
        # generate dataerror to be caught by the calling function
        raise ValueError(f'Cell: {cellID} - No 1Sq data for this cell. {pulse_index}, {freq}, {patternID}')
    # data
    pattern_response_df             = celldf[(celldf['patternList'] == patternID) & (celldf['stimFreq'] == freq)  ]
    if pattern_response_df.shape[0] == 0:
        raise ValueError(f'Cell: {cellID} - No data for this pattern {patternID} and freq {freq} Hz')
    
    # get the pattern
    constituent_spots_of_pattern    = get_patternIDlist_for_nSq_pattern(patternID) #1sq spots that make the pattern in the patternID
    numSq                           = len(constituent_spots_of_pattern)

    obs_col = 'PSC_' + str(pulse_index)
    obs_col_field = 'pfn_' + str(pulse_index)

    # select only the columns needed for expected response calculation
    _keep = (['cellID', 'numSq', 'stimFreq', 'patternList'] +
             [f'PSC_{i}' for i in range(9)] +
             [f'pfn_{i}' for i in range(9)])
    celldf = celldf[[c for c in _keep if c in celldf.columns]].copy()
    celldf.loc[:, 'patternList']    = celldf['patternList'].astype('int32')
    
    # step 0: get the observed response from the pattern_response_df
    observed_response_cell      = pattern_response_df.loc[:, obs_col].values
    if field_data:
        observed_response_field     = pattern_response_df.loc[:, obs_col_field].values
        observed_response_scaled    = observed_response_cell / observed_response_field
    else:
        observed_response_scaled    = observed_response_cell * np.nan
    
    # expected response calculation
    # step 1: slice the dataframe to get only those rows where 'patternList' is in the list 'constituent_spots_of_pattern'
    df1sq = celldf.loc[celldf['patternList'].isin(constituent_spots_of_pattern), :].copy()
    
    # step 2: get the peaks for each row between columns probePulseStart and probePulseStart+ipi
    # here i am taking the mean of all the trials of the constituent patterns and then summing those means
    expected_response = df1sq.loc[:,('patternList','PSC_0')].groupby(by='patternList').mean().sum()['PSC_0']
    
    return numSq, freq, patternID, pulse_index, field_data, observed_response_cell, observed_response_scaled, expected_response

def sdn_fits(xdata, ydata, f,s,t):
    # # if xdata and ydata lenghts are not same, 
    if len(xdata) != len(ydata):
        raise ValueError('Length of xdata and ydata are not same')
        return np.nan, np.nan, np.nan, np.nan
        
    # if xdata or yadata is empty or have length 0, 
    if len(xdata) <3 or len(ydata) <3:
        return np.nan, np.nan, np.nan, np.nan
        
    # Create an lmfit model for the sdnfunc
    model = Model(sdnfunc)
    # Create a set of parameters
    params = Parameters()
    params.add('gamma', value=5)

    # Create an lmfit model for the data
    model_linear = Model(nosdn)
    # Create a set of parameters
    params_linear = Parameters()
    params_linear.add('m', value=1)

    # Fit the sdnfunc to  data using lmfit and method = cobyla
    result = model.fit(ydata, params, expected=xdata, method='cobyla')

    # also try fitting xdata and ydata to a linear model
    result_linear = model_linear.fit(ydata, params_linear, expected=xdata, method='cobyla')

    # Extract the fitted parameters
    fitted_gamma = np.round( result.best_values['gamma'], 3)
    fitted_slope = np.round( result_linear.best_values['m'], 3)

    return fitted_gamma, fitted_slope, result.rsquared, result_linear.rsquared


def gamma_distribution(df_sdn, fitdf=None, x='expected_response', y='observed_response', first='cellID', second='pulse_index', third='freq'):
    gamma_dist = []

    # get gamma distribution for the entire dataset
    dfslice = df_sdn.dropna(subset=[x,y])
    g,m,r2g,r2lin = sdn_fits(dfslice[x], dfslice[y], 'all','all','all')
    gamma_dist.append({'expected': x, 'observed':y, first:1000, second:1000, third:1000, 'sample_size':dfslice.shape[0], 'gamma':g, 'slope':m, 'r2_sdn':r2g, 'r2_lin':r2lin})
    f,s,t = np.nan, np.nan, np.nan

    for f in np.sort(df_sdn[first].unique()):
        dfslice = df_sdn[(df_sdn[first] == f)]
        # remove nan and inf from data
        dfslice = dfslice[(np.abs(dfslice[x]) != np.inf) & (np.abs(dfslice[y]) != np.inf)].dropna(subset=[x,y])
        # if most of x and y data is 0, the model will not converge, so we need to remove all zero entries
        dfslice = dfslice[(dfslice[x] != 0) & (dfslice[y] != 0)]
        
        g,m,r2g,r2lin = sdn_fits(dfslice[x], dfslice[y], f,'all','all')
        gamma_dist.append({'expected': x, 'observed':y, first:f, second:1000, third:1000, 'sample_size':dfslice.shape[0], 'gamma':g, 'slope':m, 'r2_sdn':r2g, 'r2_lin':r2lin})
        
        for s in np.sort(df_sdn[second].unique()):
            dfslice = df_sdn[(df_sdn[first] == f) & (df_sdn[second] == s)].dropna(subset=[x,y])
            # remove np.inf from data
            dfslice = dfslice[(np.abs(dfslice[x]) != np.inf) & (np.abs(dfslice[y]) != np.inf)].dropna(subset=[x,y])
            # if most of x and y data is 0, the model will not converge, so we need to remove all zero entries
            dfslice = dfslice[(dfslice[x] != 0) & (dfslice[y] != 0)]
            
            g,m,r2g,r2lin = sdn_fits(dfslice[x], dfslice[y], f,s,'all')
            gamma_dist.append({'expected': x, 'observed':y, first:f, second:s, third:1000, 'sample_size':dfslice.shape[0], 'gamma':g, 'slope':m, 'r2_sdn':r2g, 'r2_lin':r2lin})
            
            for t in np.sort(df_sdn[third].unique()):
                dfslice = df_sdn[(df_sdn[first] == f) & (df_sdn[second] == s) & (df_sdn[third] == t)].dropna(subset=[x,y])
                # remove nan
                # remove np.inf from data
                dfslice = dfslice[(np.abs(dfslice[x]) != np.inf) & (np.abs(dfslice[y]) != np.inf)].dropna(subset=[x,y])
                # if most of x and y data is 0, the model will not converge, so we need to remove all zero entries
                dfslice = dfslice[(dfslice[x] != 0) & (dfslice[y] != 0)]
                g,m,r2g,r2lin = sdn_fits(dfslice[x], dfslice[y], f,s,t)
                gamma_dist.append({'expected': x, 'observed':y, first:f, second:s, third:t, 'sample_size':dfslice.shape[0], 'gamma':g, 'slope':m, 'r2_sdn':r2g, 'r2_lin':r2lin})


    # create a dataframe from the list of dicts
    df_gamma_dist = pd.DataFrame(gamma_dist)

    if fitdf is not None:
        fitdf2 = pd.concat([fitdf, df_gamma_dist], axis=0)
    else:
        fitdf2 = df_gamma_dist

    return fitdf2


def generate_ebyi_df(vc_shortdf):
    idvars = ['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','trialID']
    valvars = [f'PSC_{i}' for i in freq_sweep_pulses]
    df2 = vc_shortdf.melt(id_vars=idvars, 
                value_vars=valvars, 
                var_name='pulse', value_name='PSC')

    # if clampPotential=-70, remove rows with positive PSC values
    df2 = df2[((df2['clampPotential']==-70) & (df2['PSC']<0)) | ((df2['clampPotential']==0) & (df2['PSC']>0))]

    df2['pulse'] = df2['pulse'].apply(lambda x: int(x.split('_')[-1]))
    df2 = df2.dropna(subset=['PSC'])
    df2.drop(columns=['trialID'], inplace=True)
    df2['numSq'] = df2['numSq'].astype('int')
    df4 = df2.groupby(['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','pulse']).mean().reset_index()
    ebyi_df = df4.pivot(index=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','pulse'], columns='clampPotential', values='PSC').reset_index()
    ebyi_df = ebyi_df.dropna(subset=[-70,0])
    # ratio of -70 and 0
    ebyi_df['EbyI'] = ( - ebyi_df[-70] / ebyi_df[0])
    # ebyi_df = ebyi_df[(ebyi_df['EbyI'] < 20) & (ebyi_df['EbyI'] > 0)]

    print(ebyi_df.shape)
    return ebyi_df


def generate_cc_delay_df(cc_short_df):
    idvars = ['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','trialID']

    valvars2a = [f'peakdelay_{i}' for i in freq_sweep_pulses]
    valvars2b = [f'onsetdelay_{i}' for i in freq_sweep_pulses]
    valvars2c = [f'normPSC_{i}' for i in freq_sweep_pulses]

    df2a = cc_short_df.melt(id_vars=idvars, 
                value_vars=valvars2a, 
                var_name='peak', value_name='peak_delay')

    df2b = cc_short_df.melt(id_vars=idvars, 
                value_vars=valvars2b, 
                var_name='onset', value_name='onset_delay')
    
    df2c = cc_short_df.melt(id_vars=idvars, 
                value_vars=valvars2c, 
                var_name='peakres', value_name='peak_PSP')

    df2a['pulse'] = df2a['peak' ].apply(lambda x: int(x.split('_')[-1]))
    df2b['pulse'] = df2b['onset'].apply(lambda x: int(x.split('_')[-1]))
    df2c['pulse'] = df2c['peakres'].apply(lambda x: int(x.split('_')[-1]))

    # # concat df1 and df2 on axis1
    cc_delay_df  = pd.merge(df2a,  df2b, on=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','trialID','pulse'], )
    cc_delay_df  = pd.merge(cc_delay_df,  df2c, on=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','trialID','pulse'], )
    # print('after merger: ', cc_delay_df.shape)

    # # remove those trials for which peak_delay and onset_delay are negative
    cc_delay_df.drop(columns=['peak','onset','peakres'], inplace=True)
    cc_delay_df = cc_delay_df[(cc_delay_df['peak_delay']>0) &(cc_delay_df['peak_delay']<0.05) ]
    cc_delay_df = cc_delay_df[(cc_delay_df['onset_delay']>0)&(cc_delay_df['onset_delay']<0.05) ]
    cc_delay_df = cc_delay_df[cc_delay_df['peak_delay'] > cc_delay_df['onset_delay']]
    cc_delay_df = cc_delay_df[(cc_delay_df['peak_PSP']>0)&(cc_delay_df['peak_PSP']<20) ]

    # print('after removing negative onset to peak times: ', cc_delay_df.shape)

    # cc_delay_df.drop(columns=['trialID'], inplace=True)
    # # drop NaNs in time_to_peak
    cc_delay_df = cc_delay_df.dropna(subset=['peak_delay','onset_delay', 'peak_PSP'])
    # print('merged and peak and onset calculated df: ', cc_delay_df.shape)
    columnscc = ['peak_delay','onset_delay']
    cc_delay_df[columnscc] = cc_delay_df[columnscc] * 1000

    print(cc_delay_df.shape)
    return cc_delay_df


def generate_cc_psp_ltm(cc_ltm_slice):
    """Extract peak_PSP (normalised amplitude) from LTM CC long data.
    No delay columns available — those are filled with NaN so downstream
    delay-based filtering naturally excludes these rows."""
    ltm_pulses = range(8)
    idvars = ['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','sweep']
    valvars = [f'normPSC_{i}' for i in ltm_pulses]
    df = cc_ltm_slice.melt(id_vars=idvars, value_vars=valvars, var_name='peakres', value_name='peak_PSP')
    df['pulse'] = df['peakres'].apply(lambda x: int(x.split('_')[-1]))
    df.drop(columns=['peakres'], inplace=True)
    df = df[(df['peak_PSP'] > 0) & (df['peak_PSP'] < 20)]
    df.rename(columns={'sweep': 'trialID'}, inplace=True)
    df['peak_delay']  = np.nan
    df['onset_delay'] = np.nan
    print(f'LTM CC PSP rows: {df.shape}')
    return df


def generate_vc_delay_df(vc_shortdf):    
    idvars = ['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','trialID']

    valvars2a = [f'peakdelay_{i}' for i in freq_sweep_pulses]
    valvars2b = [f'onsetdelay_{i}' for i in freq_sweep_pulses]

    df2a = vc_shortdf.melt(id_vars=idvars, 
                value_vars=valvars2a, 
                var_name='peak', value_name='peak_delay')

    df2b = vc_shortdf.melt(id_vars=idvars, 
                value_vars=valvars2b, 
                var_name='onset', value_name='onset_delay')

    df2a['pulse'] = df2a['peak' ].apply(lambda x: int(x.split('_')[-1]))
    df2b['pulse'] = df2b['onset'].apply(lambda x: int(x.split('_')[-1]))

    # concat df1 and df2 on axis1
    df3  = pd.merge(df2a,  df2b, on=['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','trialID','pulse'], )

    # remove those trials from df3 where time to peak is smaller than time to valley
    # remove those trials for which peak_onset and valley_onset are negative
    df3.drop(columns=['trialID'], inplace=True)
    df3.drop(columns=['onset','peak'], inplace=True)
    df3 = df3[(df3['peak_delay']>0) &(df3['peak_delay']<0.05)]
    df3 = df3[(df3['onset_delay']>0)&(df3['onset_delay']<0.05)]
    df3 = df3[ df3['peak_delay']    > df3['onset_delay']]
    # # drop NaNs in time_to_peak
    df3 = df3.dropna(subset=['peak_delay','onset_delay'])

    df4 = df3.groupby(['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','pulse']).median().reset_index()
    # # pivot w.r.t clampPotential
    peak_delay_df  = df4.pivot(index=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','pulse'], columns='clampPotential', values='peak_delay').reset_index()
    onset_delay_df = df4.pivot(index=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','pulse'], columns='clampPotential', values='onset_delay').reset_index()

    # drop NaNs from df4pivot from columns -70 and 0
    peak_delay_df  = peak_delay_df.dropna(subset=[-70,0])
    onset_delay_df = onset_delay_df.dropna(subset=[-70,0])

    # # subtract -70 from 0
    peak_delay_df['peak_delayEI']   = (peak_delay_df[0]  - peak_delay_df[-70] )
    onset_delay_df['onset_delayEI'] = (onset_delay_df[0] - onset_delay_df[-70])

    # # rename -70 and 0 columns to exc_onset and inh_onset
    peak_delay_df.rename(columns={-70:'exc_peak', 0:'inh_peak'}, inplace=True)
    onset_delay_df.rename(columns={-70:'exc_onset', 0:'inh_onset'}, inplace=True)

    # # merge the two
    vc_delay_df = pd.merge(peak_delay_df, onset_delay_df, on=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','pulse'])
    # remove those rows where the onset delay is more than 20 or less than -20
    vc_delay_df = vc_delay_df[(vc_delay_df['onset_delayEI'] < 20) & (vc_delay_df['onset_delayEI'] > -20)]
    # multiply the delay by 1000 to convert to ms: following columns: 'exc_peak','inh_peak','peak_delayEI','exc_onset','inh_onset','onset_delayEI'
    columnsvc = ['exc_peak','inh_peak','peak_delayEI','exc_onset','inh_onset','onset_delayEI']
    vc_delay_df[columnsvc] = vc_delay_df[columnsvc] * 1000
    print(vc_delay_df.shape)
    return vc_delay_df


def generate_ebyi_df_norm(vc_shortdf):
    # normPSC_i = PSC_i / PSC_0 — ratio of two same-sign values, so always positive for both
    # clampPotentials. Keep rows where normPSC > 0 (drops noise artifacts).
    idvars  = ['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','trialID']
    valvars = [f'normPSC_{i}' for i in freq_sweep_pulses]
    df2 = vc_shortdf.melt(id_vars=idvars, value_vars=valvars, var_name='pulse', value_name='PSC')
    df2 = df2[df2['PSC'] > 0]
    df2['pulse'] = df2['pulse'].apply(lambda x: int(x.split('_')[-1]))
    df2 = df2.dropna(subset=['PSC'])
    df2.drop(columns=['trialID'], inplace=True)
    df2['numSq'] = df2['numSq'].astype('int')
    df4 = df2.groupby(['cellID','clampPotential','stimFreq','numSq','patternList','pulseWidth','intensity','pulse']).mean().reset_index()
    ebyi_df = df4.pivot(index=['cellID','stimFreq','numSq','patternList','pulseWidth','intensity','pulse'], columns='clampPotential', values='PSC').reset_index()
    ebyi_df = ebyi_df.dropna(subset=[-70,0])
    ebyi_df['EbyI'] = (ebyi_df[-70] / ebyi_df[0])
    ebyi_df = ebyi_df[(ebyi_df['EbyI'] < 20) & (ebyi_df['EbyI'] > 0)]
    print(ebyi_df.shape)
    return ebyi_df


# ── Compute or load intermediate dataframes ────────────────────────────────────

_sdn_file   = paper_figure_export_location / "Figure3_sdn_data_FS_LTM_v2.h5"
_fitdf_file = paper_figure_export_location / "Figure3_gamma_and_slope_fits_FS_LTM_v2.h5"

if _sdn_file.exists():
    sdn_df = pd.read_hdf(_sdn_file, key='data')
    print('Loaded sdn_df:', sdn_df.shape)
else:
    print('Computing sdn_df (FS patterns)...')
    sdn_data = []
    for cell in np.sort(cc_FS_shortdf_slice['cellID'].unique()):
        for patternID in patternIDs:
            for freq in [20, 30, 40, 50]:
                celldf = cc_FS_shortdf_slice[(cc_FS_shortdf_slice['cellID'] == cell)]
                for pulse_index in freq_sweep_pulses:
                    try:
                        x = calculate_expected_response(celldf, pulse_index, freq, patternID)
                        _numSq, _freq, _pid, _pidx, _fld, obs_r, obs_sc, exp_r = x
                        for obs, obs_sc_ in zip(obs_r, obs_sc):
                            sdn_data.append({'cellID': cell, 'numSq': _numSq, 'stimFreq': freq,
                                             'patternID': patternID, 'pulse': pulse_index,
                                             'obs': obs, 'obs_scaled': obs_sc_, 'exp': exp_r,
                                             'AP': 1 if obs > 20 else 0})
                    except ValueError as e:
                        print(e); continue

    # Add LTM patterns (7-sq and 15-sq) for cells that also have FS 1-sq reference data
    print('Computing sdn_df (LTM patterns)...')
    fs_cells_with_1sq = set(cc_FS_shortdf_slice[cc_FS_shortdf_slice['numSq']==1]['cellID'].unique())
    for cell in np.sort(cc_LTM_slice['cellID'].unique()):
        if cell not in fs_cells_with_1sq:
            continue
        fs_1sq = cc_FS_shortdf_slice[(cc_FS_shortdf_slice['cellID']==cell) & (cc_FS_shortdf_slice['numSq']==1)]
        ltm_msq = cc_LTM_slice[cc_LTM_slice['cellID']==cell]
        celldf_combined = pd.concat([fs_1sq, ltm_msq], axis=0, ignore_index=True)
        ltm_patternIDs = np.sort(ltm_msq['patternList'].unique())
        for patternID in ltm_patternIDs:
            for freq in [20, 30, 40, 50]:
                for pulse_index in range(8):  # LTM has 8 pulses (0-7)
                    try:
                        x = calculate_expected_response(celldf_combined, pulse_index, freq, patternID)
                        _numSq, _freq, _pid, _pidx, _fld, obs_r, obs_sc, exp_r = x
                        for obs, obs_sc_ in zip(obs_r, obs_sc):
                            sdn_data.append({'cellID': cell, 'numSq': _numSq, 'stimFreq': freq,
                                             'patternID': patternID, 'pulse': pulse_index,
                                             'obs': obs, 'obs_scaled': obs_sc_, 'exp': exp_r,
                                             'AP': 1 if obs > 20 else 0})
                    except ValueError as e:
                        print(e); continue

    sdn_df = pd.DataFrame(sdn_data)
    sdn_df.to_hdf(_sdn_file, key='data')
    print('Saved sdn_df:', sdn_df.shape)

if _fitdf_file.exists():
    fitdf = pd.read_hdf(_fitdf_file, key='data')
    print('Loaded fitdf:', fitdf.shape)
else:
    print('Computing fitdf (this is slow)...')
    fitdf_temp = gamma_distribution(sdn_df[sdn_df['AP']==0], fitdf=None,
                                    x='exp', y='obs',        first='cellID', second='pulse', third='stimFreq')
    fitdf      = gamma_distribution(sdn_df[sdn_df['AP']==0], fitdf=fitdf_temp,
                                    x='exp', y='obs_scaled', first='cellID', second='pulse', third='stimFreq')
    fitdf.to_hdf(_fitdf_file, key='data')
    print('Saved fitdf:', fitdf.shape)

cc_delay_df  = pd.concat([generate_cc_delay_df(cc_FS_shortdf_slice),
                          generate_cc_psp_ltm(cc_LTM_slice)],
                         ignore_index=True)
vc_delay_df  = generate_vc_delay_df(vc_FS_shortdf_slice)
ebyi_df_norm = generate_ebyi_df_norm(vc_FS_shortdf_slice)
ebyi_df_raw  = generate_ebyi_df(vc_FS_shortdf_slice)


def _exp_sat(x, a, tau, b):
    """Exponential-saturation curve: b + a*(1 - exp(-x/tau))."""
    return b + a * (1 - np.exp(-x / tau))


def _prepare_ebyi_norm_filt(ebyi_df_norm):
    """Filter to multi-square patterns and rename clamp-potential columns."""
    ebyi_df_norm_filt = ebyi_df_norm[(ebyi_df_norm['numSq'] > 1)].copy()
    ebyi_df_norm_filt.rename(columns={-70: 'Exc', 0: 'Inh'}, inplace=True)
    return ebyi_df_norm_filt


def _apply_global_formatting(ax3):
    """Apply uniform tick/label/spine formatting to every panel in ax3."""
    for a in ax3.keys():
        ax3[a].tick_params(axis='both', which='major', labelsize=16)
        ax3[a].set_xlabel(ax3[a].get_xlabel(), fontsize=16)
        ax3[a].set_ylabel(ax3[a].get_ylabel(), fontsize=16)
        ax3[a].spines['left'].set_linewidth(2)
        ax3[a].spines['bottom'].set_linewidth(2)


def panel_A(ax, cc_delay_df):
    """Panel A — normalised PSP amplitude across pulse indices (20 Hz, multi-sq)."""
    ax.text(-0.1, 1.10, 'A', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    cc_delay_df_filt = cc_delay_df[
        (cc_delay_df['numSq'] > 1) & (cc_delay_df['stimFreq'] == 20)
    ]
    sns.pointplot(data=cc_delay_df_filt, x='pulse', y='peak_PSP',
                  ax=ax, color=rocket_r(0.5), errorbar='ci')
    pulsewise_responses = cc_delay_df_filt.pivot_table(
        columns='pulse', index='trialID', values='peak_PSP')
    kruskal(*[pulsewise_responses[col] for col in pulsewise_responses.columns])
    ax.legend([], [], frameon=False)
    ax.set_ylabel('PSP (norm.)', fontsize=16)
    ax.set_xlabel('Pulse Index', fontsize=16)
    sns.despine(ax=ax, top=True, right=True, offset=10, trim=True)

    pm   = cc_delay_df_filt.groupby(['cellID', 'pulse'])['peak_PSP'].mean().reset_index()
    p2v  = pm[pm['pulse'] == 2]['peak_PSP'].dropna().values
    p8v  = pm[pm['pulse'] == 8]['peak_PSP'].dropna().values
    _, pval2 = wilcoxon(p2v - 1.0, alternative='greater')
    _, pval8 = wilcoxon(p8v - 1.0, alternative='less')
    print(f'Panel A  pulse 2 vs 1.0 (greater): p={pval2:.4f}')
    print(f'Panel A  pulse 8 vs 1.0 (less):    p={pval8:.4f}')
    ytop = ax.get_ylim()[1]
    for xi, pv in [(2, pval2), (8, pval8)]:
        star = ('***' if pv < 0.001 else '**' if pv < 0.01
                else '*' if pv < 0.05 else 'ns')
        ax.text(xi, ytop, star, ha='center', va='bottom',
                fontsize=14, color='dimgray')


def panel_B(ax, ebyi_df_norm_filt):
    """Panel B — Exc and Inh normalised PSC amplitudes across pulse indices."""
    ax.text(-0.1, 1.10, 'B', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    sns.pointplot(data=ebyi_df_norm_filt, x='pulse', y='Exc',
                  ax=ax, color=flare(0.5), errorbar='ci', label='Exc')
    sns.pointplot(data=ebyi_df_norm_filt, x='pulse', y='Inh',
                  ax=ax, color=crest(0.5), errorbar='ci', label='Inh')
    ax.legend(loc='upper right', ncols=1, fontsize='small', frameon=False)
    ax.set_ylabel('PSC Amplitude (norm.)', fontsize=16)
    ax.set_xlabel('Pulse Index', fontsize=16)
    sns.despine(ax=ax, top=True, right=True, offset=10, trim=True)

    # Holm-Bonferroni corrected Wilcoxon per pulse (skip pulse 0)
    pulses  = sorted(ebyi_df_norm_filt['pulse'].unique())
    pv_raw  = []
    for bp in pulses:
        if bp == 0:
            pv_raw.append(np.nan)
            continue
        bdf = ebyi_df_norm_filt[ebyi_df_norm_filt['pulse'] == bp]
        be  = bdf.groupby('cellID')['Exc'].mean()
        bi  = bdf.groupby('cellID')['Inh'].mean()
        cc  = be.index.intersection(bi.index)
        try:
            _, pv = (wilcoxon(be[cc].values, bi[cc].values)
                     if len(cc) >= 5 else (np.nan, np.nan))
        except ValueError:
            pv = np.nan
        pv_raw.append(pv)
    valid = [(i, p) for i, p in enumerate(pv_raw) if not np.isnan(p)]
    if valid:
        vi, vp = zip(*valid)
        _, pc, _, _ = multipletests(vp, method='holm')
        ytop = ax.get_ylim()[1]
        for xi, xp in zip(vi, pc):
            star = ('***' if xp < 0.001 else '**' if xp < 0.01
                    else '*' if xp < 0.05 else '')
            if star:
                ax.text(xi, ytop, star, ha='center', va='bottom',
                        fontsize=14, color='dimgray')


def panel_C(ax, ebyi_df_norm_filt):
    """Panel C — E/I ratio across pulse indices with exponential-saturation fit."""
    ax.text(-0.1, 1.10, 'C', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    sns.pointplot(data=ebyi_df_norm_filt, x='pulse', y='EbyI',
                  ax=ax, color=rocket_r(0.5), errorbar='ci')
    ax.legend([], [], frameon=False)
    ax.set_ylabel('E / I', fontsize=16)
    ax.set_xlabel('Pulse Index', fontsize=16)
    ax.set_ylim(top=3.0)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    sns.despine(ax=ax, top=True, right=True, offset=10, trim=True)

    pm = ebyi_df_norm_filt.groupby(['cellID', 'pulse'])['EbyI'].mean().reset_index()
    b0 = float(pm[pm['pulse'] == 0]['EbyI'].mean())
    p0 = [pm['EbyI'].max() - b0, 1.5, b0]
    try:
        popt, pcov = curve_fit(
            _exp_sat, pm['pulse'].values, pm['EbyI'].values,
            p0=p0, maxfev=10000,
            bounds=([0, 0.1, -np.inf], [np.inf, 20, np.inf]))
        tau_se = np.sqrt(np.diag(pcov))[1]
        xfit   = np.linspace(0, 8, 200)
        ax.plot(xfit, _exp_sat(xfit, *popt), color='black', lw=1.5, ls='--')
        print(f'Panel C  exp-sat fit: τ={popt[1]:.3f} ± {tau_se:.3f} pulses,  '
              f'a={popt[0]:.3f},  b={popt[2]:.3f}')
    except RuntimeError:
        print('Panel C: exponential saturation fit did not converge')


def panel_D(ax, vc_delay_df):
    """Panel D — Exc and Inh PSC onset delays across pulse indices."""
    ax.text(-0.1, 1.10, 'D', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    vc_delay_df_filt = vc_delay_df[(vc_delay_df['numSq'] > 1)]
    sns.pointplot(data=vc_delay_df_filt, x='pulse', y='exc_onset',
                  ax=ax, color=flare(0.5), errorbar='ci', label='Exc')
    sns.pointplot(data=vc_delay_df_filt, x='pulse', y='inh_onset',
                  ax=ax, color=crest(0.5), errorbar='ci', label='Inh')
    ax.set_xticks(np.arange(0, 9))
    ax.set_ylim(0, 10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
    ax.set_xlabel('Pulse Index', fontsize=16)
    ax.set_ylabel('PSC Onset Delay (ms)')
    ax.legend(loc='lower right', ncols=1, fontsize='small', frameon=False)
    [ax.spines[place].set_visible(False) for place in ['top', 'right']]
    sns.despine(ax=ax, offset=10, trim=True)

    exc_m = vc_delay_df_filt.groupby('pulse')['exc_onset'].mean().dropna()
    inh_m = vc_delay_df_filt.groupby('pulse')['inh_onset'].mean().dropna()
    se, ie, _, pe, _ = linregress(exc_m.index.values, exc_m.values)
    si, ii, _, pi, _ = linregress(inh_m.index.values, inh_m.values)
    ax.plot(exc_m.index.values, ie + se * exc_m.index.values,
            color=flare(0.8), lw=1.5, ls='--')
    ax.plot(inh_m.index.values, ii + si * inh_m.index.values,
            color=crest(0.8), lw=1.5, ls='--')

    cs_e, cs_i = [], []
    for cid in vc_delay_df_filt['cellID'].unique():
        cdf = vc_delay_df_filt[vc_delay_df_filt['cellID'] == cid]
        ce  = cdf.groupby('pulse')['exc_onset'].mean().dropna()
        ci  = cdf.groupby('pulse')['inh_onset'].mean().dropna()
        if len(ce) >= 3 and len(ci) >= 3:
            cs_e.append(linregress(ce.index.values, ce.values)[0])
            cs_i.append(linregress(ci.index.values, ci.values)[0])
    if len(cs_e) >= 5:
        _, p_comp = wilcoxon(np.array(cs_e), np.array(cs_i))
    else:
        p_comp = np.nan
    print(f'Panel D  Exc slope: {se:.3f} ms/pulse  (p={pe:.4f})')
    print(f'Panel D  Inh slope: {si:.3f} ms/pulse  (p={pi:.4f})')
    print(f'Panel D  slope diff (Exc vs Inh): p={p_comp:.4f}')


def panel_E(ax, sdn_df, fitdf, selected_cell):
    """Panel E — scatter + SDN/linear fits for a single example cell."""
    ax.text(-0.1, 1.10, 'E', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    fitdf_slice = fitdf[~fitdf['gamma'].isna()]
    dftemp = sdn_df[
        (sdn_df['cellID'] == selected_cell) &
        (sdn_df['pulse'] == 0) &
        (sdn_df['AP'] == 0)
    ]
    ax = sns.scatterplot(data=dftemp, x='exp', y='obs', hue='numSq',
                         size='numSq', sizes=[150], ax=ax,
                         palette=color_squares)
    gammatemp    = fitdf_slice[
        (fitdf_slice['cellID'] == selected_cell) &
        (fitdf_slice['observed'] == 'obs') &
        (fitdf_slice['pulse'] == 0) &
        (fitdf_slice['stimFreq'] == 1000)
    ]['gamma'].values
    slopetemp    = fitdf_slice[
        (fitdf_slice['cellID'] == selected_cell) &
        (fitdf_slice['observed'] == 'obs') &
        (fitdf_slice['pulse'] == 0) &
        (fitdf_slice['stimFreq'] == 1000)
    ]['slope'].values
    r2_gammatemp = fitdf_slice[
        (fitdf_slice['cellID'] == selected_cell) &
        (fitdf_slice['observed'] == 'obs') &
        (fitdf_slice['pulse'] == 0) &
        (fitdf_slice['stimFreq'] == 1000)
    ]['r2_sdn'].values
    r2_slopetemp = fitdf_slice[
        (fitdf_slice['cellID'] == selected_cell) &
        (fitdf_slice['observed'] == 'obs') &
        (fitdf_slice['pulse'] == 0) &
        (fitdf_slice['stimFreq'] == 1000)
    ]['r2_lin'].values
    print(gammatemp, slopetemp, r2_gammatemp, r2_slopetemp)
    xline = np.linspace(0, 20, 20)
    ax.plot(xline, sdnfunc(xline, gammatemp),
            color='purple', linewidth=3, label=f'γ = {gammatemp[0]:.2f}')
    ax.plot(xline, nosdn(xline, slopetemp),
            color='green', linewidth=3, label=f'm = {slopetemp[0]:.2f}')
    ax.plot([0, 15], [0, 15], color='grey', linestyle='--')
    ax.set_xlabel('Expected response (mV)')
    ax.set_ylabel('Observed response (mV)')
    handles_E, labels_E = ax.get_legend_handles_labels()
    ax.legend(handles_E[::-1], labels_E[::-1], loc='upper left',
              bbox_to_anchor=(0, 1.1), fontsize=16, frameon=False)
    ax.set_xlim([0, 15])
    ax.set_ylim([0, 10])
    sns.despine(bottom=False, left=False, trim=True, ax=ax)


def panel_F(ax, fitdf):
    """Panel F — cumulative histogram of gamma at pulse 0 vs pulse 8."""
    ax.text(-0.1, 1.10, 'F', fontweight='bold', fontsize=20,
            ha='center', transform=ax.transAxes)
    gammadist0 = fitdf[
        (fitdf['cellID'] != 1000) & (fitdf['pulse'] == 0) &
        (fitdf['stimFreq'] != 1000) & (fitdf['observed'] == 'obs') &
        (fitdf['sample_size'] != 0)
    ].dropna(subset=['gamma', 'slope'])
    gammadist8 = fitdf[
        (fitdf['cellID'] != 1000) & (fitdf['pulse'] == 8) &
        (fitdf['stimFreq'] != 1000) & (fitdf['observed'] == 'obs') &
        (fitdf['sample_size'] != 0)
    ].dropna(subset=['gamma', 'slope'])
    cap = 50
    gammadist0 = gammadist0.copy()
    gammadist8 = gammadist8.copy()
    gammadist0['gamma'] = gammadist0['gamma'].clip(upper=cap)
    gammadist8['gamma'] = gammadist8['gamma'].clip(upper=cap)
    gammadist0['pulse'] = 0
    gammadist8['pulse'] = 8
    gammadist = pd.concat([gammadist0, gammadist8], axis=0)

    sns.histplot(
        data=gammadist, x='gamma', hue='pulse',
        palette={0: 'purple', 8: 'pink'}, kde=True, ax=ax, alpha=0.5,
        line_kws={'lw': 2}, binwidth=1, element='step',
        cumulative=True, stat='density', common_norm=False, legend=True)
    sns.move_legend(ax, 'lower right', bbox_to_anchor=(1, 0),
                    title='Pulse', fontsize=16)
    ax.set_xlabel('Gamma (γ)', fontsize=16)
    ax.set_ylabel('Count', fontsize=16)
    sns.despine(bottom=False, left=False, ax=ax)
    ax.tick_params(axis='both', which='major', labelsize=16)

    _, pval_gamma = mannwhitneyu(gammadist0['slope'], gammadist8['slope'])
    ax.text(0.1, 0.9, f'p = {pval_gamma:.3f}',
            transform=ax.transAxes, fontsize=16, color='black')


def panel_G(ax, Fig3, fitdf):
    """Panel G — heatmap of median gamma by pulse × stimFreq.

    Returns the new axes object created by ax_to_partial_dist_heatmap_ax,
    which replaces the original axes in Fig3.
    """
    fitdf_slice = fitdf[
        (fitdf['cellID'] != 1000) & (fitdf['pulse'] != 1000) &
        (fitdf['stimFreq'] != 1000) & (fitdf['observed'] == 'obs') &
        (fitdf['sample_size'] != 0)
    ].dropna(subset=['gamma', 'slope']).copy()
    fitdf_slice.drop(columns=['expected', 'observed', 'cellID'], inplace=True)

    x = fitdf_slice.groupby(['pulse', 'stimFreq']).median().reset_index()
    n = fitdf_slice.groupby(['pulse', 'stimFreq']).count().reset_index()
    gammapivot   = x.pivot(index='stimFreq', columns='pulse', values='gamma')
    gammapivot_n = n.pivot(index='stimFreq', columns='pulse', values='gamma')

    new_ax, _, _, axc_G, cbar_G = ax_to_partial_dist_heatmap_ax(
        gammapivot, gammapivot_n, Fig3, ax,
        barw=0.03, pad=0.01, shrink=0.8,
        palette='Purples_r', annotate=False,
        show_marginals=False, cbar_label='Gamma')
    cbar_G.set_label('Gamma', fontsize=14)
    cbar_G.ax.tick_params(labelsize=14)
    new_ax.text(-0.1, 1.10, 'G', fontweight='bold', fontsize=20,
                ha='center', transform=new_ax.transAxes)
    return new_ax


def main():
    plt.close('all')
    Fig3, ax3 = plt.subplot_mosaic(
        [['A', 'B', 'C', 'D'], ['E', 'F', 'G', 'G']],
        figsize=(21, 10))
    plt.subplots_adjust(wspace=0.3, hspace=0.5)

    selected_cell = 3402
    ebyi_norm_filt = _prepare_ebyi_norm_filt(ebyi_df_norm)

    panel_A(ax3['A'], cc_delay_df)
    panel_B(ax3['B'], ebyi_norm_filt)
    panel_C(ax3['C'], ebyi_norm_filt)
    panel_D(ax3['D'], vc_delay_df)
    panel_E(ax3['E'], sdn_df, fitdf, selected_cell)
    panel_F(ax3['F'], fitdf)
    ax3['G'] = panel_G(ax3['G'], Fig3, fitdf)

    _apply_global_formatting(ax3)

    plt.show()

# make dataset
# sdn_df, fitdf, cc_delay_df, vc_delay_df, ebyi_df = make_dataset()


main()