import os
from pathlib import Path

import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt
import seaborn           as sns
import pandas            as pd
from scipy.stats import kruskal

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

color_E             = "rocket"
color_I             = "mako"
color_freq = {1:magma(0.05), 5:magma(0.1), 10:magma(0.2), 20:magma(.4), 30:magma(.5), 40:magma(.6), 50:magma(.7), 100:magma(.9)}
color_squares = {1:viridis(0.2), 5:viridis(.4), 7:viridis(.6), 15:viridis(.8), 20:viridis(1.0)}

freq_sweep_pulses = np.arange(9)

def expand_list_column(df, col, prefix):
    expanded = pd.DataFrame(df[col].tolist(), index=df.index)
    expanded.columns = [f'{prefix}{i}' for i in range(expanded.shape[1])]
    return pd.concat([df, expanded], axis=1)


# Paths
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
cc_LTM_slice = expand_list_column(cc_LTM_slice, 'peaks_cell_norm', 'normPSC_')
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


# ── Compute intermediate dataframes ───────────────────────────────────────────

cc_delay_df  = pd.concat([generate_cc_delay_df(cc_FS_shortdf_slice),
                          generate_cc_psp_ltm(cc_LTM_slice)],
                         ignore_index=True)
vc_delay_df  = generate_vc_delay_df(vc_FS_shortdf_slice)
ebyi_df_norm = generate_ebyi_df_norm(vc_FS_shortdf_slice)
ebyi_df_raw  = generate_ebyi_df(vc_FS_shortdf_slice)


def extended_figure2():
    Fig3, ax3 = plt.subplot_mosaic([['a','b','c'],['d','e','f'],['g','h','i'],['j','k','l'],['m','n','o'],['p','q','r']], figsize=(18,21),)
    plt.subplots_adjust(wspace=0.3, hspace=0.5)

    color_pulses_lin   = mpl.colormaps['Greens']
    color_pulses_gamma = mpl.colormaps['Purples']

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # ############ Voltage clamp plots #################
    # E by I ratio vs pulse index across frequencies
    ax3['a'].text(-0.1, 1.05, 'a', fontweight='bold', fontsize=16, ha='center', transform=ax3['a'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='peak_PSP', hue='stimFreq', ax=ax3['a'], palette=color_freq, errorbar='ci',)
    # run a kruskal wallis test across the pulsewise responses
    pulsewise_responses = cc_delay_df.pivot_table(columns='pulse', index='trialID', values='peak_PSP', )
    # now run KW test across columns
    kw_res = kruskal(*[pulsewise_responses[col] for col in pulsewise_responses.columns])
    ax3['a'].set_ylim([0, 2])
    ax3['a'].set_yticks([0,0.5,1.0,1.5,2.0])
    ax3['a'].set_ylabel('PSP (mV)', fontsize=12)
    ax3['a'].set_xlabel('Pulse Index', fontsize=12)
    ax3['a'].legend([],[], frameon=False)
    sns.despine(ax=ax3['a'], top=True, right=True, offset=10, trim=True)

    ax3['d'].text(-0.1, 1.05, 'd', fontweight='bold', fontsize=16, ha='center', transform=ax3['d'].transAxes)
    sns.pointplot(data=ebyi_df_raw, x='pulse', hue='stimFreq', y=-70, ax=ax3['d'], palette=color_E, errorbar='ci', )
    sns.pointplot(data=ebyi_df_raw, x='pulse', hue='stimFreq', y = 0, ax=ax3['d'], palette=color_I, errorbar='ci',  )
    # ax3['d'].set_ylim([0,1.5])
    # ax3['d'].set_yticks([0,0.5,1.0,1.5])
    ax3['d'].set_ylabel('PSC Amplitude (norm.)', fontsize=12)
    ax3['d'].set_xlabel('Pulse Index', fontsize=12)
    ax3['d'].legend([],[], frameon=False)
    sns.despine(ax=ax3['d'], top=True, right=True, offset=10, trim=True)

    ax3['g'].text(-0.1, 1.05, 'g', fontweight='bold', fontsize=16, ha='center', transform=ax3['g'].transAxes)
    sns.pointplot(data=ebyi_df_norm, x='pulse', y='EbyI', ax=ax3['g'], hue='stimFreq', palette=color_freq, errorbar='ci',)
    # ax3['g'].set_ylim([0, 5])
    ax3['g'].set_yticks(np.arange(0,5.1,1))
    ax3['g'].set_ylabel('E / I', fontsize=12)
    ax3['g'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['g'], top=True, right=True, offset=10, trim=True)
    # legend outside
    ax3['g'].legend( [],[], frameon=False)

    ### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ax3['b'].text(-0.1, 1.05, 'b', fontweight='bold', fontsize=16, ha='center', transform=ax3['b'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='onset_delay', ax=ax3['b'], hue='stimFreq', palette=color_freq, errorbar='ci',)
    ax3['b'].set_ylim([0, 15])
    ax3['b'].set_yticks([0,5,10,15])
    ax3['b'].set_ylabel('PSP Onset Delay (ms)', fontsize=12)
    ax3['b'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['b'], top=True, right=True, offset=10, trim=True)
    ax3['b'].legend([],[], frameon=False)

    ax3['e'].text(-0.1, 1.05, 'e', fontweight='bold', fontsize=16, ha='center', transform=ax3['e'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='exc_onset', hue='stimFreq', ax=ax3['e'], palette=color_E, errorbar='ci',)
    sns.pointplot(data=vc_delay_df, x='pulse', y='inh_onset', hue='stimFreq', ax=ax3['e'], palette=color_I, errorbar='ci',)
    ax3['e'].set_xticks( np.arange(0,9))
    ax3['e'].set_yticks( np.arange(0,16,5))
    ax3['e'].set_xlabel('Pulse Index', fontsize=12)
    ax3['e'].set_ylabel('PSC Onset Delay (ms)')
    [ax3['e'].spines[place].set_visible(False) for place in ['top', 'right', ] ]
    sns.despine(ax=ax3['e'], offset=10, trim=True)
    ax3['e'].legend([],[], frameon=False)

    ax3['h'].text(-0.1, 1.05, 'h', fontweight='bold', fontsize=16, ha='center', transform=ax3['h'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='onset_delayEI', hue='stimFreq', ax=ax3['h'], palette=color_freq, errorbar='ci',)
    ax3['h'].set_ylim([-2,6])
    ax3['h'].set_yticks([-2,0,2,4,6])
    ax3['h'].set_ylabel('Onset Delay (E-I) (ms)', fontsize=12)
    ax3['h'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['h'], top=True, right=True, offset=10, trim=True)
    ax3['h'].legend( [],[], frameon=False)


    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Voltage clamp plots
    # onset delay vs pulse index across frequencies
    # E by I ratio vs pulse index across frequencies
    ax3['c'].text(-0.1, 1.05, 'c', fontweight='bold', fontsize=16, ha='center', transform=ax3['c'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='peak_delay', hue='stimFreq', ax=ax3['c'], palette=color_freq, errorbar='ci',)
    ax3['c'].set_ylim([0, 30])
    ax3['c'].set_yticks(np.arange(0,31,5))
    ax3['c'].set_ylabel('PSP Peak Delay (ms)', fontsize=12)
    ax3['c'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['c'], top=True, right=True, offset=10, trim=True)
    ax3['c'].legend(loc='lower left', fontsize='small', bbox_to_anchor=(1.0, 0.0), title='stimulus frequency (Hz)')

    sns.pointplot(data=vc_delay_df, x='pulse', y='exc_peak', hue='stimFreq', ax=ax3['f'], palette=color_E, errorbar='ci',)
    sns.pointplot(data=vc_delay_df, x='pulse', y='inh_peak', hue='stimFreq', ax=ax3['f'], palette=color_I, errorbar='ci',)

    ax3['f'].text(-0.1, 1.05, 'f', fontweight='bold', fontsize=16, ha='center', transform=ax3['f'].transAxes)
    ax3['f'].legend(loc='lower left', ncols=2, fontsize='small', bbox_to_anchor=(1.0, 0.0), columnspacing=0.5, title='stimulus frequency (Hz) \n Exc vs Inh')
    # ax3['f'].set_ylim([0, 20])
    ax3['f'].set_xticks( np.arange(0,9))
    ax3['f'].set_yticks( np.arange(0,31,5))
    ax3['f'].set_xlabel('Pulse Index', fontsize=12)
    ax3['f'].set_ylabel('PSC Peak Delay (ms)')
    [ax3['f'].spines[place].set_visible(False) for place in ['top', 'right', ] ]
    sns.despine(ax=ax3['f'], offset=10, trim=True)

    ax3['i'].text(-0.1, 1.05, 'i', fontweight='bold', fontsize=16, ha='center', transform=ax3['i'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='peak_delayEI', ax=ax3['i'], hue='stimFreq', palette=color_freq, errorbar='ci',)
    ax3['i'].set_ylim([-2,6])
    ax3['i'].set_yticks([-2,0,2,4,6])
    ax3['i'].set_ylabel('Peak Delay (E-I) (ms)', fontsize=12)
    ax3['i'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['i'], top=True, right=True, offset=10, trim=True)
    ax3['i'].legend( loc='lower left', fontsize='small', bbox_to_anchor=(1.0, 0.0), title='stimulus frequency (Hz)')

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # numSq plots

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # ############ Voltage clamp plots #################
    # E by I ratio vs pulse index across frequencies
    ax3['j'].text(-0.1, 1.05, 'j', fontweight='bold', fontsize=16, ha='center', transform=ax3['j'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='peak_PSP', hue='numSq', ax=ax3['j'], palette=color_squares, errorbar='ci',)
    # run a kruskal wallis test across the pulsewise responses
    pulsewise_responses = cc_delay_df.pivot_table(columns='pulse', index='trialID', values='peak_PSP', )
    # now run KW test across columns
    kw_res = kruskal(*[pulsewise_responses[col] for col in pulsewise_responses.columns])
    ax3['j'].set_ylim([0, 2])
    ax3['j'].set_yticks([0,0.5,1.0,1.5,2.0])
    ax3['j'].set_ylabel('PSP (mV)', fontsize=12)
    ax3['j'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['j'], top=True, right=True, offset=10, trim=True)
    ax3['j'].legend([],[], frameon=False)

    ax3['m'].text(-0.1, 1.05, 'm', fontweight='bold', fontsize=16, ha='center', transform=ax3['m'].transAxes)
    sns.pointplot(data=ebyi_df_raw, x='pulse', hue='numSq', y=-70, ax=ax3['m'], palette=color_E, errorbar='ci', )
    sns.pointplot(data=ebyi_df_raw, x='pulse', hue='numSq', y = 0, ax=ax3['m'], palette=color_I, errorbar='ci',  )
    # ax3['m'].set_ylim([0,1.5])
    # ax3['m'].set_yticks([0,0.5,1.0,1.5])
    ax3['m'].legend([],[], frameon=False)
    ax3['m'].set_ylabel('PSC Amplitude (norm.)', fontsize=12)
    ax3['m'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['m'], top=True, right=True, offset=10, trim=True)

    ax3['p'].text(-0.1, 1.05, 'p', fontweight='bold', fontsize=16, ha='center', transform=ax3['p'].transAxes)
    sns.pointplot(data=ebyi_df_norm, x='pulse', y='EbyI', ax=ax3['p'], hue='numSq', palette=color_squares, errorbar='ci',)
    # ax3['p'].set_ylim([0, 5])
    ax3['p'].set_yticks(np.arange(0,5.1,1))
    ax3['p'].set_ylabel('E / I', fontsize=12)
    ax3['p'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['p'], top=True, right=True, offset=10, trim=True)
    ax3['p'].legend( [],[], frameon=False)

    ### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ax3['k'].text(-0.1, 1.05, 'k', fontweight='bold', fontsize=16, ha='center', transform=ax3['k'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='onset_delay', ax=ax3['k'], hue='numSq', palette=color_squares, errorbar='ci',)
    ax3['k'].set_ylim([0, 15])
    ax3['k'].set_yticks([0,5,10,15])
    ax3['k'].legend([],[], frameon=False)
    ax3['k'].set_ylabel('PSP Onset Delay (ms)', fontsize=12)
    ax3['k'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['k'], top=True, right=True, offset=10, trim=True)

    ax3['n'].text(-0.1, 1.05, 'n', fontweight='bold', fontsize=16, ha='center', transform=ax3['n'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='exc_onset', hue='numSq', ax=ax3['n'], palette=color_E, errorbar='ci',)
    sns.pointplot(data=vc_delay_df, x='pulse', y='inh_onset', hue='numSq', ax=ax3['n'], palette=color_I, errorbar='ci',)
    ax3['n'].set_xticks( np.arange(0,9))
    ax3['n'].set_yticks( np.arange(0,16,5))
    ax3['n'].set_xlabel('Pulse Index', fontsize=12)
    ax3['n'].set_ylabel('PSC Onset Delay (ms)')
    [ax3['n'].spines[place].set_visible(False) for place in ['top', 'right', ] ]
    sns.despine(ax=ax3['n'], offset=10, trim=True)
    ax3['n'].legend( [],[], frameon=False)

    ax3['q'].text(-0.1, 1.05, 'q', fontweight='bold', fontsize=16, ha='center', transform=ax3['q'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='onset_delayEI', hue='numSq', ax=ax3['q'], palette=color_squares, errorbar='ci',)
    ax3['q'].set_ylim([-6,6])
    ax3['q'].set_yticks([-6,-4,-2,0,2,4,6])
    ax3['q'].set_ylabel('Onset Delay (E-I) (ms)', fontsize=12)
    ax3['q'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['q'], top=True, right=True, offset=10, trim=True)
    ax3['q'].legend( [],[], frameon=False)


    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Voltage clamp plots
    # onset delay vs pulse index across frequencies
    # E by I ratio vs pulse index across frequencies
    ax3['l'].text(-0.1, 1.05, 'l', fontweight='bold', fontsize=16, ha='center', transform=ax3['l'].transAxes)
    sns.pointplot(data=cc_delay_df, x='pulse', y='peak_delay', hue='numSq', ax=ax3['l'], palette=color_squares, errorbar='ci',)
    ax3['l'].legend(loc='lower left', fontsize='small', bbox_to_anchor=(1.0, 0.0), title='Num Squares')
    ax3['l'].set_ylim([0, 30])
    ax3['l'].set_yticks(np.arange(0,31,5))
    ax3['l'].set_ylabel('PSP Peak Delay (ms)', fontsize=12)
    ax3['l'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['l'], top=True, right=True, offset=10, trim=True)

    sns.pointplot(data=vc_delay_df, x='pulse', y='exc_peak', hue='numSq', ax=ax3['o'], palette=color_E, errorbar='ci',)
    sns.pointplot(data=vc_delay_df, x='pulse', y='inh_peak', hue='numSq', ax=ax3['o'], palette=color_I, errorbar='ci',)

    ax3['o'].text(-0.1, 1.05, 'o', fontweight='bold', fontsize=16, ha='center', transform=ax3['o'].transAxes)
    ax3['o'].legend(loc='lower left', ncols=2, fontsize='small', bbox_to_anchor=(1.0, 0.0), title='Num Squares \n Exc vs Inh', columnspacing=0.5)
    # ax3['o'].set_ylim([0, 20])
    ax3['o'].set_xticks( np.arange(0,9))
    ax3['o'].set_yticks( np.arange(0,31,5))
    ax3['o'].set_xlabel('Pulse Index', fontsize=12)
    ax3['o'].set_ylabel('PSC Peak Delay (ms)')
    [ax3['o'].spines[place].set_visible(False) for place in ['top', 'right', ] ]
    sns.despine(ax=ax3['o'], offset=10, trim=True)

    ax3['r'].text(-0.1, 1.05, 'r', fontweight='bold', fontsize=16, ha='center', transform=ax3['r'].transAxes)
    sns.pointplot(data=vc_delay_df, x='pulse', y='peak_delayEI', ax=ax3['r'], hue='numSq', palette=color_squares, errorbar='ci',)
    ax3['r'].legend( loc='lower left', fontsize='small', bbox_to_anchor=(1.0, 0.0), title='Num Squares')
    # ax3['r'].set_ylim([-2,6])
    # ax3['r'].set_yticks([-2,0,2,4,6])
    ax3['r'].set_ylabel('Peak Delay (E-I) (ms)', fontsize=12)
    ax3['r'].set_xlabel('Pulse Index', fontsize=12)
    sns.despine(ax=ax3['r'], top=True, right=True, offset=10, trim=True)


    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    for a in ax3.keys():
        ax3[a].tick_params(axis='both', which='major', labelsize=12)
        # axis label fontsize
        ax3[a].set_xlabel(ax3[a].get_xlabel(), fontsize=12)
        ax3[a].set_ylabel(ax3[a].get_ylabel(), fontsize=12)
        # spine width
        ax3[a].spines['left'].set_linewidth(1)
        ax3[a].spines['bottom'].set_linewidth(1)

    # text label for part A and B
    ax3['a'].text(-0.2, 1.1, 'A', fontweight='bold', fontsize=20, ha='center', transform=ax3['a'].transAxes)
    ax3['j'].text(-0.2, 1.1, 'B', fontweight='bold', fontsize=20, ha='center', transform=ax3['j'].transAxes)

    for a in ['j','k','l','m','n','o','p','q','r']:
        pos = ax3[a].get_position()
        # shift the entire plot to down
        new_pos = [pos.x0, pos.y0-0.075, pos.width, pos.height]
        ax3[a].set_position(new_pos)


    plt.show()


extended_figure2()
