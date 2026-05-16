# Short-term dynamics of Excitation-Inhibition Balance in Hippocampal CA3-CA1 circuit
# Aditya Asopa, Upinder Singh Bhalla, NCBS
# Figure 2
# September 2024
# Imports -----------------------------------------------------------------------------------------------
from   pathlib      import Path

import numpy                as np
import matplotlib           as mpl
import matplotlib.pyplot    as plt
import seaborn              as sns
import pandas               as pd

import scipy.optimize as sci
import math

# ---------------------------------------------------------------------------
# Deconvolution helpers (from original plotFig2 / foo.py, U. S. Bhalla)
# ---------------------------------------------------------------------------
sampleRate = 20000.0

def tauFit(kernel, baseline):
    y = kernel[int(round(0.25 * len(kernel))):]
    pk = y[0]
    x = np.linspace(0, len(y) / sampleRate, len(y), endpoint=False)
    ret, cov = sci.curve_fit(lambda t, a, tau: a * np.exp(-t / tau), x, y, p0=(pk - baseline, 0.02))
    return ret

def calcKernel(dat):
    startIdx = int(round(startT * sampleRate))
    endIdx   = int(round(0.8 * endT * sampleRate))
    baseline = np.mean(dat.iloc[startIdx - int(0.005 * sampleRate):startIdx])
    rawKernel = np.array(dat.iloc[startIdx:endIdx])
    try:
        kmax = max(rawKernel)
        kmin = min(rawKernel)
    except:
        raise FloatingPointError("calcKernel: kmax or kmin is a nan")
    if math.isnan(kmax) or math.isnan(kmin):
        raise FloatingPointError("calcKernel: kmax or kmin is a nan")
    if abs(kmax) > abs(kmin):
        return kmax, rawKernel, baseline, tauFit(rawKernel, baseline)
    else:
        return kmin, rawKernel, baseline, tauFit(rawKernel, baseline)

def findStpScale(kernel, kpk, ret, si, stimWidth, tau, ax):
    if ret[0] < endT and si < (endT * sampleRate):
        return 1.0
    if kpk < 0:
        kpkIdx = np.argmin(kernel[:-stimWidth])
    else:
        kpkIdx = np.argmax(kernel[:-stimWidth])
    riseIdx   = int(round((ret[2] - ret[0]) * sampleRate))
    riseDelta = ret[1] - ret[1] * np.exp(-(ret[2] - ret[0]) / tau[1])
    if ax:
        label = "Min to Max" if (si < 11000 and kpk > 0) else None
        ax.plot([ret[2], ret[2]], [-riseDelta + ret[1], ret[3]], "ro-", label=label)
    if ret[0] < endT + 0.01:
        riseTotal = ret[3] - ret[1]
    else:
        riseTotal = riseDelta + ret[3] - ret[1]
    return riseTotal / kpk

def findPkVal(dat, freq, startIdx, isExc):
    stimWidth = int(round(0.7 * sampleRate / freq))
    d2 = np.array(dat.iloc[startIdx:startIdx + stimWidth])
    if isExc:
        imin = np.argmin(d2)
        d3   = np.array(dat.iloc[startIdx + imin - stimWidth:imin + startIdx])
        imax = np.argmax(d3)
        return [(imax + startIdx + imin - stimWidth) / sampleRate, d3[imax],
                (startIdx + imin) / sampleRate, d2[imin]]
    else:
        imax = np.argmax(d2)
        d3   = np.array(dat.iloc[startIdx + imax - stimWidth:imax + startIdx])
        try:
            imin = np.argmin(d3)
        except:
            raise FloatingPointError("findPkVal: imin is a nan")
        return [(imin + startIdx + imax - stimWidth) / sampleRate, d3[imin],
                (startIdx + imax) / sampleRate, d2[imax]]

def plotFromKernel(scaleList, stimIdx, kernel, freq, npv, label, ax):
    ret = np.zeros(int(round(sampleRate * 1.5)), dtype=np.float64)
    ret[int(round(sampleRate * endT)):] += npv[1][1]
    for ii in range(len(scaleList)):
        ss  = scaleList[ii]
        idx = stimIdx[ii]
        if idx > 0:
            ks = kernel * ss
            if label == "Inh":
                offset = npv[3, ii] - max(ks + ret[idx:len(kernel) + idx])
            else:
                offset = npv[3, ii] - min(ks + ret[idx:len(kernel) + idx])
            ret[idx:len(kernel) + idx] = ret[idx:len(kernel) + idx] + ks + offset
    t = np.arange(0.0, 1.0 - 1e-6, 1.0 / sampleRate)
    if ax:
        el1 = None if label == "Inh" else "Troughs"
        el2 = None if label == "Inh" else "Peaks"
        ax.plot(npv[0], npv[1], "c*-", label=el1)
        ax.plot(npv[2], npv[3], "y.-", label=el2)
    return ret[:len(t)]

def deconv(dat, freq, start_time, end_time, ax, noprobepulse=False):
    global startT, endT
    startT = start_time
    endT   = end_time

    stimWidth = int(round(sampleRate / freq))
    stimIdx   = [int(startT * sampleRate)] + [int(round(sampleRate * (endT + i / freq))) for i in range(8)]

    if noprobepulse:
        stimIdx = [int(startT * sampleRate)] + [int(round(sampleRate * (endT + i / freq))) for i in range(9)]
        stimIdx = stimIdx[1:]
        startT  = stimIdx[0] / 2e4
        endT    = stimIdx[1] / 2e4

    kpk, kernel, baseline, tau = calcKernel(dat)
    kpkidx = np.argmax(kernel) if kpk > 0 else np.argmin(kernel)

    pv    = []
    scaleList = []
    for si in stimIdx:
        ret   = findPkVal(dat, freq, si + kpkidx // 2, (kpk < 0))
        pv.append(ret)
        scale = findStpScale(kernel, kpk, ret, si, stimWidth, tau, None)
        scaleList.append(scale)
    label     = "Inh" if kpk > 0 else "Exc"
    npv       = np.array(pv).transpose()
    synthPlot = plotFromKernel(scaleList, stimIdx, kernel, freq, npv, label, None)
    return np.array(scaleList), synthPlot, npv, stimIdx
# ---------------------------------------------------------------------------

plt.rcParams['font.size'] = 14
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['svg.fonttype'] = 'none'

flare = mpl.colormaps["flare"]
crest = mpl.colormaps["crest"]


# Load data -----------------------------------------------------------------------------------------------
data_path = Path("../../DATA_Nov2025")

cc_FS_shortdf = pd.read_hdf(data_path / "all_cells_FreqSweep_CC_kernelfit_response_measurements.h5", key='data')
print(cc_FS_shortdf.shape)

vc_FS_shortdf = pd.read_hdf(data_path / "all_cells_FreqSweep_VC_kernelfit_response_measurements.h5", key='data')
print(vc_FS_shortdf.shape)

# CC data screening based on dataflag_fields
cc_FS_shortdf_slice = cc_FS_shortdf[
            (cc_FS_shortdf['location'] == 'CA1') &
            (cc_FS_shortdf['numSq'].isin([1,5,15])) &
            (cc_FS_shortdf['stimFreq'].isin([20,30,40,50])) &
            (cc_FS_shortdf['condition'] == 'Control') &
            (cc_FS_shortdf['ch0_response']==1) &
            (cc_FS_shortdf['IR'] >50) & (cc_FS_shortdf['IR'] < 300) &
            (cc_FS_shortdf['tau'] < 40) &
            (cc_FS_shortdf['intensity'] == 100) &
            (cc_FS_shortdf['pulseWidth'] == 2) &
            (cc_FS_shortdf['spike_in_baseline_period'] == 0) &
            (cc_FS_shortdf['ac_noise_power_in_ch0'] < 40)
        ]
print(cc_FS_shortdf.shape, '--screened-->', cc_FS_shortdf_slice.shape)
screened_cc_trialIDs = cc_FS_shortdf_slice['trialID'].unique()

print(f"Unique cells in screened data: { cc_FS_shortdf_slice['cellID'].nunique()}")
print(f"Unique sweeps in screened data: {cc_FS_shortdf_slice['trialID'].nunique()}")

np.savetxt(data_path / "Figure2_screened_trialIDs_CC_FS.txt", screened_cc_trialIDs, fmt='%s')

# VC data screening based on dataflag_fields
vc_FS_shortdf_slice = vc_FS_shortdf[
            (vc_FS_shortdf['location'] == 'CA1') &
            (vc_FS_shortdf['numSq'].isin([1,5,15])) &
            (vc_FS_shortdf['stimFreq'].isin([20,30,40,50])) &
            (vc_FS_shortdf['condition'] == 'Control') &
            (vc_FS_shortdf['ch0_response']==1) &
            (vc_FS_shortdf['intensity'] == 100) &
            (vc_FS_shortdf['pulseWidth'] == 2) &
            (vc_FS_shortdf['probePulseStart']==0.2) &
            (vc_FS_shortdf['IR'] >50) & (vc_FS_shortdf['IR'] < 300) &
            (vc_FS_shortdf['tau'] < 40) &
            (vc_FS_shortdf['ac_noise_power_in_ch0'] < 40)&
            (vc_FS_shortdf['valley_0'].notnull())
        ]
print(vc_FS_shortdf.shape, '--screened-->', vc_FS_shortdf_slice.shape)
screened_vc_trialIDs = vc_FS_shortdf_slice['trialID'].unique()

print(f"Unique cells in screened data: { vc_FS_shortdf_slice['cellID'].nunique()}")
print(f"Unique sweeps in screened data: {vc_FS_shortdf_slice['trialID'].nunique()}")

np.savetxt(data_path / "Figure2_screened_trialIDs_VC_FS.txt", screened_vc_trialIDs, fmt='%s')

# combine short dataframes slice and delete the original ones
xc_FS_shortdf_slice = pd.concat([cc_FS_shortdf_slice, vc_FS_shortdf_slice], axis=0)
del cc_FS_shortdf, vc_FS_shortdf
del cc_FS_shortdf_slice, vc_FS_shortdf_slice

### Load the Longform data and keep the screened trials only to save space
vc_FS_datapath =  data_path / "all_cells_FreqSweep_VC_long.h5"
vc_FS_longdf = pd.read_hdf(vc_FS_datapath, key='data')

vc_FS_longdf_slice = vc_FS_longdf[ vc_FS_longdf['trialID'].isin(screened_vc_trialIDs) ]
print('VC: ', vc_FS_longdf.shape, '--screened-->', vc_FS_longdf_slice.shape)
del vc_FS_longdf

# ---------------------------------------------------------------------------
# Heatmap helpers (from plot_tools.py)
# ---------------------------------------------------------------------------
def ax_to_partial_dist_heatmap_ax(pivotdf, numdf, fig, ax, barw=0.03, pad=0.01, shrink=0.8, palette='viridis', annotate=False, show_marginals=True, cbar_label='', y_offset=0.0):
    bboxA = ax.get_position()
    x0,x1 = bboxA.x0,bboxA.x1
    y0,y1 = bboxA.y0 + y_offset, bboxA.y1 + y_offset
    w, h  = bboxA.width,bboxA.height

    ax.remove()

    if show_marginals:
        ax  = fig.add_axes([x0, y0, shrink*w, shrink*h])
        axx = fig.add_axes([x0, y0+shrink*h+pad, shrink*w, barw], aspect='auto')
        axy = fig.add_axes([x0+shrink*w+pad, y0, barw, shrink*h], aspect='auto')
        axc = fig.add_axes([x0+shrink*w+barw+2*pad, y0, barw, shrink*h], aspect='auto')
    else:
        main_w = w - barw - 2*pad
        ax  = fig.add_axes([x0, y0, main_w, h])
        axx = None
        axy = None
        axc = fig.add_axes([x0+main_w+pad, y0, barw, h], aspect='auto')

    partial_pulse_wise   = pivotdf.mean(axis=0).values.reshape(1,-1)
    partial_freq_wise    = pivotdf.mean(axis=1).values.reshape(-1,1)
    partial_pulse_wise_n = numdf.sum(axis=0).values.reshape(1,-1)
    partial_freq_wise_n  = numdf.sum(axis=1).values.reshape(-1,1)
    maxlim = np.round(np.max(pivotdf.values),2)
    minlim = np.round(np.min(pivotdf.values),2)

    ax.imshow(pivotdf, cmap=palette, vmin=minlim, vmax=maxlim, aspect='auto')
    if show_marginals:
        axx.imshow(partial_pulse_wise, cmap=palette, vmin=minlim, vmax=maxlim)
        axy.imshow(partial_freq_wise,  cmap=palette, vmin=minlim, vmax=maxlim, origin='lower')

    if annotate:
        for i in range(partial_freq_wise_n.shape[0]):
            if show_marginals:
                axy.text(0, i, f'{partial_freq_wise[i,0]:.2f}', ha='center', va='center', color='white', fontsize=12)
                axy.text(0-0.2, i-0.2, f'{partial_freq_wise_n[i,0]:.0f}', ha='center', va='center', color='yellow', fontsize=10)
            for j in range(partial_pulse_wise_n.shape[1]):
                ax.text(j, i, f'{pivotdf.values[i,j]:.2f}', ha='center', va='center', color='white', fontsize=12)
                ax.text(j-0.2, i-0.2, f'{numdf.values[i,j]:.0f}', ha='center', va='center', color='yellow', fontsize=10)
                if i == 0 and show_marginals:
                    axx.text(j, 0, f'{partial_pulse_wise[0,j]:.2f}', ha='center', va='center', color='white', fontsize=12)
                    axx.text(j-0.2, 0-0.2, f'{partial_pulse_wise_n[0,j]:.0f}', ha='center', va='center', color='yellow', fontsize=10)

    ax.set_xticks(np.arange(9), labels=np.arange(9))
    ax.set_ylim([-0.5,3.5])
    ax.set_yticks([0,1,2,3], labels=[20,30,40,50])
    ax.set_xlabel('Pulse Index')
    ax.set_ylabel('Frequency (Hz)')
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(ax.get_images()[0], cax=axc)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=14)

    if show_marginals:
        axx.xaxis.set_tick_params(labelbottom=False)
        axy.yaxis.set_tick_params(labelleft=False)
        axx.get_xaxis().set_visible(False)
        axx.get_yaxis().set_visible(False)
        axy.get_xaxis().set_visible(False)
        axy.get_yaxis().set_visible(False)
        for ax_ in [axx, axy, axc]:
            for spine in ax_.spines.values():
                spine.set_visible(False)
        axx.set_aspect('auto')
        axy.set_aspect('auto')
    else:
        for spine in axc.spines.values():
            spine.set_visible(False)

    return ax, axx, axy, axc, cbar

def plot_response_heatmaps(datadf, feature='PSC', skip1sq=True, include_spike_trials=False, Fig=None, figlabels=[], clampMode='VC', heatmap_title=True, annot=False, show_marginals=True, cbar_label='', row_y_offsets=None):
    if feature == 'spike_':
        include_spike_trials = True
        datadf = datadf[datadf['AP'] == 1]

    if include_spike_trials == False and clampMode == 'CC':
        datadf = datadf[datadf['spike_in_stim_period'] == 0]

    print('data shape:', datadf.shape)

    freq_sweep_pulses = range(9)
    to_plot = [f'{feature}{i}' for i in freq_sweep_pulses]
    df_melt = pd.melt(datadf, id_vars=['cellID', 'clampPotential','stimFreq','numSq','patternList'], value_vars=to_plot, var_name='pulseIndex', value_name='peak_response')
    df_melt['pulse'] = df_melt.apply(lambda x: x['pulseIndex'][-1], axis=1)
    df_melt['pulse'] = df_melt['pulse'].astype(int)
    df_melt['numSq'] = df_melt['numSq'].astype(int)
    df_melt['clampPotential'] = df_melt['clampPotential'].astype(int)
    df_melt['stimFreq'] = df_melt['stimFreq'].astype(int)
    df_melt['patternList'] = df_melt['patternList'].apply(lambda x: int(x))
    df_melt.drop(columns=['pulseIndex'], inplace=True)
    df_melt['peak_response'] = df_melt['peak_response'].abs()

    sqs = np.sort(df_melt['numSq'].unique())
    if skip1sq:
        sqs = np.delete(sqs, np.where(sqs == 1))
    clamps = df_melt['clampPotential'].unique()

    colors_EI = {-70: flare, 0: crest}
    palette = {-70: 'viridis'} if clampMode == 'CC' else colors_EI

    if Fig is not None:
        Fig.clear()
        ax2 = Fig.subplots(len(sqs), len(clamps), sharex=False, sharey=False)
    else:
        Fig, ax2 = plt.subplots(len(sqs), len(clamps), figsize=(15,10), sharex=False, sharey=False)

    ax2 = np.array(ax2).reshape(len(sqs), len(clamps))
    Fig.subplots_adjust(hspace=0.5, wspace=0.5)
    if not figlabels:
        figlabels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    assert len(figlabels) >= len(sqs)*len(clamps)

    axes = []
    counter = 0
    for s,sq in enumerate(sqs):
        for c,clamp in enumerate(clamps):
            xpscdf = df_melt[(df_melt['numSq'] == sq) & (df_melt['clampPotential'] == clamp)]
            x = xpscdf.groupby(['pulse', 'stimFreq']).mean().reset_index()
            x_matrix = x.pivot(index='stimFreq', columns='pulse', values='peak_response')
            num_trials = xpscdf.groupby(['pulse', 'stimFreq']).count().reset_index()
            num_trials_matrix = num_trials.pivot(index='stimFreq', columns='pulse', values='peak_response')

            if x_matrix.shape[0] == 0:
                x_matrix = pd.DataFrame(np.nan, index=[20,30,40,50], columns=np.arange(9))
                num_trials_matrix = pd.DataFrame(0, index=[20,30,40,50], columns=np.arange(9))
            for freq in [20,30,40,50]:
                if freq not in x_matrix.index:
                    x_matrix.loc[freq] = np.nan
                if freq not in num_trials_matrix.index:
                    num_trials_matrix.loc[freq] = 0

            x_matrix = x_matrix.sort_index()
            num_trials_matrix = num_trials_matrix.sort_index()

            y_off = row_y_offsets[s] if (row_y_offsets is not None and s < len(row_y_offsets)) else 0.0
            axs = ax_to_partial_dist_heatmap_ax(x_matrix, num_trials_matrix, Fig, ax2[s,c], barw=0.03, pad=0.01, shrink=0.8, palette=palette[clamp], annotate=annot, show_marginals=show_marginals, cbar_label=cbar_label, y_offset=y_off)
            axs[0].text(-0.1, 1.1, figlabels[counter], transform=axs[0].transAxes, size=20, weight='bold')
            if heatmap_title:
                axs[0].text(0.05, 1.04, f'{sq} Sq', transform=axs[0].transAxes, fontsize=15)
            axes.append(axs)
            counter += 1

    return Fig, axes

# ---------------------------------------------------------------------------

def main():
    # Setup the figure
    w,h = 15,25
    fig = plt.figure(layout='constrained', figsize=(w,h))

    [Fig2Top, Fig2Mid, Fig2Bottom] = fig.subfigures(3,1, wspace=0.03, hspace=0.02, height_ratios=[2, 1.2, 2])

    [subfigsA, subfigsB] = Fig2Top.subfigures(1,2)
    [subfigsE, subfigsF] = Fig2Bottom.subfigures(1,2)

    [ax2C, ax2D] = Fig2Mid.subplots(1, 2)

    ## -----------------------------------------------------------------------------------------------
    # Fig2A: CC heatmap of normPSC
    dftemp = xc_FS_shortdf_slice[(xc_FS_shortdf_slice['clampMode']=='CC')]
    f,a = plot_response_heatmaps(dftemp[dftemp['AP']==0], feature='normPSC_', Fig=subfigsA, figlabels=['Ai','Aii'], clampMode='CC', annot=False, show_marginals=False, cbar_label='Norm. PSP', row_y_offsets=[0, -0.08])
    for axs in a:
        axs[0].set_ylabel('Frequency (Hz)', labelpad=15)

    # spike likelihood heatmap
    dftemp = xc_FS_shortdf_slice[(xc_FS_shortdf_slice['clampMode']=='CC')]
    dftemp['numspikes'] = dftemp[[f'spike_{i}' for i in range(9)]].sum(axis=1)
    dftemp= dftemp[(dftemp['AP']==1) ]
    f,a = plot_response_heatmaps(dftemp, feature='spike_', Fig=subfigsB, figlabels=['Bi','Bii'], clampMode='CC', annot=False, show_marginals=False, cbar_label='Spike Probability', row_y_offsets=[0, -0.08])

    ## -----------------------------------------------------------------------------------------------
    # Fig2C: kernel fit example for both E and I
    cell = 7492
    pattern = 52
    trial = 0
    exc_sweep = vc_FS_longdf_slice[(vc_FS_longdf_slice['cellID']==cell) & (vc_FS_longdf_slice['patternList']==pattern) & (vc_FS_longdf_slice['clampPotential']==-70)& (vc_FS_longdf_slice['stimFreq']==20)]
    inh_sweep = vc_FS_longdf_slice[(vc_FS_longdf_slice['cellID']==cell) & (vc_FS_longdf_slice['patternList']==pattern) & (vc_FS_longdf_slice['clampPotential']==0)  & (vc_FS_longdf_slice['stimFreq']==20)]

    row = exc_sweep.iloc[0, :]
    _,_, npv_exc, _ = exc_results = deconv(row[49:80049], row['stimFreq'], row['probePulseStart'], row['pulseTrainStart'], None, noprobepulse=(row['probePulseStart']==0.5))
    ax2C.plot(np.linspace(0,1,20000), row[49:20049], color='#e46d5dff', label='PSC Exc')
    ax2C.plot(npv_exc[0], npv_exc[1], color='#ecab7d9d', marker=".", linewidth=2, label='Peak Exc')
    ax2C.plot(npv_exc[2], npv_exc[3], color='#c23f69bf', marker="*", linewidth=2, label='Valley Exc')

    row = inh_sweep.iloc[0, :]
    _,_, npv_inh, _ = inh_results = deconv(row[49:80049], row['stimFreq'], row['probePulseStart'], row['pulseTrainStart'], None, noprobepulse=(row['probePulseStart']==0.5))
    ax2C.plot(np.linspace(0,1,20000), row[49:20049], color='#2f818dff', label='PSC Inh')
    ax2C.plot(npv_inh[0], npv_inh[1], color='#9dca929a', marker="^", linewidth=2, label='Valley Inh')
    ax2C.plot(npv_inh[2], npv_inh[3], color='#2c3071b3', marker="s", linewidth=2, label='Peak Inh')

    ax2C.set_xlabel('Time (s)')
    ax2C.set_ylabel('PSC (pA)')
    ax2C.set_xlim([0,1])
    ax2C.set_ylim([-500,1500])
    ax2C.yaxis.set_major_locator(plt.MultipleLocator(500))
    handles, labels = ax2C.get_legend_handles_labels()
    # PSC E, PE, VE, PSCI, VI, PI, 
    handles = [handles[i] for i in [3, 0, 4, 2, 5, 1]]
    labels  = [labels[i]  for i in [3, 0, 4, 2, 5, 1]]
    ax2C.legend(handles, labels, loc='upper left', fontsize=14, ncols=3, frameon=False)
    ax2C.text(-0.1, 1.05, 'C', fontsize=18, fontweight='bold', transform=ax2C.transAxes)
    sns.despine(ax=ax2C, top=True, right=True, trim=True, offset=10)

    ## -----------------------------------------------------------------------------------------------
    # Fig2D: E and I kernel fits
    ax2D.plot(range(9), exc_results[0], linewidth=2, color='#e46d5dff', label='Exc')
    ax2D.plot(range(9), inh_results[0], linewidth=2, color='#2f818dff', label='Inh')
    ax2D.set_xticks(range(9))
    ax2D.set_ylim([0, 1.5])
    ax2D.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax2D.set_xlabel('Pulse #')
    ax2D.set_ylabel('Normalized PSC')
    ax2D.legend(loc='upper right', fontsize=14, frameon=False)
    ax2D.text(-0.1, 1.05, 'D', fontsize=18, fontweight='bold', transform=ax2D.transAxes)
    sns.despine(ax=ax2D, top=True, right=True, trim=True, offset=10)

    ## -----------------------------------------------------------------------------------------------
    # Fig2E: VC heatmap of normPSC, excitatory (-70 mV)
    dftemp = xc_FS_shortdf_slice[(xc_FS_shortdf_slice['clampMode']=='VC')]
    f,a = plot_response_heatmaps(dftemp[dftemp['clampPotential']==-70], feature='normPSC_', Fig=subfigsE, figlabels=['Ei','Eii'], clampMode='VC', annot=False, show_marginals=False, cbar_label='Norm. EPSC', row_y_offsets=[0.08, 0])
    for axs in a:
        axs[0].set_ylabel('Frequency (Hz)', labelpad=15)

    # Fig2F: VC heatmap of normPSC, inhibitory (0 mV)
    dftemp = xc_FS_shortdf_slice[(xc_FS_shortdf_slice['clampMode']=='VC')]
    f,a = plot_response_heatmaps(dftemp[dftemp['clampPotential']==0], feature='normPSC_', Fig=subfigsF, figlabels=['Fi','Fii'], clampMode='VC', annot=False, show_marginals=False, cbar_label='Norm. IPSC', row_y_offsets=[0.08, 0])

    plt.show()

if __name__ == "__main__":
    main()
