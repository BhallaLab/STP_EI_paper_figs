import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import pandas
import copy
import numpy as np
import math
import argparse
import scipy
from scipy.stats import linregress

freq = 80.0 # Hz
width = 0.002

SAMPLE_FREQ = 20000
SAMPLE_TIME = 5
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16

def rippleSpikeRate( dcell, spikeCriterion = -0.03, windowSize = 0.02 ):
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_OFF = int( 1.3 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_NUM_SAMPLES = STIM_OFF - STIM_ON
    df = dcell.loc[(dcell['stimFreq'] == 100)]
    time = np.array( np.arange( STIM_NUM_SAMPLES, dtype = float) ) / SAMPLE_FREQ
    sumTrain = np.zeros( STIM_NUM_SAMPLES )
    binSize = int( 0.01 * SAMPLE_FREQ ) # 10 ms bins for this.
    n_bins = STIM_NUM_SAMPLES // binSize
    trialBins = np.zeros( (len(df), n_bins), dtype=float )
    for sweep in range( len(df) ):
        epsp = np.array(df.iloc[sweep, STIM_ON:STIM_OFF ])
        spikeTrain = np.zeros_like(epsp, dtype=float)
        crossings = np.where(np.diff(epsp > spikeCriterion) > 0)
        spikeTrain[crossings] = 1
        sumTrain += spikeTrain
        trialBins[sweep] = np.sum(spikeTrain.reshape(-1, binSize), axis=1)

    # Convolve the summed spike trains with a rectangular window
    window = np.ones(int(windowSize * SAMPLE_FREQ ))
    firingRate = np.convolve(sumTrain, window, mode='same') / (windowSize * len( df ) )

    return time, firingRate, np.sum(sumTrain.reshape(-1, binSize), axis=1), trialBins

def spikeRate( dcell, spikeCriterion = -0.03, windowSize = 0.02 ):
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_OFF = int( 1.3 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_NUM_SAMPLES = STIM_OFF - STIM_ON
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    time = np.array( np.arange( STIM_NUM_SAMPLES, dtype = float) ) / SAMPLE_FREQ
    sumTrain = np.zeros( STIM_NUM_SAMPLES )
    binSize = 400   # 20 ms * 20KHz
    n_bins = STIM_NUM_SAMPLES // binSize
    trialBins = np.zeros( (len(df), n_bins), dtype=float )
    for sweep in range( len(df) ):
        epsp = np.array(df.iloc[sweep, STIM_ON:STIM_OFF ])
        spikeTrain = np.zeros_like(epsp, dtype=float)
        crossings = np.where(np.diff(epsp > spikeCriterion) > 0)
        spikeTrain[crossings] = 1
        sumTrain += spikeTrain
        trialBins[sweep] = np.sum(spikeTrain.reshape(-1, binSize), axis=1)

    # Convolve the summed spike trains with a rectangular window to plot
    window = np.ones(int(windowSize * SAMPLE_FREQ ))
    firingRate = np.convolve(sumTrain, window, mode='same') / (windowSize * len( df ) )

    # Bin the spikes into 20 ms bins, to match the 50Hz input train.
    return time, firingRate, np.sum(sumTrain.reshape(-1, binSize), axis=1), trialBins

def sig_stars( p ):
    if np.isnan(p): return ''
    if p < 0.001:   return '***'
    if p < 0.01:    return '**'
    if p < 0.05:    return '*'
    return ''

##########################################################################
### Here we have functions for the panels for the figurel


def panelM_ThetaSchematic( ax ):
    print( "PANEL = M: Theta schematic" )
    duration = 0.52
    thetaFreq = 100/13.0
    t = np.linspace(0, duration, int(SAMPLE_FREQ * duration), 
        endpoint=False)
    y = np.sin(2 * np.pi * thetaFreq * t)
    TRIG = np.zeros_like( t )
    TRIG[int( round( 0.2*SAMPLE_FREQ ) )] = 1.0
    for ii in range(50):
        idx = int( round (ii*SAMPLE_FREQ/100) )
        if not(ii < 4 or ii in range(9,17) or ii in range(22,30) or ii in range( 35,43 ) or ii > 47):
            TRIG[idx] = 1.0
    ax.plot( t, y )
    ax.plot( t, TRIG * 0.5 - 1.7, color = "green" )
    ax.scatter( [0.13, 0.26, 0.39],[-2,-2,-2], marker = '^', color='red' )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim( -0.02, 0.7 )
    ax.set_ylim( -2, 1.3 )
    ax.set_yticklabels([])
    ax.set_xlabel( "Time (s)" )
    ax.text( -0.15, 1.06, "M", fontsize = 20, weight = "bold", transform=ax.transAxes )
    ax.text( 0.1, 0.99, "Schematic of burst stimulus and theta", fontsize = 14, transform=ax.transAxes )

def panelL_FreqSweep( ax ):
    print( "PANEL = L: spk freq dep" )
    dcell = pandas.read_hdf( "fig8L_spk_freq_sweep_orig_0.h5" )
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    spikeCriterion = -0.03
    freqs = [8, 20, 50, 80, 100, 150]
    transitions = [
        ( slice(5,8),  slice(8,11)  ),
        ( slice(13,16), slice(16,19) ),
        ( slice(21,24), slice(24,27) ),
    ]
    means = []
    sems = []
    pvals = []

    for ff in freqs:
        df = dcell.loc[ dcell['stimFreq'] == ff ]
        binSize = round( SAMPLE_FREQ / ff )
        trialBins = np.zeros( (len(df), 32) )
        for sweep in range( len(df) ):
            epsp = np.array( df.iloc[sweep, STIM_ON : STIM_ON + 32 * binSize] )
            spikeTrain = np.zeros_like( epsp, dtype=float )
            crossings = np.where( np.diff( epsp > spikeCriterion ) > 0 )
            spikeTrain[crossings] = 1
            trialBins[sweep] = spikeTrain.reshape( 32, binSize ).sum( axis=1 )

        ratios = []
        pre_total = 0
        post_total = 0
        for pre_sl, post_sl in transitions:
            pre  = trialBins[:, pre_sl].sum()
            post = trialBins[:, post_sl].sum()
            pre_total += pre
            post_total += post
            if pre > 0:
                ratios.append( post / pre )
        means.append( np.mean(ratios) if ratios else np.nan )
        sems.append( np.std(ratios) / np.sqrt(len(ratios)) if len(ratios) > 1 else np.nan )
        n_total = int( pre_total + post_total )
        if n_total > 0:
            pvals.append( scipy.stats.binomtest( int(post_total), n_total, p=0.5,
                alternative='greater' ).pvalue )
        else:
            pvals.append( np.nan )

    ax.errorbar( freqs, means, yerr=sems, fmt='o-', capsize=4 )
    ax.axhline( y=1.0, color='gray', linestyle='--', linewidth=1 )

    def sig_stars( p ):
        if np.isnan(p):  return ''
        if p < 0.001:    return '***'
        if p < 0.01:     return '**'
        if p < 0.05:     return '*'
        return 'ns'

    ylo, yhi = ax.get_ylim()
    star_y = yhi + 0.12 * (yhi - ylo)
    for x, p, m in zip( freqs, pvals, means ):
        label = sig_stars( p )
        if label and not np.isnan(m):
            ax.text( x, star_y, label, ha='center', va='top', fontsize=12 )
    print( "Panel L p-values: ",
        ["{}: {:.3g}".format(f, p) for f, p in zip(freqs, pvals)] )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Frequency (Hz)" )
    ax.set_ylabel( "post/pre" )
    ax.set_xlim( 0, 180 )
    ax.set_ylim( 0, 9 )
    ax.set_xticks( freqs )
    ax.set_xticklabels( [str(f) for f in freqs] )
    ax.text( -0.15, 1.06, "L", fontsize = 20, weight = "bold", transform=ax.transAxes )

def panelA_Vm( ax, fname ):
    dcell = pandas.read_hdf( fname )
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_OFF = int( 1.3 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_NUM_SAMPLES = STIM_OFF - STIM_ON
    sweep0 = 0
    sweep1 = 1
    sweep2 = 3

    time = np.array( np.arange( STIM_NUM_SAMPLES, dtype = float) ) / SAMPLE_FREQ
    epsp0 = np.array(df.iloc[sweep0, STIM_ON:STIM_OFF ])
    epsp1 = np.array(df.iloc[sweep1, STIM_ON:STIM_OFF ])
    epsp2 = np.array(df.iloc[sweep2, STIM_ON:STIM_OFF ])

    print( "Panel A LEN = ", len( time ), len( epsp0 ), len( epsp1 ) )
    ax.plot( time, epsp0, "b" )
    ax.plot( time, epsp1, "y" )
    ax.plot( time, epsp2, "g" )
    ax.scatter( [0.16, 0.32, 0.48],[-75,-75,-75], marker = '^', color='red' )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "Vm (mV)" )
    ax.set_xticklabels([])
    #ax.legend( loc = "upper right", frameon = False, fontsize = 14 )
    ax.set_xlim( -0.02, 0.82 )
    ax.set_ylim( -80, 20 )
    ax.text( -0.07, 1.05, "A", fontsize = 20, weight = "bold", transform=ax.transAxes )

def panelB_raster( ax, fname ):
    spikeCriterion = -0.03
    dcell = pandas.read_hdf( fname )
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_OFF = int( 1.3 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_NUM_SAMPLES = STIM_OFF - STIM_ON
    time = np.array( np.arange( STIM_NUM_SAMPLES, dtype = float) ) / SAMPLE_FREQ
    for sweep in range( len(df) ):
        epsp = np.array(df.iloc[sweep, STIM_ON:STIM_OFF ])
        t = np.where(np.diff(epsp > spikeCriterion) > 0)[0]
        t = t/SAMPLE_FREQ
        y = np.ones_like( t, dtype=float ) * sweep
        ax.scatter( t, y, color="blue", marker = '.' )

    print( "Panel B: ", len(t), len(y) )
    ax.scatter( [0.16, 0.32, 0.48],[-5,-5,-5], marker = '^', color='red' )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "Trial #\n" )
    #ax.set_xlabel( "Time (s)" )
    ax.set_xticklabels([])
    ax.set_xlim( -0.02, 0.82 )
    ax.text( -0.07, 1.05, "B", fontsize = 20, weight = "bold", transform=ax.transAxes )

def printStats( title, panel, bins ):
    # Each bin is 20 ms. Binomial test: under H0 (ratio=1) each spike is equally
    # likely to be pre or post (windows are equal width), so p_null=0.5.
    def binom_p( pre_sl, post_sl ):
        pre  = int( bins[pre_sl].sum() )
        post = int( bins[post_sl].sum() )
        n = pre + post
        if n == 0:
            return np.nan
        return scipy.stats.binomtest( post, n, p=0.5, alternative='greater' ).pvalue

    w1 = binom_p( slice(5,8),  slice(8,11)  )
    w2 = binom_p( slice(13,16), slice(16,19) )
    w3 = binom_p( slice(21,24), slice(24,27) )
    denom = sum(bins[5:8])
    r1 = sum(bins[8:11])/denom if denom > 0 else -1
    denom = sum(bins[13:16])
    r2 = sum(bins[16:19])/denom if denom > 0 else -1
    denom = sum(bins[21:24])
    r3 = sum(bins[24:27])/denom if denom > 0 else -1
    print( "p{}: w1={:12.4g}, w2={:12.4g}, w3={:12.4g} {}".format( panel, w1, w2, w3, title ) )
    print( "{} : r1={:12.4g}, r2={:12.4g}, r3={:12.4g} {}".format( panel, r1, r2, r3, title ) )

    # Kendall's tau: one-sided test for monotonic decline over 32 pulses.
    tau, p_decline = scipy.stats.kendalltau( np.arange(32), bins[:32], alternative='less' )
    print( "{} : tau={:8.4f}, p_decline={:12.4g} {}".format( panel, tau, p_decline, title ) )
    return w1, w2, w3


def panelCK_SampleTrace( ax, dcell, panel, title ):
    print( "PANEL = ", panel )
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    time, rate, binned, trialBins = spikeRate( df, spikeCriterion=-0.03, windowSize = 0.01 )
    ax.plot( time, rate, "b" )
    if panel not in ['L','M']:
        ax.scatter( [0.16, 0.32, 0.48],[-5,-5,-5], marker='^', color='red' )
        w1, w2, w3 = printStats( title, panel, binned )
        trans = blended_transform_factory( ax.transData, ax.transAxes )
        for x, p in zip( [0.16, 0.32, 0.48], [w1, w2, w3] ):
            stars = sig_stars( p )
            if stars:
                ax.text( x, 0.90, stars, transform=trans, ha='center', va='top', fontsize=14 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if panel == 'C':
        ax.text( -0.07, 0.90, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
    elif panel in ['D', 'E']:
        ax.text( -0.15, 0.99, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
        ax.set_xticklabels([])
    else:
        ax.text( -0.15, 1.05, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
        if not panel in ['J', 'K']:
            ax.set_xticklabels([])
    if panel in ['H', 'D']:
        ax.set_ylabel( "Mean Spike Rate (Hz)\n" )
    if panel in ['C', 'J','K']:
        ax.set_xlabel( "Time (s)" )
    ax.set_xlim( -0.02, 0.82 )
    ax.text( 0.1, 0.99, title, fontsize = 14, transform=ax.transAxes )

def printThetaBurstStats( title, panel, rate ):
    # Binomial test comparing burst 1 vs bursts 2, 3, 4 (two-sided: rate could go
    # up or down). Each window is 4 x 10ms bins; equal widths so p_null=0.5.
    # Burst positions (100Hz bin indices): burst1: 6:10, burst2: 19:23, etc.
    def binom_p( ref_sl, cmp_sl ):
        ref = int( rate[ref_sl].sum() )
        cmp = int( rate[cmp_sl].sum() )
        n = ref + cmp
        if n == 0:
            return np.nan
        return scipy.stats.binomtest( cmp, n, p=0.5, alternative='two-sided' ).pvalue

    w1 = binom_p( slice(6,10), slice(19,23) )
    w2 = binom_p( slice(6,10), slice(32,36) )
    w3 = binom_p( slice(6,10), slice(45,49) )

    r1 = np.mean(rate[19:23])/ np.mean(rate[6:10])
    r2 = np.mean(rate[32:36])/ np.mean(rate[6:10])
    r3 = np.mean(rate[45:49])/ np.mean(rate[6:10])
    print( "{} : w1={:12.4g}, w2={:12.4g}, w3={:12.4g} {}".format( panel, w1, w2, w3, title ) )
    print( "{} : r1={:12.4g}, r2={:12.4g}, r3={:12.4g} {}".format( panel, r1, r2, r3, title ) )
    return w1, w2, w3

def panelPQRS_ThetaSampleTrace( ax, dcell, panel, title ):
    print( "PANEL = ", panel )
    df = dcell.loc[(dcell['stimFreq'] == 100)]
    time, rate, binned, trialBins = rippleSpikeRate( df, spikeCriterion=-0.03, windowSize = 0.005 )
    w1, w2, w3 = printThetaBurstStats( title, panel, binned )
    ax.plot( time, rate, "b" )
    trans = blended_transform_factory( ax.transData, ax.transAxes )
    for x, p in zip( [0.13, 0.26, 0.39], [w1, w2, w3] ):
        stars = sig_stars( p )
        if stars:
            ax.text( x, 0.80, stars, transform=trans, ha='center', va='top', fontsize=14 )
    ax.scatter( [0.13, 0.26, 0.39],[-10,-10,-10], marker = '^', color='red')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    if panel == 'N':
        ax.set_ylabel ( "Rate (Hz)" )
    ax.set_xlim( -0.02, 0.54 )
    ax.set_ylim( -20, 180 )
    ax.text( -0.15, 1.05, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
    ax.text( 0.1, 0.89, title, fontsize = 14, transform=ax.transAxes )


def main():
    path = "SIMDATA"
    fnames = [  f"{path}/fig8_reference_orig_0.h5",
                f"{path}/fig8_zeroIndices_32.h5",
                f"{path}/fig8_zeroIndices_240.h5",
                f"{path}/fig8_oddball_orig_0.h5",
                f"{path}/fig8_gap_orig_0.h5",
                f"{path}/fig8_modelName_BothPresyn90_noGlu_STP.g.h5",
                f"{path}/fig8_modelName_BothPresyn90_noGABA_STP.g.h5",
                f"{path}/fig8_uniform_orig_0.h5",
                f"{path}/fig8_random_orig_0.h5",
                f"{path}/fig8_theta_orig_0.h5",
                f"{path}/fig8_theta_uniform_orig_0.h5",
    ]

    titles = [  
                "Reference",
                "Dense stimulus",
                "Sparse stimulus",
                "Oddball",
                "Gap",
                "No STP in Glu",
                "No STP in GABA",
                "Uniform pattern",
                "Random pattern",
                "Theta burst mismatch",
                "Theta burst uniform",
    ]

    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (12,21), layout='constrained' )
    gs_outer = fig.add_gridspec( 3, 1, height_ratios=[3, 4, 2] )

    # Section 1: Panels A, B, C — full-width, A and B suppress x labels
    gs1 = gs_outer[0].subgridspec( 3, 1, hspace=0.41 )
    ax = fig.add_subplot( gs1[0, 0] )
    panelA_Vm( ax, fnames[0] )
    ax = fig.add_subplot( gs1[1, 0] )
    panelB_raster( ax, fnames[0] )
    ax = fig.add_subplot( gs1[2, 0] )
    df = pandas.read_hdf( fnames[0] )
    panelCK_SampleTrace( ax, df, "C", titles[0] )

    # Section 2: Panels D–K — 4 rows × 2 cols, D–I suppress x labels
    gs2 = gs_outer[1].subgridspec( 4, 2, hspace=0.35 )
    for idx, fname in enumerate( fnames[1:9] ):
        print( "loading: ", fname )
        panel = chr( ord("D") + idx )
        df = pandas.read_hdf( fname )
        ax = fig.add_subplot( gs2[idx//2, idx%2] )
        panelCK_SampleTrace( ax, df, panel, titles[idx+1] )

    # Section 3: Panels L, M, N, O — 2 rows × 2 cols, both rows carry xlabels
    gs3 = gs_outer[2].subgridspec( 2, 2, hspace=0.85 )
    ax = fig.add_subplot( gs3[0, 0] )
    panelL_FreqSweep( ax )
    ax = fig.add_subplot( gs3[0, 1] )
    panelM_ThetaSchematic( ax )
    for idx, fname in enumerate( fnames[9:] ):
        print( "loading: ", fname )
        panel = chr( ord("N") + idx )
        df = pandas.read_hdf( fname )
        ax = fig.add_subplot( gs3[1, idx] )
        panelPQRS_ThetaSampleTrace( ax, df, panel, titles[idx+9] )

    # Panel letters sit at y=1.05 in axes coords, outside the axes bounding box.
    # constrained_layout accounts for them and adds whitespace at the figure top.
    # Excluding them from layout computation removes that whitespace; they still render.
    for ax in fig.get_axes():
        for txt in ax.texts:
            if txt.get_transform() is ax.transAxes:
                _, y = txt.get_position()
                if y > 1.0:
                    txt.set_in_layout( False )

    # x range 0 to 1. Tweaked ylabel_x so labels sit aligned.
    ylabel_x = [0.045, 0.045, 0.075, 0.075, 0.075, 0.045 ]
    for ax in fig.get_axes():
        if ax.get_ylabel():
            ax.yaxis.set_label_coords(
                ylabel_x.pop(), 0.5,
                transform=blended_transform_factory( fig.transFigure, ax.transAxes )
            )
    plt.show()

if __name__ == "__main__":
    main()


