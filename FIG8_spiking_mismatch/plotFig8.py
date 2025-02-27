import matplotlib.pyplot as plt
import pandas
#import pylab
import numpy as np
import math
import argparse
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
    for sweep in range( len(df) ):
        epsp = np.array(df.iloc[sweep, STIM_ON:STIM_OFF ])
        spikeTrain = np.zeros_like(epsp, dtype=float)
        crossings = np.where(np.diff(epsp > spikeCriterion) > 0)
        spikeTrain[crossings] = 1
        sumTrain += spikeTrain

    # Convolve the summed spike trains with a rectangular window
    window = np.ones(int(windowSize * SAMPLE_FREQ ))
    firingRate = np.convolve(sumTrain, window, mode='same') / (windowSize * len( df ) )

    return time, firingRate

def spikeRate( dcell, spikeCriterion = -0.03, windowSize = 0.02 ):
    STIM_ON = int( 0.5 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_OFF = int( 1.3 * SAMPLE_FREQ ) + SAMPLE_START
    STIM_NUM_SAMPLES = STIM_OFF - STIM_ON
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    time = np.array( np.arange( STIM_NUM_SAMPLES, dtype = float) ) / SAMPLE_FREQ
    sumTrain = np.zeros( STIM_NUM_SAMPLES )
    for sweep in range( len(df) ):
        epsp = np.array(df.iloc[sweep, STIM_ON:STIM_OFF ])
        spikeTrain = np.zeros_like(epsp, dtype=float)
        crossings = np.where(np.diff(epsp > spikeCriterion) > 0)
        spikeTrain[crossings] = 1
        sumTrain += spikeTrain

    # Convolve the summed spike trains with a rectangular window
    window = np.ones(int(windowSize * SAMPLE_FREQ ))
    firingRate = np.convolve(sumTrain, window, mode='same') / (windowSize * len( df ) )

    # Optionally, Downsample to match the desired step size
    # time_points = time[::int(step_size / dt)]
    # firing_rates = firing_rate[::int(step_size / dt)]
    # return time_points, firing_rates
    return time, firingRate

##########################################################################
### Here we have functions for the panels for the figurel


def panelO_ThetaSchematic( ax ):
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
    ax.set_xlim( -0.02, 0.54 )
    ax.set_yticklabels([])
    ax.set_xlabel( "Time (s)" )
    ax.text( -0.15, 1.05, "O", fontsize = 20, weight = "bold", transform=ax.transAxes )
    ax.text( 0.1, 1.05, "Schematic of burst stimulus and theta", fontsize = 14, transform=ax.transAxes )

def panelN_Freq( ax ):
    df = pandas.read_hdf( "SURP39_Freq_Sweep/surprise_param_dep2.h5" )
    params = df['param'].unique()
    print( "Panel N freqs" )
    for idx, pp in enumerate( params ):
        dp = df.loc[df["param"]==pp]
        flist = dp["freq"].unique()
        maxy = 0.0
        freqs = np.array( flist ) 
        ay1 = []
        ay2 = []
        ay3 = []
        for ii, ff in enumerate( reversed(flist) ):
            dfreq = dp.loc[dp["freq"]==ff]
            y1 = np.array(dfreq["m2"])/np.array(dfreq["m1"]).flatten()
            y2 = np.array(dfreq["m4"])/np.array(dfreq["m3"]).flatten()
            y3 = np.array(dfreq["m6"])/np.array(dfreq["m5"]).flatten()
            ay1.append( np.mean( y1 ) )
            ay2.append( np.mean( y2 ) )
            ay3.append( np.mean( y3 ) )
            maxy = max( [maxy, max( ay1 ), max( ay2 ), max( ay3 )] )
    y = np.array( [np.array( ay1 ), np.array(ay2), np.array(ay3) ] ) * 100
    ym = np.mean( y, axis = 0 )
    ys = np.std( y, axis = 0 )
    ax.errorbar( freqs, ym, yerr = ys, fmt='o-' )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Frequency (Hz)" )
    ax.set_ylabel( "% base" )
    ax.text( -0.15, 1.05, "N", fontsize = 20, weight = "bold", transform=ax.transAxes )

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

def panelCK_SampleTrace( ax, dcell, panel, title ):
    print( "PANEL = ", panel )
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    time, rate = spikeRate( df, spikeCriterion=-0.03, windowSize = 0.01 )
    ax.plot( time, rate, "b" )
    if panel not in ['L','M']:
        ax.scatter( [0.16, 0.32, 0.48],[-5,-5,-5], marker = '^', color='red' )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if panel == 'C':
        ax.text( -0.07, 1.05, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
    else:
        ax.text( -0.15, 1.05, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
        if not panel in ['L', 'M']:
            ax.set_xticklabels([])
    if panel in ['J', 'D']:
        ax.set_ylabel( "Mean Spike Rate (Hz)\n" )
    if panel in ['L','M']:
        ax.set_xlabel( "Time (s)" )
    ax.set_xlim( -0.02, 0.82 )
    ax.text( 0.1, 1.05, title, fontsize = 14, transform=ax.transAxes )

def panelPQRS_ThetaSampleTrace( ax, dcell, panel, title ):
    print( "PANEL = ", panel )
    df = dcell.loc[(dcell['stimFreq'] == 100)]
    time, rate = rippleSpikeRate( df, spikeCriterion=-0.03, windowSize = 0.005 )
    ax.plot( time, rate, "b" )
    ax.scatter( [0.13, 0.26, 0.39],[-10,-10,-10], marker = '^', color='red')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    if panel == 'N':
        ax.set_ylabel ( "Rate (Hz)" )
    ax.set_xlim( -0.02, 0.54 )
    ax.set_ylim( -20, 180 )
    ax.text( -0.15, 1.05, panel, fontsize = 20, weight = "bold", transform=ax.transAxes )
    ax.text( 0.1, 1.05, title, fontsize = 14, transform=ax.transAxes )


def main():
    fnames = [  "SURP38_fixed/simData_zeroIndices_192.h5",
                "SURP38_fixed/simData_zeroIndices_64.h5",
                "SURP38_fixed/simData_zeroIndices_224.h5",
                "SURP40_deviant/fixwt_dev_orig_0.h5",
                "SURP41_GAP/fixwt_gap_orig_0.h5",
                "SURP48_noSTP_fixed/simData_modelName_BothPresyn90_noGlu_STP.g.h5",
                "SURP48_noSTP_fixed/simData_modelName_BothPresyn90_noGABA_STP.g.h5",
                "SURP48_noSTP_fixed/simData_modelName_BothPresyn90_noSTP.g.h5",
                "SURP45_noGABA/simData_wtGlu5.h5",
                "SURP43_uniform_pattern/fixwt_unif_orig_0.h5",
                "SURP44_rand_pattern/fixwt_rand_orig_0.h5",
                "SURP46_theta_pat/fixwt_theta_pat_orig_0.h5",
                "SURP47_theta_uniform/fixwt_theta_unif_orig_0.h5",
                "SURP50_theta_gamma_precess/prec_precession_13.h5",
                "SURP50_theta_gamma_precess/prec_precession_12.h5"
    ]

    titles = [  
                "Reference",
                "Dense connectivity",
                "Sparse connectivity",
                "Oddball",
                "Gap",
                "No STP in Glu",
                "No STP in GABA",
                "No STP in either",
                "No GABA",
                "Uniform pattern",
                "Random pattern",
                "Theta burst mismatch",
                "Theta burst uniform",
                "Theta burst no precess",
                "Theta burst precess"
    ]

    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (12,17) )
    gs = fig.add_gridspec( 8, 2 ) # 8 rows, 2 cols
    ax = fig.add_subplot( gs[0,:] ) # Do a timeseries voltage plot
    panelA_Vm( ax, fnames[0] )
    ax = fig.add_subplot( gs[1,:] ) # Do a raster plot.
    panelB_raster( ax, fnames[0] )
    ax = fig.add_subplot( gs[2,:] ) # Do reference spike rate plot
    df = pandas.read_hdf( fnames[0] )
    panelCK_SampleTrace( ax, df, "C", titles[0] )
    for idx, fname in enumerate( fnames[1:11] ):
        print( "loading: ", fname )
        panel = chr( ord("D") + idx )
        df = pandas.read_hdf( fname )
        ax = fig.add_subplot( gs[3+idx//2,idx%2] )
        panelCK_SampleTrace( ax, df, panel, titles[idx+1] )

    fig2 = plt.figure( figsize = (11.65, 8) )
    gs = fig2.add_gridspec( 3, 2 ) # 3 rows, 2 cols
    ax = fig2.add_subplot( gs[0,0] ) # Do freq plot
    panelN_Freq( ax )
    ax = fig2.add_subplot( gs[0,1] ) # Do schematic
    panelO_ThetaSchematic( ax )
    for idx, fname in enumerate( fnames[11:] ):
        print( "loading: ", fname )
        panel = chr( ord("P") + idx )
        df = pandas.read_hdf( fname )
        ax = fig2.add_subplot( gs[1 + idx//2,idx%2] )
        panelPQRS_ThetaSampleTrace( ax, df, panel, titles[idx+11] )


    fig.tight_layout()
    fig2.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()


