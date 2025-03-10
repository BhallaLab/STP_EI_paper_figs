import pandas
import pylab
import numpy as np
import math
import argparse
from scipy.stats import linregress
import scipy.stats as stats

import matplotlib.pyplot as plt

#freq = 80.0 # Hz
#settleTime = 0.1    # seconds
#stimDuration = 0.002   # seconds
#postStim = 0.4
#stimAmpl = 5e-2     # mM
#basalCa = 0.08e-3   # mM
#GABAdelay = 5.0e-3  # seconds
#width = 0.002

pulseTrig = []
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 1
NUM_SAMPLES = int(SAMPLE_FREQ * SAMPLE_TIME)
SAMPLE_START = 49
PulseTrain = {}

def panelA_SampleTrace( ax, dcell, dsim, freq, label ):
    df = dcell.loc[(dcell['stimFreq'] == freq) & (dcell['numSq'] == 5)]
    cellid = df['cellID'].unique()[0]
    shift = int( SAMPLE_FREQ * 0.3 ) - 600 if cellid in [5611, 1941] else 0
    print( cellid, shift )
    PLOTLEN = 2.5
    epsp = np.mean(np.array(df.iloc[:, SAMPLE_START:SAMPLE_START+NUM_SAMPLES+shift]), axis = 0 )
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    df2 = dsim.loc[dsim['stimFreq'] == freq]
    simepsp = np.mean(np.array(df2.iloc[:, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ]), axis = 0 )
    #print( "SHAPES = ", df.shape, df2.shape, epsp.shape, simepsp.shape )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    pulses = np.arange(0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]

    baseline = np.median( epsp[:int(0.2*SAMPLE_FREQ)] )
    epsp -= baseline # hack to handle traces with large ipsps.
    baseline = np.median( simepsp[:int(0.2*SAMPLE_FREQ)] )
    simepsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, SAMPLE_TIME, SAMPLE_FREQ )
    #tepsp = np.linspace( 0, SAMPLE_TIME, len(epsp) )
    #tepsp = tepsp[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsp)]
    ax.plot( tepsp, epsp[shift:len(tepsp)+shift], "m", label = "Expt" )
    ax.plot( tepsp, simepsp[:len(tepsp)], "y", label = "Sim" )
    ax.plot( tepsp, pt * 0.5 - 1, "g", label = "Trig" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "EPSP (mV)" )
    ax.legend( loc = "upper right", frameon = False, fontsize = 14 )
    ax.set_ylim( -1.2, max(2, max(epsp)*1.05 ) )
    ax.text( 0.05, 0.90, str(freq)+" Hz", fontsize = 16, transform=ax.transAxes )
    ax.text( -0.28, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    if freq == 50:
        ax.set_xlabel("Time (s)")
    return( max( simepsp ) )

def computeEPSPmax( dcell, freq ):
    df = dcell.loc[(dcell['stimFreq'] == freq) & (dcell['numSq'] == 5)]
    epsp = np.mean(np.array(df.iloc[:, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ]), axis = 0 )
    baseline = np.median( epsp[:int(0.2*SAMPLE_FREQ)] )
    epsp -= baseline # hack to handle traces with large ipsps.
    return( max( epsp ) )

def panelEPSPmax( ax, df, simEPSPmax, label ):
    cellEmax = []
    freqs = [20, 30, 40, 50]
    for cc in [1941, 5611, 5291, 5212, 5211, 3791, 3531, 3402, 3373]:
        dcell = df.loc[df['cellID']==cc]
        emax = np.array( [ computeEPSPmax(dcell, ff) for ff in freqs ] )
        emax /= np.mean( emax )
        cellEmax.append( emax )
        ax.scatter( freqs, emax, alpha = 0.6 )
    ne = np.array( cellEmax )
    print( "SHAPE NE = ",ne.shape )
    means = np.mean( ne, axis = 0 )
    std = np.std( ne, axis = 0 )
    ax.errorbar( freqs, means, yerr=std, fmt='o', color="maroon", capsize=5 )
    ax.plot( freqs, means, marker='o', color="m", linestyle='-', 
        linewidth=2, label = "Expt")
    simEPSPmax /= np.mean( simEPSPmax )
    ax.plot( freqs, simEPSPmax, marker='+', color="y", linestyle='-', 
        linewidth=2, label = "Sim")
    ax.set_xlabel( "Burst Frequency (Hz)" )
    ax.set_ylabel('Norm. EPSP')
    ax.set_ylim( 0.4, 1.7 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text( -0.28, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    ax.legend( loc = "upper center", frameon = False, fontsize=14, ncol=2 )
    ft = np.tile( freqs, ne.shape[0])
    slope, intercept, r_value, p_value, std_err = linregress(ft, ne.flatten())

    print("Linear regression results:")
    print(f"Slope: {slope}")
    print(f"Intercept: {intercept}")
    print(f"R-squared: {r_value**2}")
    print(f"P-value: {p_value}")

    spearman_corr, spearman_p_value = stats.spearmanr(ft, ne.flatten())
    print("\nSpearman's rank correlation results:")
    print(f"Spearman correlation: {spearman_corr}")
    print(f"P-value: {spearman_p_value}")

    kendall_tau, kendall_p_value = stats.kendalltau(ft, ne.flatten())
    print("\nKendall's tau correlation results:")
    print(f"Kendall's tau: {kendall_tau}")
    print(f"P-value: {kendall_p_value}")

def panelEI( ax, dcell, freq, label ):
    df = dcell.loc[(dcell['stimFreq'] == freq)]
    PLOTLEN = 2.5
    # EPSC data is in and IPSC data is in block 1 and 3 resp.
    epsc = np.array(df.iloc[0, SAMPLE_START+NUM_SAMPLES*1:SAMPLE_START+NUM_SAMPLES*2 ])
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    ipsc = np.array(df.iloc[0, SAMPLE_START+NUM_SAMPLES*3:SAMPLE_START+NUM_SAMPLES*4 ])
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    pulses = np.arange(0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    #print( PulseTrain )

    baseline = np.median( epsc[:int(0.2*SAMPLE_FREQ)] )
    '''
    if dcell['cellID'][0] == 0:
        baseline = min(epsp )
    else:
        baseline = min( np.percentile(epsp, 25 ), 0.0 )
    '''
    epsc -= baseline # hack to handle traces with large ipsps.
    tepsc = np.linspace( 0, SAMPLE_TIME, len(epsc) )

    tepsc = tepsc[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsc)]
    ax.plot( tepsc, -epsc[:len(tepsc)], "b", label = "epsc " )
    ax.plot( tepsc, -ipsc[:len(tepsc)], "r", label = "ipsc " )
    low = max( epsc )
    ax.plot( tepsc, pt * low/10 - low*1.2, "g", label = "Trig" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "E/IPSC (pA)" )
    ax.legend( loc = "lower right", frameon = False, fontsize = 14 )
    #ax.set_ylim( -1.2, max(2, max(epsc)+1 ) )
    ax.text( 0.05, 0.90, str(freq)+" Hz", fontsize = 16, transform=ax.transAxes )
    ax.text( -0.28, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return -max( epsc ), -min( ipsc ) 

def panelCond( ax, dcell, freq, label ):
    numExc = 100
    numInh = 200
    df = dcell.loc[(dcell['stimFreq'] == freq)]
    PLOTLEN = 2.5
    # EPSC data is in and IPSC data is in block 1 and 3 resp.
    epsg = np.array(df.iloc[0, SAMPLE_START+NUM_SAMPLES*4:SAMPLE_START+NUM_SAMPLES*5 ])
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    ipsg = np.array(df.iloc[0, SAMPLE_START+NUM_SAMPLES*5:SAMPLE_START+NUM_SAMPLES*6 ])
    #field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    pulses = np.arange(0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    #print( PulseTrain )

    #baseline = np.median( epsc[:int(0.2*SAMPLE_FREQ)] )
    #epsc -= baseline # hack to handle traces with large ipsps.
    tepsg = np.linspace( 0, SAMPLE_TIME, len(epsg) )

    tepsg = tepsg[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsg)]
    ax.plot( tepsg, 1e-3*epsg[:len(tepsg)], "b", label = "epsG " )
    ax.plot( tepsg, 1e-3*ipsg[:len(tepsg)], "r", label = "ipsG " )
    low = max( epsg*1e-3 )
    ax.plot( tepsg, pt * low/10 - low*1.2, "g", label = "Trig" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "E/IPSG (nS)" )
    ax.legend( loc = "upper right", frameon = False, fontsize = 14 )
    #ax.set_ylim( -1.2, max(2, max(epsc)+1 ) )
    ax.text( 0.05, 0.90, str(freq)+" Hz", fontsize = 16, transform=ax.transAxes )
    ax.text( -0.28, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return max( epsg )/1000, max( ipsg )/1000 

def panelJ( ax, ie, ii, f ):
    ax.plot( f, -np.array(ie), "b", label = "epsc" )
    ax.plot( f, ii, "r", label = "ipsc" )
    ax.plot( f, -np.array(ie) / np.array(ii), "k", label = "ratio" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "E/IPSC pk (pA)" )
    ax.set_xlabel( "Burst Freq (Hz)" )
    ax.set_ylim( 1, 6.6 )
    ax.legend( loc = "upper left", frameon = False, fontsize = 14 )
    ax.text( -0.28, 1.05, "J", fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelO( ax, ie, ii, f ):
    ax.plot( f, np.array(ie), "b", label = "epsG" )
    ax.plot( f, ii, "r", label = "ipsG" )
    ax.plot( f, np.array(ie) / np.array(ii), "k", label = "ratio" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "Gmax (nS)" )
    ax.set_xlabel( "Burst Freq (Hz)" )
    ax.set_ylim( 0, 0.8 )
    ax.legend( loc = "upper left", frameon = False, fontsize = 14 )
    ax.text( -0.28, 1.05, "O", fontsize = 22, weight = "bold", transform=ax.transAxes )



def main():
    global pulseTrig
    parser = argparse.ArgumentParser( description = "Read and display surprise protocol data from specified file" )
    parser.add_argument( "-e", "--exptFname", type = str, help = "Optional: File with exptdata. Required to be hdf5.", default = "../../../2022/VC_DATA/all_cells_FreqSweep_CC_long.h5"  )
    parser.add_argument( "-s", "--simFname", type = str, help = "Optional: File with simdata. Required to be hdf5.", default = "simData_orig_0.h5" )
    parser.add_argument( "-c", "--cell", type = int, help = "Optional: Cell number to use. default = 3531", default = 3531 )
    args = parser.parse_args()
    cellid = 1941
    cellid = 5611
    cellid = 5291
    cellid = 5212
    cellid = 5211
    cellid = 3791
    cellid = 3531
    cellid = 3402
    cellid = 3373
    #cellid = 3781 This one is too small epsp.

    cellid = 3872 # spikes
    cellid = 3871 # spikes
    cellid = 3531
    cellid = args.cell


    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (14,18) )
    fig.suptitle( "Cell = {}".format( cellid ), fontsize = 16 )
    gs = fig.add_gridspec( 5, 3 ) # 5 rows, 3 cols
    #fig, ax = plt.subplots( nrows = 3, ncols=3, figsize = (18, 15) )

    # Set up the stimulus timings
    exptdf = pandas.read_hdf( args.exptFname )
    #cells = df['cellID'].unique()
    celldf = exptdf.loc[exptdf['cellID']== cellid]
    simdf = pandas.read_hdf( args.simFname )

    vmax = [0.0]*4
    iemax = [0.0]*4
    iimax = [0.0]*4
    gemax = [0.0]*4
    gimax = [0.0]*4
    vmax[0] = panelA_SampleTrace( fig.add_subplot( gs[0,0]), celldf, simdf, 20, "A" )
    vmax[1] = panelA_SampleTrace( fig.add_subplot( gs[1,0]), celldf, simdf, 30, "B" )
    vmax[2] = panelA_SampleTrace( fig.add_subplot( gs[2,0]), celldf, simdf, 40, "C" )
    vmax[3] = panelA_SampleTrace( fig.add_subplot( gs[3,0]), celldf, simdf, 50, "D" )
    #panelI( fig.add_subplot( gs[4,0] ), vmax, [20, 30, 40, 50] )
    panelEPSPmax( fig.add_subplot( gs[4,0] ), exptdf, vmax, "E" )
    iemax[0], iimax[0] = panelEI( fig.add_subplot( gs[0,1]), simdf, 20, "F")
    iemax[1], iimax[1] = panelEI( fig.add_subplot( gs[1,1]), simdf, 30, "G")
    iemax[2], iimax[2] = panelEI( fig.add_subplot( gs[2,1]), simdf, 40, "H")
    iemax[3], iimax[3] = panelEI( fig.add_subplot( gs[3,1]), simdf, 50, "I")
    panelJ( fig.add_subplot( gs[4,1] ), iemax, iimax, [20, 30, 40, 50] )

    gemax[0],gimax[0] = panelCond(fig.add_subplot( gs[0,2]), simdf, 20, "K")
    gemax[1],gimax[1] = panelCond(fig.add_subplot( gs[1,2]), simdf, 30, "L")
    gemax[2],gimax[2] = panelCond(fig.add_subplot( gs[2,2]), simdf, 40, "M")
    gemax[3],gimax[3] = panelCond(fig.add_subplot( gs[3,2]), simdf, 50, "N")
    panelO( fig.add_subplot( gs[4,2] ), gemax, gimax, [20, 30, 40, 50] )
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
