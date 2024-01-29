import pandas
import pylab
import numpy as np
import math
import matplotlib.pyplot as plt

freq = 80.0 # Hz
settleTime = 0.1    # seconds
stimDuration = 0.002   # seconds
postStim = 0.4
stimAmpl = 5e-2     # mM
basalCa = 0.08e-3   # mM
GABAdelay = 5.0e-3  # seconds
width = 0.002
doPlot = False

gluStimStr = "8e-5"
GABAStimStr = "8e-5"
gluR_clamp_potl = "-0.07"
GABAR_clamp_potl = "0.0"
GABAR_clamp_offset = 0.1    # nA
gluConductanceScale = 0.5   # Relative to default value in the spine proto
gluTau2Scale = 4   # Relative to default value in the spine proto

numCA1Exc = 100
numCA1Inh = 200
pCA3_CA1 = 0.0002
pCA3_Inter = 0.0008
pInter_CA1 = 1.0/256.0
interState = 0
thresh_CA3_Inter = 0.9999   # Avoid doing exact float comparisons to cross thresh.
thresh_CA3_CA1 = 0.9999
thresh_Inter_CA1 = 0.9999
repeatPatterns = False
inputs = []
stimList = []
pulseTrig = []
patternData = "../../../2022/VC_DATA/all_cells_SpikeTrain_CC_long.h5"
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 11
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
PKDELAY = int( 0.020 * SAMPLE_FREQ )
PKTHRESH = 0.0005 # Threshold for a distinct EPSP pk. In Volts
ALPHATAU1 = 0.005 * SAMPLE_FREQ
ALPHATAU2 = 0.042 * SAMPLE_FREQ
ALPHADELAY = int( 0.010 * SAMPLE_FREQ )
ALPHAWINDOW = int( ALPHATAU2 * 8 )
MINIMUM_PULSE_THRESH = 4e-4

def alphaFunc( t, tp ):
    return (t/tp) * np.exp(1-t/tp)

def dualAlphaFunc( t, t1, t2 ):
    if t < 0:
        return 0.0
    if abs( t1 - t2 ) < 1e-6:
        return alphaFunc( t, t1 )
    return (1.0/(t1-t2)) * (np.exp(-t/t1) - np.exp(-t/t2))

ALPHA = np.array( [ dualAlphaFunc(t, ALPHATAU1, ALPHATAU2 ) for t in range( ALPHAWINDOW ) ] )
ALPHA = ALPHA / max( ALPHA )
longAlpha = np.zeros(NUM_SAMPLES)
longAlpha[:ALPHAWINDOW] += ALPHA

def findPeaks( pulseTrig, pulseThresh, pkDelay, Vm, width = 0.002 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    lastPP = 0.0
    pks = []
    pkIdx = []
    for idx, pp in enumerate( pulseTrig ):
        if pp > pulseThresh and lastPP < pulseThresh:
            idx2 = idx + pkDelay - widthSamples
            pkIdx.append( idx )
            pks.append(np.median( Vm[idx2:idx2 + widthSamples*2] ))
        lastPP = pp

    return pkIdx, pks

def doStats( df, alphaTab, alphaDelay, pkDelay ):
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    print( "baseline={:.3f}, pulseThresh = {:.6f}".format( 
        baseline, pulseThresh ) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, SAMPLE_TIME, len(epsp) )
    pkIdx, pks = findPeaks( pulseTrig, pulseThresh, pkDelay, epsp )
    fitEPSP = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    foundPks = []
    foundIdx = []
    for idx, pp in zip( pkIdx, pks ):
        ii = idx + alphaDelay
        ascale = pp - lastPP*longAlpha[int(ALPHATAU1*2) + idx-lastIdx]
        if ascale > 0:
            fitEPSP[ii:ii+ALPHAWINDOW] += ascale * alphaTab
            lastIdx = idx
            lastPP = pp
            foundPks.append( ascale )
            foundIdx.append( idx )

    runtime = settleTime + SAMPLE_TIME + postStim

    if doPlot:
        plt.rcParams.update( {"font.size": 24} )
        fig = plt.figure( "cell = ", figsize = (32, 3.5) )


        tpt = np.arange( 0.0, len( pulseTrig ), 1 )/SAMPLE_FREQ
        plt.plot( tpt, pulseTrig * 100, "g", label = "Trigger" )
        plt.ylabel( "Vm (mV)" )
        plt.xlim( settleTime - 0.1, runtime + 0.1 )
        plt.plot( tepsp, epsp, "b", label = "Data EPSP" )
        plt.plot( tepsp, fitEPSP[:len(tepsp)], "r", label = "Fit EPSP" )
        plt.ylabel( "EPSP (mV)" )
        plt.legend()
        plt.xlim( settleTime - 0.1, runtime + 0.1 )
        plt.ylim( -5, 10 )
        plt.xlabel("Time (s)")
    
        fig.tight_layout()
        plt.show()

    return foundIdx, foundPks, len( pkIdx )

def analyze( ax, cell, num5, pkVal5, num15, pkVal15 ):
    prob5 = len( pkVal5 ) / num5
    prob15 = len( pkVal15 ) / num15
    print( "prob5={:.3f},{:.3f}, prob15={:.3f},{:.3f}".format( 
        prob5, prob5 * sum( np.array(pkVal5) > 0.5 )/len( pkVal5 ),
        prob15, prob15 * sum( np.array(pkVal15) > 0.5 )/len( pkVal15 ),
        ) )
    ax.hist( pkVal5, bins = 40, alpha = 0.5, label = "5 Sq", histtype =     "step", linewidth = 2, color = None )
    ax.hist( pkVal15, bins = 40, alpha = 0.5, label = "15 Sq", histtype =     "step", linewidth = 2, color = None )
    ax.set_xlabel( "EPSP Amplitude (mV)" )
    ax.set_ylabel( "#")
    ax.text( 0.10, 1.05, "p5={:.3f}".format( prob5 ), fontsize = 16, transform=ax.transAxes )
    ax.text( 0.50, 1.05, "p15={:.3f}".format( prob15 ), fontsize = 16, transform=ax.transAxes )
    ax.text( -0.10, 1.05, str(cell), fontsize = 22, weight = "bold", transform=ax.transAxes )
    ax.legend( loc = 'upper right', frameon = False, title = None )

def main():
    global pulseTrig
    global pulseThresh
    global doPlot

    plt.rcParams.update( {"font.size": 24} )
    fig, ax = plt.subplots( nrows = 3, ncols=1, figsize = (8, 15) )

    # Set up the stimulus timings
    df = pandas.read_hdf( patternData )
    print( df.columns[:50] )
    print( "SIZEOF = ", len( df ) )
    cellList = df['cellID'].unique()
    idx = 0
    for cellIdx, cell in enumerate( cellList ):
        if cell == 4001:
            continue
        if cell == 521:
            alphaTab = np.array( [ dualAlphaFunc(t, ALPHATAU1, 0.018*SAMPLE_FREQ ) for t in range( ALPHAWINDOW ) ] )
            alphaTab = alphaTab / max( alphaTab )
            alphaDelay = int( 0.005 * SAMPLE_FREQ )
            pkDelay = int( 0.015 * SAMPLE_FREQ )
            doPlot = False
        else:
            alphaTab = ALPHA
            alphaDelay = ALPHADELAY
            pkDelay = PKDELAY
            doPlot = False

        dcell = df.loc[df["cellID"] == cell]
        ipat = dcell["patternList"].astype(int)
        patList = ipat.unique()
        pkTime5 = []
        pkTime15 = []
        pkVal5 = []
        pkVal15 = []
        num5 = 0
        num15 = 0
        for pp in patList:
            dpat = dcell.loc[ipat == pp]
            sweepList = dpat['sweep'].unique()
            for ss in sweepList:
                dsweep = dpat.loc[ dpat['sweep'] == ss ]
                seqList = dsweep['exptSeq'].unique()
                for seq in seqList:
                    print( "{}, cell{}, pat{}, sweep{}, seq{}".format( 
                        idx, cell, pp, ss, seq ) )
                    idx += 1
                    foundIdx, foundPks, totNumPulse = doStats( 
                        dsweep.loc[dsweep['exptSeq'] == seq],
                        alphaTab, alphaDelay, pkDelay )
                    if totNumPulse == 0:
                        continue
                    if pp in [46,47,48,49,50]:
                        pkTime5.extend( foundIdx )
                        pkVal5.extend( foundPks )
                        num5 += totNumPulse
                    else:
                        pkTime15.extend( foundIdx )
                        pkVal15.extend( foundPks )
                        num15 += totNumPulse
        analyze( ax[cellIdx], cell, num5, pkVal5, num15, pkVal15 )
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
