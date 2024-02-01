import pandas
import pylab
import numpy as np
import math
from scipy.stats import linregress

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
doFieldPlot = False
doEpspVsFieldPlot = True
doFieldFitPlot = False

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

def findMinima( pulseTrig, pulseThresh, pkDelay, Vm, width = 0.002 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    lastPP = 0.0
    pks = []
    pkIdx = []
    for idx, pp in enumerate( pulseTrig ):
        if pp > pulseThresh and lastPP < pulseThresh:
            idx2 = idx + pkDelay - widthSamples
            pkIdx.append( idx )
            pks.append(-min( Vm[idx2:idx2 + widthSamples*2] ))
        lastPP = pp

    return pkIdx, pks

def medianWindowFilter(data):
    WINDOW_SIZE = 201
    PAD_WIDTH = WINDOW_SIZE // 2
    filtered = doFilter( data )
    padData = np.pad( data, PAD_WIDTH, mode = 'edge' )
    filtered = []
    for i in range( len( data ) ):
        window = padData[i:i+WINDOW_SIZE]
        filtered.append(np.median(window))
    data = np.array(data) - np.array(filtered)
    return data

def fftFilter( data ):
    cut_off_frequency = 0.0001
    fft_data = np.fft.fft(data)
    N = len(data)
    frequencies = np.fft.fftfreq(N, d=1)  # d is the sample spacing; assumed to be 1 for simplicity
    filter_mask = np.abs(frequencies) > cut_off_frequency
    fft_data[filter_mask] = 0
    filtered_data = np.real(np.fft.ifft(fft_data))
    '''
    plt.figure(figsize=(15, 5))
    #plt.plot(filtered_data, label='Filtered Data', linewidth=2)
    #plt.plot(data, label='Original Data')
    data = np.array( data - filtered_data )
    plt.plot(data, label='Original Data')
    plt.legend()
    plt.show()
    '''
    return np.array(data - filtered_data)

def parseRow( df, cell ):
    alphaTab, alphaDelay, pkDelay = setFittingParams( cell )
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pkIdx, pks = findPeaks( pulseTrig, pulseThresh, pkDelay, epsp )

    if cell == 4041:
        field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
        field = fftFilter( field )
        fpkIdx, fpks = findMinima( pulseTrig, pulseThresh, 100, field, width = 0.002 )
    else:
        field = np.zeros( len( epsp ) )
        fpkIdx = pkIdx
        fpks = np.zeros(len(fpkIdx))


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

    if doFieldFitPlot:
        plt.figure( figsize = (30, 5 ))
        plt.plot( tepsp, field )
        print( "LENS = ", len( fpks ), len(pkIdx), np.median( field ) )
        plt.scatter( np.array(pkIdx)/SAMPLE_FREQ, 
                -np.array(fpks) -np.median(field), c = "red", s = 20 )
        plt.show()

    return pkIdx, fpks, foundIdx, foundPks 

def panelBC_pkVsTime( ax1, ax2, df ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    for cellIdx, cell in enumerate( cellList ):
        if cell == 4001:
            continue
        alphaTab, alphaDelay, pkDelay = setFittingParams( cell )
        dcell = df.loc[df["cellID"] == cell]
        ipat = dcell["patternList"].astype(int)
        patList = ipat.unique()
        pk5 = {}
        pk15 = {}
        for pp in patList:
            dpat = dcell.loc[ipat == pp]
            sweepList = dpat['sweep'].unique()
            for ss in sweepList:
                dsweep = dpat.loc[ dpat['sweep'] == ss ]
                seqList = dsweep['exptSeq'].unique()
                for seq in seqList:
                    dseq = dsweep.loc[dsweep['exptSeq'] == seq]
                    print( "{}, cell{}, pat{}, sweep{}, seq{}".format( 
                        idx, cell, pp, ss, seq ) )
                    idx += 1
                    pkIdx, fpks, foundIdx, foundPks = parseRow( dseq, cell )
                    #print( "SUM = ", sum( foundIdx ) )
                    print( foundIdx[1:10] )
                    if len( pkIdx ) == 0:
                        continue
                    if pp in [46,47,48,49,50]:
                        pk5[pp] = [foundIdx, foundPks, pkIdx]
                    else:
                        pk15[pp] = [foundIdx, foundPks, pkIdx]
    
        cellStats[cell] = [pk5, pk15]
    return

def panelA_SampleTrace( ax, dcell ):
    ipat = dcell["patternList"].astype(int)
    df = dcell.loc[(ipat == 46) & (dcell['sweep'] == 12)]
    alphaTab, alphaDelay, pkDelay = setFittingParams( 4041 )
    PLOTLEN = 2.0
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )

    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pkIdx, pks = findPeaks( pulseTrig, pulseThresh, pkDelay, epsp )
    fpkIdx, fpks = findPeaks( pulseTrig, pulseThresh, 90, np.median(field)-field, width = 0.001 )
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

    tepsp = tepsp[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsp)] - 0.06
    ax.plot( tepsp, epsp[:len(tepsp)], "b", label = "Data EPSP" )
    ax.plot( tepsp, fitEPSP[:len(tepsp)], "r", label = "Fit EPSP" )
    ax.plot( tepsp, (field[:len(tepsp)] - 0.2) * 10, "m", label = "Field" )
    ax.plot( tepsp, pt * 100, "g", label = "Trigger" )
    ax.plot( [0,0,0.25], [10,5,5], color="black", linewidth=2.5 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    #ax.set_ylabel( "EPSP (mV)" )
    ax.legend( loc = "upper right", frameon = False, fontsize = 14 )
    ax.set_xlim( 0, PLOTLEN )
    ax.set_ylim( -7, 10 )
    #ax.set_xlabel("Time (s)")

def doStats( df, alphaTab, alphaDelay, pkDelay, pattern = None ):
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    if doFieldPlot and pattern == 46:
        plt.figure( figsize = (30, 5 ))
        plt.plot( field )
        plt.plot( pulseTrig * 5 + 0.1 )
        plt.show()

    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    print( "baseline={:.3f}, pulseThresh = {:.6f}".format( 
        baseline, pulseThresh ) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, SAMPLE_TIME, len(epsp) )
    pkIdx, pks = findPeaks( pulseTrig, pulseThresh, pkDelay, epsp )
    fpkIdx, fpks = findPeaks( pulseTrig, pulseThresh, 90, np.median(field)-field, width = 0.001 )
    if doEpspVsFieldPlot and pattern == 46:
        sweep = int( df['sweep'] )
        print( sweep )
        if sweep == 4:
            plt.figure( figsize = (8, 8 ))
        plt.scatter( fpks, pks )
        ret = linregress(fpks,pks)
        for rr in ret:
            # slope, intercept, r-value, p-value, stderr of p_value
            print( "{:.3f}  ".format( rr ), end = "" )
        print()
        if sweep == 20:
            plt.ylabel( "EPSP (mV)" )
            plt.xlabel("fEPSP (mV)")
            plt.show()
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

    if doPlot and pattern == 46:
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

    return np.array(foundIdx), np.array(foundPks), np.array(pkIdx)

def analyze( ax1, ax2, ax3, cell, pks, label ):
    totPks = np.array( [] )
    leftPks = np.array( [] )
    rightPks = np.array( [] )
    leftProb = 0.0
    rightProb = 0.0
    totProb = 0.0
    # Total duration is 11 sec. We split halfway.
    for key in pks:
        ii, pp, nn = pks[key]
        leftPks = np.append( leftPks,  pp[ ii < 5.5*SAMPLE_FREQ] )
        rightPks = np.append( rightPks,  pp[ ii >= 5.5*SAMPLE_FREQ] )
        totPks = np.append( totPks,  pp )
        # 0.5 mV is cutoff for using EPSP.
        leftProb += 2*len( pp[(pp>0.5) & (ii<5.5*SAMPLE_FREQ)] ) / nn
        rightProb += 2*len( pp[(pp>0.5) & (ii>=5.5*SAMPLE_FREQ)] ) / nn
        totProb += len(pp) / nn
        if key <= 50:
            ax3.hist( pp, bins = 20, alpha = 0.5, label = str(key) + " " + label, histtype = "step", linewidth = 2, color = None )
    print("PROBS = tot{:.3f} L{:.3f} R{:.3f}".format(
        totProb / len(pks), leftProb / len(pks), rightProb/len(pks)
    ) )
   
    ax1.hist( totPks, bins = 32, alpha = 0.5, label = label, histtype = "step", linewidth = 2, color = None )
    ax2.hist( leftPks, bins = 32, alpha = 0.5, label = "L " + label, histtype = "step", linewidth = 2, color = None )
    ax2.hist( rightPks, bins = 32, alpha = 0.5, label = "R " + label, histtype = "step", linewidth = 2, color = None )
    ax1.set_xlabel( "EPSP Amplitude (mV)" )
    ax1.set_ylabel( "#")
    #ax1.text( 0.10, 1.05, "p5={:.3f}".format( prob5 ), fontsize = 16, transform=ax.transAxes )
    #ax1.text( 0.50, 1.05, "p15={:.3f}".format( prob15 ), fontsize = 16, transform=ax.transAxes )
    ax1.text( -0.10, 1.05, str(cell), fontsize = 22, weight = "bold", transform=ax1.transAxes )
    ax1.legend( loc = 'upper right', frameon = False, title = None )
    ax2.legend( loc = 'upper right', frameon = False, title = None )
    ax3.legend( loc = 'upper right', frameon = False, title = None )
    


    '''
    # pk5[pp] = [foundIdx, foundPks, totNumPulse]
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
    '''

def setFittingParams( cell ):
    if cell == 521:
        alphaTab = np.array( [ dualAlphaFunc(t, ALPHATAU1, 0.018*SAMPLE_FREQ ) for t in range( ALPHAWINDOW ) ] )
        alphaTab = alphaTab / max( alphaTab )
        alphaDelay = int( 0.005 * SAMPLE_FREQ )
        pkDelay = int( 0.015 * SAMPLE_FREQ )
    else:
        alphaTab = ALPHA
        alphaDelay = ALPHADELAY
        pkDelay = PKDELAY
    return alphaTab, alphaDelay, pkDelay

def main():
    global pulseTrig
    global pulseThresh

    plt.rcParams.update( {"font.size": 24} )
    fig = plt.figure( figsize = (10,15) )
    gs = fig.add_gridspec( 4, 2 ) # 4 rows, 2 cols
    #fig, ax = plt.subplots( nrows = 3, ncols=3, figsize = (18, 15) )

    # Set up the stimulus timings
    df = pandas.read_hdf( patternData )
    cellList = df['cellID'].unique()
    dcell4041 = df.loc[df["cellID"] == 4041]
    ax = fig.add_subplot( gs[0,:] )
    panelA_SampleTrace( ax, dcell4041 )
    ax1 = fig.add_subplot( gs[1,0] )
    ax2 = fig.add_subplot( gs[1,1] )
    panelBC_pkVsTime( ax1, ax2, df )
    '''
    idx = 0
    cellStats = {}
    for cellIdx, cell in enumerate( cellList ):
        if cell == 4001:
            continue
        alphaTab, alphaDelay, pkDelay = setFittingParams( cell )
        dcell = df.loc[df["cellID"] == cell]
        ipat = dcell["patternList"].astype(int)
        patList = ipat.unique()
        pk5 = {}
        pk15 = {}
        for pp in patList:
            dpat = dcell.loc[ipat == pp]
            sweepList = dpat['sweep'].unique()
            for ss in sweepList:
                dsweep = dpat.loc[ dpat['sweep'] == ss ]
                seqList = dsweep['exptSeq'].unique()
                for seq in seqList:
                    dseq = dsweep.loc[dsweep['exptSeq'] == seq]
                    if cell == 4041 and pp == 46 and ss == 4:
                        ax = fig.add_subplot( gs[0,:] )
                        plotSampleTrace( ax, dseq, alphaTab, alphaDelay, pkDelay )
                    print( "{}, cell{}, pat{}, sweep{}, seq{}".format( 
                        idx, cell, pp, ss, seq ) )
                    break
                    idx += 1
                    foundIdx, foundPks, pkIdx = doStats( 
                        dseq, alphaTab, alphaDelay, pkDelay, pattern = pp)
                    if len( pkIdx ) == 0:
                        continue
                    if pp in [46,47,48,49,50]:
                        pk5[pp] = [foundIdx, foundPks, pkIdx]
                    else:
                        pk15[pp] = [foundIdx, foundPks, pkIdx]
    
        cellStats[cell] = [pk5, pk15]
        #analyze( ax[0][cellIdx], ax[1][cellIdx], ax[2][cellIdx], cell, pk5, "5 Sq" )
        #analyze( ax[0][cellIdx], ax[1][cellIdx], ax[2][cellIdx], cell, pk15, "15 Sq")
    '''
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
