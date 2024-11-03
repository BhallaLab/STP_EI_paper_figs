import pandas
import pylab
import numpy as np
import argparse
import math
from scipy.stats import linregress
from scipy.stats import kstest

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
exptFile = "../../../2022/VC_DATA/all_cells_SpikeTrain_CC_long.h5"
#simFile = "SIMDATA47/simData_orig_0.h5"
simFile = "SIMDATA49/simData_wtGABA_10.h5"
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 11
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
EPSPTHRESH = 0.0004 # Threshold for a distinct EPSP pk. In Volts
'''
PKDELAY = int( 0.020 * SAMPLE_FREQ )
ALPHATAU1 = 0.005 * SAMPLE_FREQ
ALPHATAU2 = 0.042 * SAMPLE_FREQ
ALPHADELAY = int( 0.010 * SAMPLE_FREQ )
'''
ALPHAWINDOW = int( 0.32 * SAMPLE_FREQ )
MINIMUM_PULSE_THRESH = 4e-4

PulseTrain = np.array([4001,10684,11276,11603,13433,15914,16193,17131,19457,19827,20561,21153,21578,
    22460,24407,24665,25093,25667,26213,26726,27343,28046,28625,29322,29608,31223,31729,32400,32756,
    33317,33897,35890,36496,36986,37267,37484,38755,39890,40495,41873,42970,43399,45768,46100,46695,
    46931,47430,47639,47877,48568,49189,51579,52910,53373,53643,56169,56686,57112,57467,57834,58721,
    59254,60261,60473,61816,63607,64798,66090,66291,69446,70416,70666,70898,71145,71821,72805,73201,
    74279,74777,75520,76181,77447,77966,78309,79050,79331,80383,81575,82380,82991,85548,87622,88515,
    88839,89510,89866,90977,91257,91841,92837,93249,94872,95549,96164,96975,98498,99152,99545,99795,
    100493,101582,102149,103757,107075,107600,107969,108705,109143,109875,110347,110856,113988,114470,
    115634,116946,117489,118060,119694,121243,122078,122580,124326,125053,127211,128234,128814,129380,
    129945,130884,131133,131550,132432,133262,133560,134345,134707,135065,135938,136529,137450,137806,
    139055,140234,141304,143221,143573,144296,145640,145984,146846,147856,148671,150909,152493,152852,
    153268,153931,155048,155690,156475,157345,158850,159443,159768,160600,160919,161424,161660,161956,
    163448,163758,164107,165661,166052,166540,167119,168032,169773,170130,171780,172502,173106,174142,
    174728,175182,175694,176340,177236,178437,179524,180446,183258,183781,185319,187213,189396,190365,
    190837,191267,191619,192282,192848,193144,193689,194521,195822,196751,197884,199981,200689,201095,
    202108,203280,204018,205585,206552,207234,207796,209126,209832])


def alphaFunc( t, tp ):
    return (t/tp) * np.exp(1-t/tp)

def dualAlphaFunc( t, t1, t2 ):
    if t < 0:
        return 0.0
    if abs( t1 - t2 ) < 1e-6:
        return alphaFunc( t, t1 )
    return (1.0/(t1-t2)) * (np.exp(-t/t1) - np.exp(-t/t2))

def findPeaks( pkDelay, Vm, width = 0.002, threshold = 0.0 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in PulseTrain:
        idx2 = idx + pkDelay - widthSamples
        vv = np.median( Vm[idx2:idx2 + widthSamples*2] )
        if vv > threshold:
            pks.append( vv )
        else:
            pks.append( 0 )
    return pks

def findMinima( pkDelay, Vm, width = 0.002 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    lastPP = 0.0
    pks = []
    for idx in PulseTrain:
        idx2 = idx + pkDelay - widthSamples
        pks.append(-min( Vm[idx2:idx2 + widthSamples*2] ))

    return pks

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

def parseRow( df, cell, args ):
    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    epsp += np.random.normal(0.0, args.noise*1e-3, epsp.shape )
    if cell == 0:
        epsp *= 1000     # Scale to mV.
        args.threshold = 0.2 # Assign for sim.
    baseline = min( np.percentile(epsp, 25 ), np.mean( epsp[100:1000 ]) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pks = findPeaks( pkDelay, epsp, threshold=args.threshold )

    if cell == 4041:
        field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
        field = fftFilter( field )
        fpks = findMinima( 100, field, width = 0.002 )
    else:
        field = np.zeros( len( epsp ) )
        fpks = np.zeros(len( PulseTrain ))

    fitEPSP = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    epks = []
    for idx, pp in zip( PulseTrain, pks ):
        ii = idx + alphaDelay
        ascale = pp - lastPP*longAlpha[int(alphaTau1*2) + idx-lastIdx]
        if ascale > 0:
            fitEPSP[ii:ii+ALPHAWINDOW] += ascale * alphaTab
            lastIdx = idx
            lastPP = pp
            epks.append( ascale )
        else:
            epks.append( 0.0 )

    if doFieldFitPlot:
        plt.figure( figsize = (30, 5 ))
        plt.plot( tepsp, field )
        print( "LENS = ", len( fpks ), len(PulseTrain), np.median( field ) )
        plt.scatter( np.array(PulseTrain)/SAMPLE_FREQ, 
                -np.array(fpks) -np.median(field), c = "red", s = 20 )
        plt.show()

    assert( len( fpks ) == len( epks ) )
    assert( len( fpks ) == len( PulseTrain ) )
    return fpks, epks 

def panelC_probVsTime( ax, column, pk5, pk15 ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    pk5 = (pk5 > EPSPTHRESH)
    prob5 = np.sum( pk5, axis = 0 ) * 100 / len( pk5 )
    pk15 = (pk15 > EPSPTHRESH)
    prob15 = np.sum( pk15, axis = 0 ) * 100 / len( pk15 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob5, color="blue", s = 10, label = "5 Sq" )
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob15, color="orange", s = 10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "Probability (%)" )
    ax.set_ylim( -5, 105 )
    label = chr( ord("C") + column )
    ax.text( -0.20, 1.06, label, fontsize = 22, weight="bold", transform=ax.transAxes )
    subtitle = "Experiment" if column == 0 else "Model"
    ax.text( 0.20, 1.06, subtitle, fontsize = 16, transform=ax.transAxes )

def panelE_epspVsTime( ax, column, pk5, pk15 ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    mean5 = np.mean( pk5, axis = 0 )
    mean15 = np.mean( pk15, axis = 0 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_ylim( -0.2, 5.1 )
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    label = chr( ord("E") + column )
    ax.text( -0.20, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelG_epspVsISI( ax, column, pk5, pk15 ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    mean5 = np.mean( pk5, axis = 0 )
    mean15 = np.mean( pk15, axis = 0 )
    padt = np.pad( PulseTrain, 1)
    isi = PulseTrain - padt[:len( PulseTrain )]
    ax.scatter( isi / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( isi / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_ylim( -0.2, 5.1 )
    ax.set_xlabel( "ISI (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    ax.set_xlim( -0.01, 0.2 )
    label = chr( ord("G") + column )
    ax.text( -0.20, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelI_epkHisto( ax, column, pk5, pk15 ):
    pk5 = np.array( pk5 ).flatten()
    pk5 = pk5[pk5 > EPSPTHRESH]
    pk15 = np.array( pk15 ).flatten()
    pk15 = pk15[pk15 > EPSPTHRESH]
    ax.hist( pk5, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "blue" )
    ax.hist( pk15, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "orange" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim( -0.5, 12 )
    ax.set_xlabel( "EPSP (mV )" )
    ax.set_ylabel( "#" )
    label = chr( ord("I") + column )
    ax.text( -0.20, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return [pk5, pk15]

def panelA_SampleTrace( ax, dcell, column, args ):
    ipat = dcell["patternList"].astype(int)
    #df = dcell.loc[(ipat == 46) & (dcell['sweep'] == 12)]
    df = dcell.loc[(ipat == 46)]
    cell = df['cellID'].unique()[0]
    #print( df.shape )
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    PLOTLEN = 6.0
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    epsp += np.random.normal(0.0, args.noise*1e-3, epsp.shape )
    if cell == 0:
        epsp *= 1000
    else:
        pulseTrig *= 100
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )

    #baseline = min( np.percentile(epsp, 25 ), 0.0 )
    baseline = min( np.percentile(epsp, 25 ), np.mean( epsp[100:1000 ]) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pks = findPeaks( pkDelay, epsp, threshold=args.threshold )
    fitEPSP = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    foundPks = []
    foundIdx = []
    for idx, pp in zip( PulseTrain, pks ):
        ii = idx + alphaDelay
        ascale = pp - lastPP*longAlpha[int(alphaTau1*2) + idx-lastIdx]
        if ascale > 0:
            fitEPSP[ii:ii+ALPHAWINDOW] += ascale * alphaTab
            lastIdx = idx
            lastPP = pp
            foundPks.append( ascale )
            foundIdx.append( idx )

    runtime = settleTime + SAMPLE_TIME + postStim

    tepsp = tepsp[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsp)] -5
    print( "MEAN = ", np.mean( epsp ) )
    ax.plot( tepsp, epsp[:len(tepsp)], "b", label = "Data EPSP   " )
    ax.plot( tepsp, fitEPSP[:len(tepsp)], "r", label = "Fit EPSP" )
    #ax.plot( tepsp, (field[:len(tepsp)] - 0.2) * 10, "m", label = "Field" )
    ax.plot( tepsp, pt, "g", label = "Trigger" )
    #ax.plot( [0,0,0.25], [10,5,5], color="black", linewidth=2.5 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    '''
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    '''
    ax.set_ylabel( "EPSP (mV)" )
    ax.set_xlabel( "Time (s)" )
    ax.legend( loc = "upper right", ncol=3, frameon=False, fontsize = 14 )
    ax.set_xlim( 0, PLOTLEN )
    #ax.set_ylim( -7, 10 )
    label = chr( ord("A") + column )
    ax.text( -0.10, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    #ax.set_xlabel("Time (s)")

def scanData( df, args ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    pk5 = []
    pk15 = []
    fpk5 = []
    fpk15 = []
    patDict = { pp:[] for pp in [46,47,48,49,50,52,53,55] } # 8 patterns.
    for cellIdx, cell in enumerate( cellList ):
        if cell == 4001:
            continue
        alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
        dcell = df.loc[df["cellID"] == cell]
        ipat = dcell["patternList"].astype(int)
        patList = ipat.unique()
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
                    fpks, epks = parseRow( dseq, cell, args )
                    #print( "SUM = ", sum( foundIdx ) )
                    patDict[pp].append( [fpks, epks] )
                    if pp in [46,47,48,49,50]:
                        pk5.append(epks)
                        if cell == 4041: # The only one with field data
                            print( "Appending 4041 ", len( fpks ), len( epks ) )
                            fpk5.append( [fpks, epks] )
                    else:
                        pk15.append(epks)
                        if cell == 4041: # The only one with field data
                            fpk15.append( [fpks, epks] )

    return pk5, pk15, fpk5, fpk15, patDict

def panelK_varianceHisto( ax, column, patDict ):
    totevar = []
    for pattern in [46,47,48,49,50]:
        if len( patDict[pattern] ) > 2:
            ee = np.array([ pp[1] for pp in patDict[pattern] ])
            evar = np.std( ee, axis = 0 )
            totevar.extend( evar )
    ax.hist( totevar, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "blue" )
    totevar5 = totevar

    totevar = []
    for pattern in [52, 53, 55]:
        if len( patDict[pattern] ) > 2:
            ee = np.array([ pp[1] for pp in patDict[pattern] ])
            evar = np.std( ee, axis = 0 )
            totevar.extend( evar )
    ax.hist( totevar, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "orange" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim( -0.5, 4.2 )
    ax.set_xlabel( "std:EPSP" )
    ax.set_ylabel( "#" )
    label = chr( ord( "K" ) + column )
    ax.text( -0.20, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return [totevar5, totevar]

def setFittingParams( cell ):
    alphaTau1 = 0.005 * SAMPLE_FREQ
    alphaTau2 = 0.042 * SAMPLE_FREQ
    alphaDelay = int( 0.01 * SAMPLE_FREQ )
    pkDelay = int( 0.020 * SAMPLE_FREQ )
    if cell == 521:
        alphaTau2 = 0.018 * SAMPLE_FREQ
        alphaDelay = int( 0.005 * SAMPLE_FREQ )
        pkDelay = int( 0.015 * SAMPLE_FREQ )
    elif cell == 0: # Simulated
        alphaTau1 = 0.004 * SAMPLE_FREQ
        alphaTau2 = 0.008 * SAMPLE_FREQ
        alphaDelay = int( 0.002 * SAMPLE_FREQ )
        pkDelay = int( 0.008 * SAMPLE_FREQ )
    else: # Use defaults from above
        pass

    alphaTab = np.array( [ dualAlphaFunc(t, alphaTau1, alphaTau2 ) for t in range( ALPHAWINDOW ) ] )
    alphaTab = alphaTab / max( alphaTab )
    return alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2

def main():
    parser = argparse.ArgumentParser( description = "Read and plot sim data" )
    parser.add_argument( "--fname", type = str, help = "Required: Name of hdf5 pandas file with data.", default = exptFile )
    parser.add_argument( "-f2", "--fname2", type = str, help = "Optional: Name of another hdf5 pandas file with data.", default = simFile )
    parser.add_argument( "-n", "--noise", type = float, help = "Optional: Noise to add to EPSP trace, in mV. Default = 0.", default = 0 )
    parser.add_argument( "-t", "--threshold", type = float, help = "Optional: Threshold for classifying an EPSP trace as an event. In mV. Default = 0.", default = 0 )
    args = parser.parse_args()
    df = pandas.read_hdf( args.fname )
    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (10,24) )
    gs = fig.add_gridspec( 7, 2 ) # 7 rows, 2 cols
    ehisto, vhisto = plotFrame( gs, fig, args, df, 0, 0 )
    if args.fname2:
        df2 = pandas.read_hdf( args.fname2 )
        # Here we check if it is real or synth data
        eh2, vh2 = plotFrame( gs, fig, args, df2, 1, 0 )

        stat5, pval5 = kstest( ehisto[0], eh2[0] )
        stat15, pval15 = kstest( ehisto[1], eh2[1] )
        print ( "KS Tests: EPSPs. 5 Sq Stat = {:.3f}, pval = {:.4f}, 15Sq: {:.3f}, {:.4f}".format( stat5, pval5, stat15, pval15 ) )

        stat5, pval5 = kstest( vhisto[0], vh2[0] )
        stat15, pval15 = kstest( vhisto[1], vh2[1] )
        print ( "KS Tests: EPSPVars. 5sq: Stat = {:.3f}, pval = {:.4f}, 15Sq: {:.3f}, {:.4f}".format( stat5, pval5, stat15, pval15 ) )
    fig.tight_layout()
    plt.show()


def plotFrame(gs, fig, args, df, column = 0, cell = 0):
    global pulseTrig
    global pulseThresh

    # Set up the stimulus timings
    ax = fig.add_subplot( gs[column,:] ) # Hack: col 0 comes in row 0, col 1 in row 2
    panelA_SampleTrace( ax, df, column, args )
    #dcell4041 = df.loc[df["cellID"] == 4041]
    pk5, pk15, fpk5, fpk15, patDict = scanData( df, args )
    print( "LENGTHS = ", len( pk5 ), len( pk15 ), len( fpk5 ), len( fpk15 ) )
    subtitle = "Experiment" if column == 0 else "Model"
    panelC_probVsTime( fig.add_subplot(gs[2,column]), column, pk5, pk15 )
    panelE_epspVsTime( fig.add_subplot(gs[3,column]), column, pk5, pk15 )
    #panelD_fepspVsISI(  fig.add_subplot(gs[2,0]), fpk5, fpk15 )
    panelG_epspVsISI(  fig.add_subplot(gs[4,column]), column, pk5, pk15 )
    #panelFG_epspVsField( fig.add_subplot(gs[3,0]), fpk5, "F" )
    #panelFG_epspVsField( fig.add_subplot(gs[3,1]), fpk15, "G" )
    #panelH_fpkHisto( fig.add_subplot(gs[4,0]), fpk5, fpk15 )
    ehisto = panelI_epkHisto( fig.add_subplot(gs[5,column]), column, pk5, pk15 )
    vhisto = panelK_varianceHisto( fig.add_subplot(gs[6,column]), column, patDict )
    return ehisto, vhisto
    

if __name__ == "__main__":
    main()
