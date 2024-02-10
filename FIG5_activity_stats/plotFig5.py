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
EPSPTHRESH = 0.0004 # Threshold for a distinct EPSP pk. In Volts
ALPHATAU1 = 0.005 * SAMPLE_FREQ
ALPHATAU2 = 0.042 * SAMPLE_FREQ
ALPHADELAY = int( 0.010 * SAMPLE_FREQ )
ALPHAWINDOW = int( ALPHATAU2 * 8 )
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

ALPHA = np.array( [ dualAlphaFunc(t, ALPHATAU1, ALPHATAU2 ) for t in range( ALPHAWINDOW ) ] )
ALPHA = ALPHA / max( ALPHA )
longAlpha = np.zeros(NUM_SAMPLES)
longAlpha[:ALPHAWINDOW] += ALPHA

def findPeaks( pkDelay, Vm, width = 0.002 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in PulseTrain:
        idx2 = idx + pkDelay - widthSamples
        pks.append(np.median( Vm[idx2:idx2 + widthSamples*2] ))
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

def parseRow( df, cell ):
    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    alphaTab, alphaDelay, pkDelay = setFittingParams( cell )
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pks = findPeaks( pkDelay, epsp )

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
        ascale = pp - lastPP*longAlpha[int(ALPHATAU1*2) + idx-lastIdx]
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

def panelB_fepspVsTime( ax, fpk5, fpk15 ):
    all5 = np.array( [ ff for ff, ee in fpk5 ] )
    all15 = np.array( [ ff for ff, ee in fpk15 ] )
    mean5 = np.mean( all5, axis = 0 )
    mean15 = np.mean( all15, axis = 0 )
    #pk5 = np.array( pk5 )
    #pk15 = np.array( pk15 )
    #mean5 = np.mean( pk5, axis = 0 )
    #mean15 = np.mean( pk15, axis = 0 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "fEPSP (mV )" )
    ax.set_ylim( -0.01, 0.4 )
    ax.text( -0.20, 1.10, "B", fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelC_epspVsTime( ax, pk5, pk15 ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    mean5 = np.mean( pk5, axis = 0 )
    mean15 = np.mean( pk15, axis = 0 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    ax.text( -0.20, 1.10, "C", fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelD_fepspProbVsTime( ax, fpk5, fpk15 ):
    all5 = np.array( [ ff for ff, ee in fpk5 ] )
    all15 = np.array( [ ff for ff, ee in fpk15 ] )
    #pk5 = np.array( pk5 )
    #pk15 = np.array( pk15 )
    pk5 = (all5 > EPSPTHRESH)
    prob5 = np.sum( pk5, axis = 0 ) * 100 / len( pk5 )
    pk15 = (all15 > EPSPTHRESH)
    prob15 = np.sum( pk15, axis = 0 ) * 100 / len( pk15 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob5, color="blue", s = 10, label = "5 Sq" )
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob15, color="orange", s = 10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "Probability (%)" )
    ax.set_ylim( -5, 105 )
    ax.text( -0.20, 1.1, "D", fontsize = 22, weight="bold", transform=ax.transAxes )

def panelE_epspProbVsTime( ax, pk5, pk15 ):
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
    ax.text( -0.20, 1.10, "E", fontsize = 22, weight="bold", transform=ax.transAxes )

def panelF_fepspVsISI( ax, fpk5, fpk15 ):
    all5 = np.array( [ ff for ff, ee in fpk5 ] )
    all15 = np.array( [ ff for ff, ee in fpk15 ] )
    mean5 = np.mean( all5, axis = 0 )
    mean15 = np.mean( all15, axis = 0 )
    padt = np.pad( PulseTrain, 1)
    isi = PulseTrain - padt[:len( PulseTrain )]
    ax.scatter( isi/ SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( isi/ SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "ISI (s)" )
    ax.set_ylabel( "fEPSP (mV )" )
    ax.set_ylim( -0.01, 0.4 )
    ax.set_xlim( -0.01, 0.2 )
    ax.text( -0.20, 1.10, "F", fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelX_ebyfVsTime( ax, fpk5, fpk15 ):
    all5 = np.array( [ np.array(ee)/np.array(ff) for ff, ee in fpk5 ] )
    all15 = np.array( [ np.array(ee)/np.array(ff) for ff, ee in fpk15 ] )
    mean5 = np.mean( all5, axis = 0 )
    mean15 = np.mean( all15, axis = 0 )
    ax.scatter( PulseTrain/ SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq" )
    ax.scatter( PulseTrain/ SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "epsp/field" )
    ax.text( -0.20, 1.10, "X", fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelG_epspVsISI( ax, pk5, pk15 ):
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
    ax.set_xlabel( "ISI (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    ax.set_xlim( -0.01, 0.2 )
    ax.text( -0.20, 1.10, "G", fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelHI_epspVsField( ax, fpk, label ):
    color = "orange" if label == "I" else "blue"

    allf = []
    alle = []
    for f, e in fpk:
        allf.extend( f )
        alle.extend( e )
    allf = np.array( allf )
    alle = np.array( alle )
    mask = (allf > 0) & (alle > 0 )
    allf = allf[mask]
    alle = alle[mask]
    ax.scatter( allf, alle, color = color, s = 10 )
    fit = linregress( allf, alle )
    # slope, intercept, r-value, p-value, stderr of p_value
    ax.plot( [0.05, 0.35],[fit[0]*0.05+fit[1], fit[0]*0.35+fit[1]], color="black" )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Field EPSP (mV)" )
    ax.set_ylabel( "EPSP (mV)" )
    ax.set_ylim( -0.05, 8 )
    ax.set_xlim( -0.01, 0.35 )
    ax.text( 0.10, 0.95, "r={:.2f}".format(fit[2]), fontsize = 16, transform=ax.transAxes )
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelJ_fpkHisto( ax, fpk5, fpk15 ):
    allf5 = []
    allf15 = []
    for f, e in fpk5:
        allf5.extend( f )
    for f, e in fpk15:
        allf15.extend( f )
    ax.hist( allf5, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, color = "blue" )
    ax.hist( allf15, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, color = "orange" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "field EPSP (mV )" )
    ax.set_ylabel( "#" )
    ax.text( -0.20, 1.10, "J", fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelK_epkHisto( ax, pk5, pk15 ):
    pk5 = np.array( pk5 ).flatten()
    pk5 = pk5[pk5 > EPSPTHRESH]
    pk15 = np.array( pk15 ).flatten()
    pk15 = pk15[pk15 > EPSPTHRESH]
    ax.hist( pk5, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "blue" )
    ax.hist( pk15, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "orange" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "EPSP (mV )" )
    ax.set_ylabel( "#" )
    ax.text( -0.20, 1.10, "K", fontsize = 22, weight = "bold", transform=ax.transAxes )

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
    pks = findPeaks( pkDelay, epsp )
    fitEPSP = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    foundPks = []
    foundIdx = []
    for idx, pp in zip( PulseTrain, pks ):
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
    ax.plot( tepsp, epsp[:len(tepsp)], "b", label = "Data EPSP   " )
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
    ax.text( -0.08, 1.05, "A", fontsize = 22, weight = "bold", transform=ax.transAxes )
    #ax.set_xlabel("Time (s)")

def scanData( df ):
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
        alphaTab, alphaDelay, pkDelay = setFittingParams( cell )
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
                    fpks, epks = parseRow( dseq, cell )
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

def panelLM_varianceHisto( ax1, ax2, patDict ):
    totfvar = []
    totevar = []
    for pattern in [46,47,48,49,50]:
        if len( patDict[pattern] ) > 2:
            ff = np.array([ pp[0] for pp in patDict[pattern] ])
            fvar = np.std( ff, axis = 0 )
            totfvar.extend( fvar )
            ee = np.array([ pp[1] for pp in patDict[pattern] ])
            evar = np.std( ee, axis = 0 )
            totevar.extend( evar )
    ax1.hist( totfvar, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "blue" )
    ax2.hist( totevar, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "blue" )

    totfvar = []
    totevar = []
    for pattern in [52, 53, 55]:
        if len( patDict[pattern] ) > 2:
            ff = np.array([ pp[0] for pp in patDict[pattern] ])
            fvar = np.std( ff, axis = 0 )
            totfvar.extend( fvar )
            ee = np.array([ pp[1] for pp in patDict[pattern] ])
            evar = np.std( ee, axis = 0 )
            totevar.extend( evar )
    totfvar = np.array( totfvar )
    totfvar = totfvar[totfvar > 0]
    ax1.hist( totfvar, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "orange" )
    ax2.hist( totevar, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "orange" )
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_xlabel( "std:field EPSP" )
    ax1.set_ylabel( "#" )
    ax1.text( -0.20, 1.05, "L", fontsize = 22, weight = "bold", transform=ax1.transAxes )
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xlabel( "std:EPSP" )
    ax2.set_ylabel( "#" )
    ax2.text( -0.20, 1.05, "M", fontsize = 22, weight = "bold", transform=ax2.transAxes )

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

    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (10,24) )
    gs = fig.add_gridspec( 7, 2 ) # 4 rows, 2 cols
    #fig, ax = plt.subplots( nrows = 3, ncols=3, figsize = (18, 15) )

    # Set up the stimulus timings
    df = pandas.read_hdf( patternData )
    ax = fig.add_subplot( gs[0,:] )
    dcell4041 = df.loc[df["cellID"] == 4041]
    panelA_SampleTrace( ax, dcell4041 )
    pk5, pk15, fpk5, fpk15, patDict = scanData( df )
    panelB_fepspVsTime( fig.add_subplot(gs[1,0]), fpk5, fpk15 )
    panelC_epspVsTime( fig.add_subplot(gs[1,1]), pk5, pk15 )
    panelD_fepspProbVsTime( fig.add_subplot(gs[2,0]), fpk5, fpk15 )
    panelE_epspProbVsTime( fig.add_subplot(gs[2,1]), pk5, pk15 )
    panelF_fepspVsISI(  fig.add_subplot(gs[3,0]), fpk5, fpk15 )
    panelG_epspVsISI(  fig.add_subplot(gs[3,1]), pk5, pk15 )
    panelHI_epspVsField( fig.add_subplot(gs[4,0]), fpk5, "H" )
    panelHI_epspVsField( fig.add_subplot(gs[4,1]), fpk15, "I" )
    panelJ_fpkHisto( fig.add_subplot(gs[5,0]), fpk5, fpk15 )
    panelK_epkHisto( fig.add_subplot(gs[5,1]), pk5, pk15 )
    panelLM_varianceHisto( fig.add_subplot(gs[6,0]), fig.add_subplot(gs[6,1]), patDict )
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
