import pandas
import pylab
import numpy as np
import math
import argparse
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
#patternData = "simData.h5"
#patternData = "simData_12.h5"
#patternData = "simData_11.h5"
#patternData = "simData_111.h5"
#patternData = "simData_222.h5"
patternData = "simData_333.h5"
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 5
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
EPSPTHRESH = 0.0004 # Threshold for a distinct EPSP pk. In Volts
ALPHAWINDOW = int( 0.32 * SAMPLE_FREQ )
MINIMUM_PULSE_THRESH = 4e-4
PulseTrain = {}


def alphaFunc( t, tp ):
    return (t/tp) * np.exp(1-t/tp)

def dualAlphaFunc( t, t1, t2 ):
    if t < 0:
        return 0.0
    if abs( t1 - t2 ) < 1e-6:
        return alphaFunc( t, t1 )
    return (1.0/(t1-t2)) * (np.exp(-t/t1) - np.exp(-t/t2))


def findPeaks( pkDelay, Vm, pulses, width = 0.004 ):
    if len( Vm ) < NUM_SAMPLES or len( pulses ) != 33:
        return []
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in pulses:
        idx2 = idx + pkDelay - widthSamples
        #print( idx, idx2, len( Vm ) )
        pks.append(np.quantile( Vm[idx2:idx2 + widthSamples*2], 0.90 ))
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

def parseRow( df, cell, ff ):
    global PulseTrain
    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    if cell == 0:
        baseline = min(epsp )
    else:
        baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    if sum( edges ) != 33:
        return []
    PulseTrain[ff] = np.arange( 0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    pks = findPeaks( pkDelay, epsp, PulseTrain[ff] )

    fitEPSP = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    epks = []
    for idx, pp in zip( PulseTrain[ff], pks ):
        ii = idx + alphaDelay
        ascale = pp - lastPP*longAlpha[int(alphaTau1*2) + idx-lastIdx]
        if ascale > 0:
            fitEPSP[ii:ii+ALPHAWINDOW] += ascale * alphaTab
            lastIdx = idx
            lastPP = pp
            epks.append( ascale )
        else:
            epks.append( 0.0 )

    return epks 

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
    #ax.set_ylim( -0.01, 0.4 )
    ax.text( -0.20, 1.10, "B", fontsize = 22, weight = "bold", transform=ax.transAxes )


def panelC_epspVsTime( ax, pk5, ff, idx ):
    pk5 = np.array( pk5 )
    pk5[pk5>7] = np.median( pk5 )   # filter out outliers.
    mean5 = np.mean( pk5, axis = 0 )
    med5 = np.median( pk5, axis = 0 )
    t0 = PulseTrain[ff][:9]/SAMPLE_FREQ
    t1 = PulseTrain[ff][9:17]/SAMPLE_FREQ
    t2 = PulseTrain[ff][17:25]/SAMPLE_FREQ
    t3 = PulseTrain[ff][25:]/SAMPLE_FREQ
    wid = 0.5/ff
    for pp in pk5:
        ax.scatter( t0 + np.random.rand( len( t0 ) ) * wid, pp[:9], color="cyan", s=1, label = "5 Sq" )
        ax.scatter( t1 + np.random.rand( len( t1 ) ) * wid, pp[9:17], color="yellow", s=1, label = "5 Sq" )
        ax.scatter( t2 + np.random.rand( len( t2 ) ) * wid, pp[17:25], color="pink", s=1, label = "5 Sq" )
        ax.scatter( t3 + np.random.rand( len( t3 ) ) * wid, pp[25:], color="palegreen", s=1, label = "5 Sq" )
    #ax.scatter( PulseTrain[ff / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    #ax.scatter( PulseTrain[ff] / SAMPLE_FREQ, mean5, color="blue", s=10, label = "Means" )
    ax.plot( PulseTrain[ff] / SAMPLE_FREQ, mean5, color="blue", markersize=10, label = "Mean" )
    ax.plot( PulseTrain[ff] / SAMPLE_FREQ, med5, color="red", markersize=10, label = "Median" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "EPSC ratio" )
    #ax.set_ylim( 0, 2.5 )
    label = chr( ord("B")+idx )
    ax.text( 0.05, 0.90, str(ff)+" Hz", fontsize = 16, transform=ax.transAxes )
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

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
    #ax.set_ylim( -0.01, 0.4 )
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
    ax.set_ylabel( "EPSC (pA)" )
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
    #print( "ALLLEEEY, ", len( allf ), len( alle ), np.mean( allf), np.mean( alle) )
    mask = (allf > 0) & (alle > 0 )
    allf = allf[mask]
    alle = alle[mask]
    ax.scatter( allf, alle, color = color, s = 10 )
    fit = linregress( allf, alle )
    # slope, intercept, r-value, p-value, stderr of p_value
    ax.plot( [0.0, 25],[fit[1], fit[0]*25+fit[1]], color="black" )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Field EPSP (mV)" )
    ax.set_ylabel( "EPSC (pA)" )
    #ax.set_ylim( -0.05, 8 )
    #ax.set_xlim( -0.01, 0.35 )
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
    ax.set_xlabel( "EPSC (pA )" )
    ax.set_ylabel( "#" )
    ax.text( -0.20, 1.10, "K", fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelA_SampleTrace( ax, dcell ):
    #df = dcell.loc[(dcell['sweep'] == 0)]
    df = dcell.loc[(dcell['stimFreq'] == 20)]
    #print( "PPPPPP STARET ", df['probePulseStart'].unique(), df['sweepLength'].unique() )
    #print( "fct", dcell['stimFreq'] )
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams(0)
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    PLOTLEN = 2.5
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    trig2 = np.array(df.iloc[0, SAMPLE_START+NUM_SAMPLES:SAMPLE_START+2*NUM_SAMPLES ])
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    pulses = np.arange(0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    #print( PulseTrain )

    if dcell['cellID'][0] == 0:
        baseline = min(epsp )
    else:
        baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, SAMPLE_TIME, len(epsp) )
    pks = findPeaks( pkDelay, epsp, pulses )
    #print( "NUMPKS = ", len( pks ) )
    fitepsp = np.zeros( len( epsp ) + ALPHAWINDOW * 2 )
    lastPP = 0.0
    lastIdx = 0
    foundPks = []
    foundIdx = []
    for idx, pp in zip( pulses, pks ):
        ii = idx + alphaDelay
        ascale = pp - lastPP*longAlpha[int(alphaTau1*2) + idx-lastIdx]
        if ascale > 0:
            fitepsp[ii:ii+ALPHAWINDOW] += ascale * alphaTab
            lastIdx = idx
            lastPP = pp
            foundPks.append( ascale )
            foundIdx.append( idx )

    tepsp = tepsp[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsp)]
    ax.plot( tepsp, epsp[:len(tepsp)], "b", label = "Data epsp   " )
    #ax.plot( tepsp, trig2[:len(tepsp)]*10 - 20, "m", label = "trig2" )
    ax.plot( tepsp, fitepsp[:len(tepsp)], "r", label = "Fit epsp" )
    #ax.plot( tepsp, (field[:len(tepsp)] - 40), "m", label = "Field" )
    ax.plot( tepsp, pt * 0.5 - 1, "g", label = "Trigger" )
    #ax.plot( [0,0,0.25], [10,5,5], color="black", linewidth=2.5 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.get_xaxis().set_ticks([])
    #ax.get_yaxis().set_ticks([])
    #ax.set_ylabel( "epsp (mV)" )
    #ax.legend( loc = "upper right", ncol=3, frameon = False, fontsize = 14 )
    ax.legend( loc = "upper left", frameon = False, fontsize = 14 )
    #ax.set_xlim( 0, PLOTLEN )
    #ax.set_ylim( -7, 10 )
    ax.text( -0.08, 1.05, "A", fontsize = 22, weight = "bold", transform=ax.transAxes )
    #ax.set_xlabel("Time (s)")

def scanData( df ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    temp = df['stimFreq'].unique()
    freq5 = { ff:[] for ff in temp }
    pk5 = []
    for cellIdx, cell in enumerate( cellList ):
        alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
        dcell = df.loc[df["cellID"] == cell]
        freqList = dcell['stimFreq'].unique()
        for ff in freqList:
            dfreq = dcell.loc[dcell["stimFreq"] == ff]
            sweepList = dfreq['sweep'].unique()
            for ss in sweepList:
                dsweep = dfreq.loc[ dfreq['sweep'] == ss ]
                seqList = dsweep['exptSeq'].unique()
                #for seq in [seqList[0],]:
                for seq in seqList:
                    dseq = dsweep.loc[dsweep['exptSeq'] == seq]
                    idx += 1
                    epks = parseRow( dseq, cell, ff )
                    OK = (len( epks ) > 0 and epks[0] > np.mean( epks )*0.2)
                    print( "{}, cell{}, freq{}, sweep{}, seq{}, {}".format( 
                        idx, cell, ff, ss, seq, "" if OK else "bad" ) )
                    if OK:
                        norm = np.array(epks)/epks[0]
                        pk5.append( norm )
                        freq5[ff].append( norm )
                        #if np.mean( epks )> 10:
                        #    print( "MEAN = ", np.mean(epks) )

    return pk5, freq5

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
    ax2.set_xlabel( "std:EPSC" )
    ax2.set_ylabel( "#" )
    ax2.text( -0.20, 1.05, "M", fontsize = 22, weight = "bold", transform=ax2.transAxes )

def setFittingParams( cell ):
    alphaTau1 = 0.005 * SAMPLE_FREQ
    alphaTau2 = 0.042 * SAMPLE_FREQ
    alphaDelay = int( 0.010 * SAMPLE_FREQ )
    pkDelay= int( 0.020 * SAMPLE_FREQ )
    if cell == 521:
        alphaTau2 = 0.018 * SAMPLE_FREQ
        alphaTab = alphaTab / max( alphaTab )
        alphaDelay = int( 0.005 * SAMPLE_FREQ )
        pkDelay = int( 0.015 * SAMPLE_FREQ )
    elif cell == 111:
        alphaTau1 = 0.001 * SAMPLE_FREQ
        alphaTau2 = 0.005 * SAMPLE_FREQ
        alphaDelay = int( 0.005 * SAMPLE_FREQ )
        pkDelay = int( 0.006 * SAMPLE_FREQ )
    elif cell == 0:
        alphaTau1 = 0.010 * SAMPLE_FREQ
        alphaTau2 = 0.012 * SAMPLE_FREQ
        alphaDelay = int( 0.0015 * SAMPLE_FREQ )
        pkDelay = int( 0.016 * SAMPLE_FREQ )
    else:
        alphaTau1 = 0.002 * SAMPLE_FREQ
        alphaTau2 = 0.030 * SAMPLE_FREQ
        alphaDelay = int( 0.007 * SAMPLE_FREQ )
        pkDelay = int( 0.015 * SAMPLE_FREQ )

    alphaTab = np.array( [ dualAlphaFunc(t, alphaTau1, alphaTau2) for t in range( ALPHAWINDOW ) ] )
    alphaTab = alphaTab / max( alphaTab )

    return alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2

def main():
    global pulseTrig
    parser = argparse.ArgumentParser( description = "Read and display surprise protocol data from specified file" )
    parser.add_argument( "fname", type = str, help = "Required: File with data. Required to be hdf5." )
    parser.add_argument( "-c", "--cellNum", type = int, help = "Optional: Cell number to select. Default: select all cells." )
    args = parser.parse_args()

    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (10,16) )
    gs = fig.add_gridspec( 4, 2 ) # 4 rows, 2 cols
    #fig, ax = plt.subplots( nrows = 3, ncols=3, figsize = (18, 15) )

    # Set up the stimulus timings
    df = pandas.read_hdf( args.fname )
    cells = df['cellID'].unique()
    if args.cellNum:
        if args.cellNum in cells:
            df = df.loc[df['cellID']==args.cellNum]
        else:
            print( "Error: selected cell '{}' not in file. Options are: {}".format( args.cellNum, cells ) )
            quit()
    ax = fig.add_subplot( gs[0,:] )
    #celldf = df.loc[df['cellID'] == 2822]
    #celldf = df.loc[df['cellID'] == 2822]
    #celldf = df.loc[df['cellID'] == 2681]
    panelA_SampleTrace( ax, df )
    pk5, freq5 = scanData( df )
    #print( len( pk5 ) )
    #print( len( pk5[0] ) )
    #panelB_fepspVsTime( fig.add_subplot(gs[1,0]), fpk5, fpk15 )
    for idx, ff in enumerate(freq5):
        panelC_epspVsTime( fig.add_subplot(gs[idx+1,:]), freq5[ff], ff, idx)
    #panelD_fepspProbVsTime( fig.add_subplot(gs[2,0]), fpk5, fpk15 )
    #panelE_epspProbVsTime( fig.add_subplot(gs[2,1]), pk5, pk5 )
    #panelF_fepspVsISI(  fig.add_subplot(gs[3,0]), fpk5, fpk15 )
    #panelG_epspVsISI(  fig.add_subplot(gs[3,1]), pk5, pk5 )
    #panelHI_epspVsField( fig.add_subplot(gs[4,0]), fpk5, "H" )
    #panelHI_epspVsField( fig.add_subplot(gs[4,1]), fpk15, "I" )
    #panelJ_fpkHisto( fig.add_subplot(gs[5,0]), fpk5, fpk15 )
    #panelK_epkHisto( fig.add_subplot(gs[5,1]), pk5, pk5 )
    #panelLM_varianceHisto( fig.add_subplot(gs[6,0]), fig.add_subplot(gs[6,1]), patDict )
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
