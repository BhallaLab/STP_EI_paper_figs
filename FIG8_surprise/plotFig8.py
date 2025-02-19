import pandas
import pylab
import numpy as np
import math
import argparse
from scipy.stats import linregress
from scipy.stats import wilcoxon

import matplotlib.pyplot as plt
datafile = "../../../2022/VC_DATA/all_cells_surprise_CC_long.h5"
simfile = "SURP37/simData_orig_0.h5"
simfile2 = "SURP37/simData_zeroIndices_32.h5"
#simfile = "SURP_DET_34_COMPARE/surp33_gluTauDown_orig_0.h5"
#simfile2 = "SURP_DET_34_COMPARE/surp33_RMdown_orig_0.h5"
#simfile2 = "SURP33/simData_zeroIndices_224.h5"
#simfile_determ_withSTP = "SURP_DET_34_COMPARE/surp33_orig_0.h5"
simfile_determ_withSTP = "SURP49_determ/simDeterm_orig_0.h5"
simfile_determ_noSTP = "SURP49_determ/simDeterm_modelName_BothPresyn90_noSTP.g.h5"
simfile_determ_noGABA = "SURP49_determ/simDterm_wtGABA_0.h5"
simfile_determ_hi_pCA3_CA1 = "SURP49_determ/simDeterm_pCA3_CA1_0.1.h5"
#simfile2 = "SURP22/simData_orig_0.h5"
paramSweepFile = "SURP36/surprise_param_dep2.h5"
exampleCell = 2821
referenceVals = { "zeroIndices": 192, "wtGlu": 3, "wtGABA": 10,
        "pCA3_CA1": 0.02 ,"pCA3_Inter": 0.01 ,"pInter_CA1":0.01 }

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


def findPeaks( pkDelay, Vm, pulses, width = 0.008 ):
    if len( Vm ) < NUM_SAMPLES or len( pulses ) != 33:
        return []
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in pulses:
        idx1 = idx
        idx2 = idx + pkDelay - widthSamples
        #print( idx, idx2, len( Vm ) )
        val = np.quantile( Vm[idx - widthSamples:idx + widthSamples], 0.01 )
        pk = np.quantile( Vm[idx2:idx2 + widthSamples*2], 0.99 )
        pks.append( pk )
        #print( "{}  {:.3f}  {:.3f}".format( idx, val, pk ) )
    return pks

    dpks = []
    prev = prev2 = prev3 =  pks[0]
    # Pad it out.
    pks.append(pks[-1])
    for pp in pks[1:]:
        dpks.append( pp + prev - prev2 -prev3 )
        prev3 = prev2
        prev2 = prev
        prev = pp
    #print( len( pks ), len( dpks ), len( pulses ) )
    #return dpks

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
    temp = sum( pks[3:8] + pks[11:16] + pks[19:24] + pks[27:32] ) / 16.0
    pks = np.array( pks ) / temp

    #normIdx = [slice(3, 8), slice(11, 16), slice(19, 24), slice(27, 32)]
    # Calculate the mean of the specified columns for each row
    #pks /= max( pks )

    '''
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
    '''

    return pks 
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


def panelC_epspVsTime( ax, pk5, ff, label, isSim ):
    pk5 = np.array( pk5 )
    pk5[pk5>7] = np.median( pk5 )   # filter out outliers.
    mean5 = np.mean( pk5, axis = 0 )
    med5 = np.median( pk5, axis = 0 )
    print( "SHAPE =  ", pk5.shape, med5.shape )
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
    #ax.plot( PulseTrain[ff] / SAMPLE_FREQ, med5, color="red", markersize=10, label = "Median" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    if isSim == 0:
        ax.set_ylabel( "Norm EPSP" )
    if label in ["C", "D", "E"]:
        ax.set_ylim( 0, 1.5 )
    else:
        ax.set_ylim( 0, 2.1 )
    #label = chr( ord("F")+ idx + 4*isSim )
    ax.text( 0.05, 0.90, str(ff)+" Hz", fontsize = 14, transform=ax.transAxes )
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

def panelA_SampleTrace( ax, dcell, label ):
    #df = dcell.loc[(dcell['sweep'] == 0)]
    df = dcell.loc[(dcell['stimFreq'] == 50)]
    sweep = 0
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams(0)
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    PLOTLEN = 2.5
    epsp = np.array(df.iloc[sweep, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    if label in ["Fi",]:
        epsp /= 1000
    trig2 = np.array(df.iloc[sweep, SAMPLE_START+NUM_SAMPLES:SAMPLE_START+2*NUM_SAMPLES ])
    pulseTrig = np.array(df.iloc[sweep, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    field = np.array(df.iloc[sweep, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    pulses = np.arange(0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    baseline = min( np.percentile(epsp, 25 ), 0.0 )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, SAMPLE_TIME, len(epsp) )

    tepsp = tepsp[:int(PLOTLEN*SAMPLE_FREQ)]
    pt = pulseTrig[:len(tepsp)]
    pt = 0.4 * pt / (max(pt) - min(pt))
    #ax.plot( tepsp, fitepsp[:len(tepsp)], "r", label = "Fit epsp" )
    #ax.plot( tepsp, (field[:len(tepsp)] - 40), "m", label = "Field" )
    #ax.plot( [0,0,0.25], [10,5,5], color="black", linewidth=2.5 )
    if label == 'Fi':
        ax.plot( tepsp, pt - 1, color = "seagreen", label = "Trigger" )
        ax.set_ylabel( "EPSP (mV)" )
        ax.set_ylim( -1.1, 3.0 )
    else:
        ax.plot( tepsp, pt*3 - 2, color = "seagreen", label = "Trigger" )
        ax.set_ylim( -2.0, 10.5 )
    ax.plot( tepsp, epsp[:len(tepsp)], color = "blue", label = "Data" )
    if label == 'Gi':
        ax.legend( loc = "upper left", frameon = False, fontsize = 14 )
        #ax.set_ylim( -1.5, 10.0 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim( -0.1, 1.5 )
    ax.set_xlabel( "Time (s)" )
    #ax.set_xlim( 0, PLOTLEN )
    ax.text( -0.20, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    title = 'Experiment' if label == 'Fi' else 'Model'
    ax.text( 0.35, 1.05, title, fontsize = 16, transform=ax.transAxes )

def scanData( df ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    temp = df['stimFreq'].unique()
    freq5 = { ff:[] for ff in sorted(temp) }
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
                    #OK = (len( epks ) > 0 and epks[0] > np.mean( epks )*0.0)
                    OK = (len( epks ) > 0)
                    print( "{}, cell{}, freq{}, sweep{}, seq{}, {}".format( 
                        idx, cell, ff, ss, seq, "" if OK else "bad" ) )
                    if OK:
                        #norm = np.array(epks)/max(epks)
                        norm = np.array(epks)
                        pk5.append( norm )
                        freq5[ff].append( norm )

    return pk5, freq5

def panelLM_varianceHisto( ax1, ax2, patDict ):
    totfvar = []
    totevar = []
    #for pattern in [46,47,48,49,50]:
    for pattern in [52,53,54,55]:
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
    for pattern in [52, 53, 54, 55]:
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

def plotExpts( fig, gs, df, setNum ):
    # Set up the stimulus timings
    cells = df['cellID'].unique()
    if setNum == 0: # real data
        df = df.loc[df['cellID']== exampleCell]
    ax = fig.add_subplot( gs[1,setNum] )
    panelA_SampleTrace( ax, df, chr(ord("F") + setNum ) + "i" )
    if setNum:
        pk5, freq5 = scanData( df )
    else:
        pk5, freq5 = scanData( df.loc[df['numSq'] == 15] )
    label = [
        ["Fii", "Fiii", "Fiv"],
        ["Gii", "Giii", "Giv"],
        ["Hii", "Hiii", "Hiv"]
    ]
    for idx, ff in enumerate(freq5):
        print( "Plotting: ", idx, ff )
        panelC_epspVsTime( fig.add_subplot(gs[idx+2,setNum]), freq5[ff], ff, label[setNum][idx], setNum)
    
def panelSweep( dfs, fig, gs, row, paramList ):
    for idx, pp in enumerate( paramList ):
        dp = dfs.loc[dfs["param"]==pp]
        flist = dp["freq"].unique()
        #dp.info()
        #print( "PARAM = ", pp, "   FLIST = ", flist )
        ax = fig.add_subplot( gs[ row, idx ] )
        for ii, ff in enumerate( reversed( flist ) ):
            dfreq = dp.loc[dp["freq"]==ff]
            x = np.array(dfreq["val"]).flatten()
            y1 = np.array(dfreq["m2"])/np.array(dfreq["m1"])
            y2 = np.array(dfreq["m4"])/np.array(dfreq["m3"])
            y3 = np.array(dfreq["m6"])/np.array(dfreq["m5"])
            y = (y1+y2+y3)/3.0  ## Vector average
            colors = ["magenta", "seagreen", "blue"]
            ax.plot( x,y, color=colors[ii], markersize=5, label = str(ff)+ " Hz" )
        ax.scatter( [referenceVals[pp]], [1], 
            marker = '^', color = 'red', s=100 )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        xlab = "Pattern sparseness" if pp == "zeroIndices" else pp
        ax.set_xlabel( xlab )
        ax.set_ylabel( "Surprise" )
        #ax.set_ylim( 0.9, max(y) * 1.1 )
        ax.set_ylim( 0.9, 2.5 )
        if idx == 0 and row == 6:
            #ax.legend( fontsize=16, frameon=False, facecolor='none', ncol=3)
            ax.legend( fontsize=14, frameon=False, facecolor='none', 
                    loc="upper left", ncol=2)
        label = chr(ord("L") + idx + (row-6)*3 )
        ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

def plotTransition( df, ax, freq = 50, numSq = 15 ):
    chosen_cells = [3101, 2681, 2682, 2822] # What about 2821?
    df2 = df[ (df['freq']==freq) & (df['numSq']==numSq) & (df['cell'].isin(chosen_cells) ) & (df['transition'].isin( [1,2] ) ) ]
    print( "DF2 shape = ", df2.shape )
    idx = 0
    for cc in chosen_cells:
        for tt in [1,2]:    # Transitions
            sample = df2[ (df2['cell'] == cc ) & (df2['transition'] == tt ) ]
            if sample.shape[0] == 0:
                continue
            #print( "CELL NUMBER = ", cc, "  Sample shape = ", sample.shape )
            #print( "PRINT IT AS IS: ", sample['pre1'] )
            #print( "NPARRAY: ", np.array(sample['pre1']) )
            #good  print( "ILOC: ", sample.iloc[0,6] )
            #good  print( "VALUE: ", sample['pre1'].values[0] )
            #bad print( "ANOTHER NP: ", sample['pre1'][0] )
            #bad print( sample['pre1'][0] )
            #bad print( np.array(sample.loc[0,'pre1']).shape )
            #print( "###################", flush = True )
            y = np.array( [
                np.array(sample['pre1'].values[0]), 
                np.array(sample['pre2'].values[0]),
                np.array(sample['post1'].values[0]), 
                np.array(sample['post2'].values[0])
            ] )
            allMean = np.mean( y.flatten() )
            ymean = np.mean( y, axis = 1 ) + 1 - allMean
            ysem = np.std(y, axis = 1, ddof = 1) / np.sqrt( y.shape[1] )
            x = np.arange( y.shape[0] ) + 0.8 + idx/20 # A little offset for clarity
            ax.errorbar( x, ymean, yerr = ysem, fmt = 'o', 
                capsize = 5, capthick = 1, markersize = 6 )
            ax.plot(x, ymean, linestyle='-', linewidth=1.5, alpha=0.7)
            idx += 1
    ax.plot([2.5,2.5], [0.5,1.4], linestyle=':', linewidth=1.5, color="black")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Pulse number" )
    ax.set_ylabel( "EPSP ratio" )
    ax.text( -0.20, 1.10, 'I', fontsize = 22, weight = "bold", transform=ax.transAxes )

def plotTransitionResponseHisto( df, ax, numSq = 15 ):
    bin_width = 0.1
    bins = np.arange(-2, 2.5 + bin_width, bin_width)
    for ff in [8, 20, 50]:
        #trans = np.array(df['pre1'].values[:]) 
        #print( "TRANS shape = ", trans.shape )
        #print( trans )
        delta = df.loc[(df['freq'] == ff) & (df['numSq']==numSq)]['delta']
        print( "DELTA shape = ", delta.shape )
        #print( "DELTA values shape = ", delta.values.shape )
        #print( "Flattened = ", delta.values.flatten().shape )
        data = []
        for dd in delta:
            #print( "DD SHAPE = ", dd.shape )
            #print( "DD = ", dd )
            data.append( dd )
        ndata = np.concatenate( data )
        print( "NDATA shape = ", ndata.shape )
        #print( "NDATA = ", ndata )
        '''
        trans8 = np.array(df.loc[(df['freq'] == ff) & (df['numSq']==numSq)]['pre1'])
        print ("TRANS8 = ", trans8 )
        trans20 = np.array(df.loc[(df['freq'] == ff) & (df['numSq']==numSq)]['pre2'].values)
        print ("TRANS20 = ", trans20 )
        trans50 = np.array(df.loc[(df['freq'] == ff) & (df['numSq']==numSq)]['post1'])
        print ("TRANS50 = ", trans50 )
        '''
        ax.hist(ndata, bins=bins, density=False, histtype='step', linewidth=2, label=str(ff) + " Hz")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( r'$\Delta$ EPSP (normalized)' )
    ax.set_ylabel( "Frequency" )
    ax.legend( loc = "upper left", frameon = False, fontsize = 14 )
    ax.text( -0.20, 1.10, 'J', fontsize = 22, weight = "bold", transform=ax.transAxes )

def plotTransitionHeatmap( df, fig, ax ):
    data = []
    for numSq in [5, 15]:
        row = []
        for ff in [8, 20, 50]:
            numSig = 0
            pdf = df.loc[(df['numSq'] == numSq) & (df['freq'] == ff)]
            sigs = ( pdf['pval'] < 0.05 ) * 1
            #print( "SIGS:   ", numSq, ff, sigs )
            numSig = sum( pdf['pval'] < 0.05 )
            print( "HEAT = ", numSq, ff, numSig, len(pdf), pdf['pval'].shape )
            row.append( 100*numSig / len( pdf ) )
        data.append( row )
    data = np.array(data)

    x_labels = [8, 20, 50]
    y_labels = ["5 Sq", "15 Sq"]

    # Create the heatmap within ax
    cax = ax.imshow(data, cmap='viridis', aspect='auto')
    # Add color bar
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label("% selective")  # Optional label for color bar

    # Add text annotations (percentage values)
    for i in range(data.shape[0]):  # Iterate over rows (numSq)
        for j in range(data.shape[1]):  # Iterate over columns (freq)
            val = int( round( data[i,j] ) )
            color = "black" if (i == 1 and j in [1,2]) else "white"
            ax.text(j, i, f"{val}%", ha='center', va='center', color=color, fontsize=14, fontweight='bold')

    # Set x and y ticks
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_yticks(np.arange(len(y_labels)))
    
    # Set x and y tick labels
    ax.set_xticklabels(x_labels)
    ax.set_yticklabels(y_labels)

    # Set axis labels
    ax.set_xlabel("Frequency (Hz)")
    #ax.set_ylabel("Y Axis Labels")
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text( -0.25, 1.10, 'K', fontsize = 22, weight = "bold", transform=ax.transAxes )

def obtainTransitions( df ):
    cells = df['cellID'].unique()
    data = []
    for cc in cells:
        for numSq in [5, 15]:
            ndf = df.loc[(df['cellID'] == cc) & ( df['numSq'] == numSq)]
            pks, freqs = scanData( ndf )
            #print( "SHAPE pks = ", np.array( pks ).shape )
            flist = ndf['stimFreq'].unique()
            for ff in flist:
                for tt in [0, 1, 2]:
                    pre1 = np.array(freqs[ff])[:,7+tt*8].flatten()
                    pre2 = np.array(freqs[ff])[:,8+tt*8].flatten()
                    post1 = np.array(freqs[ff])[:,9+tt*8].flatten()
                    post2 = np.array(freqs[ff])[:,10+tt*8].flatten()
                    delta = (2*(post1 + post2) - (pre1 + pre2)) / (post1+post2+pre1+pre2)
                    pre = np.array(freqs[ff])[:,7+tt*8:9+tt*8].flatten()
                    post = np.array(freqs[ff])[:,9+tt*8:11+tt*8].flatten()
                    wp = wilcoxon( pre, post, alternative="less" ).pvalue
                    data.append( [cc, numSq, ff, tt, wp ,pre1, pre2, post1, post2, delta ] )
    headers = ["cell", "numSq", "freq", "transition", "pval", "pre1", "pre2", "post1", "post2", "delta"]
    print( "LENS = ", len( data[0] ), len( pre ), len( post ), len( headers) )
    df = pandas.DataFrame(data, columns=headers)
    print( "transition data frame shape = ", df.shape )
    return df


def main():
    global pulseTrig
    plt.rcParams.update( {"font.size": 16} )
    fig = plt.figure( figsize = (15,21) )
    fig.suptitle( "Fig8_v21", fontsize = 16 )
    gs = fig.add_gridspec( 8, 3 ) # 9 rows, 3 cols, but top row is empty
    # Row 1: Leave blank, for the schematics.
    # Row 2: determ runs.
    dfdet2 = pandas.read_hdf( simfile_determ_noSTP )
    pk5, freq5 = scanData( dfdet2 )
    panelC_epspVsTime( fig.add_subplot(gs[0,0]), freq5[50], 50, "C", 1)

    dfdet = pandas.read_hdf( simfile_determ_withSTP )
    pk5, freq5 = scanData( dfdet )
    panelC_epspVsTime( fig.add_subplot(gs[0,1]), freq5[50], 50, "D", 1)

    #dfdet = pandas.read_hdf( simfile_determ_noGABA )
    dfdet = pandas.read_hdf( simfile_determ_hi_pCA3_CA1 )
    pk5, freq5 = scanData( dfdet )
    panelC_epspVsTime( fig.add_subplot(gs[0,2]), freq5[50], 50, "E", 1)


    # Row 3 to 6 inclusive: Waveforms and scatter plots.
    df = pandas.read_hdf( datafile )
    plotExpts( fig, gs, df, 0 )
    df2 = pandas.read_hdf( simfile )
    plotExpts( fig, gs, df2, 1 )
    df3 = pandas.read_hdf( simfile2 )
    plotExpts( fig, gs, df3, 2 )

    # Row 7: Transition stats.
    tdf = obtainTransitions( df )
    plotTransition( tdf, ax=fig.add_subplot( gs[ 5, 0 ] ) )
    plotTransitionResponseHisto( tdf, ax=fig.add_subplot( gs[ 5, 1 ] ) )
    plotTransitionHeatmap( tdf, fig, ax=fig.add_subplot( gs[ 5, 2 ] ) )

    # Row 8: Param variations: zeroIndices, wtGlu, wtGABA
    dfs = pandas.read_hdf( paramSweepFile )
    #dfs.info()
    #print( "ORIG PARAMS = ", dfs['param'].unique() )
    panelSweep( dfs, fig, gs, 6, ["zeroIndices", "wtGlu", "wtGABA"] )
    #print( "ORIG PARAMS = ", dfs['param'].unique() )
    # Row 9: Param variations: pCA3_CA1, pCA3_Inter, pInter_CA1
    panelSweep( dfs, fig, gs, 7, ["pCA3_CA1","pCA3_Inter","pInter_CA1"])

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
