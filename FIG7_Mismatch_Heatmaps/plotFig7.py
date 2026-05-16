import pandas
import pylab
import numpy as np
import math
import argparse
from scipy.stats import linregress
from scipy.stats import wilcoxon
import scipy.signal as signal

import matplotlib.pyplot as plt
datafile = "../../2022/VC_DATA/all_cells_surprise_CC_long.h5"
simfile_37deg_runs = "fig7panelC_orig_0.h5"
simfile_ei = "fig7panelAB_orig_0.h5"
exampleCell = 2821
referenceVals = { "zeroIndices": 192, "wtGlu": 3, "wtGABA": 10,
        "pCA3_CA1": 0.02 ,"pCA3_Inter": 0.01 ,"pInter_CA1":0.01 }
jitterDir = "./JITTER"

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
SAMPLE_RATIO = 10 # Ratio between SAMPLE_FREQ and SIM_SAMPLE_FREQ
chemDt = 0.0005
SAMPLE_TIME = 5
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
EPSPTHRESH = 0.0004 # Threshold for a distinct EPSP pk. In Volts
ALPHAWINDOW = int( 0.32 * SAMPLE_FREQ )
MINIMUM_PULSE_THRESH = 4e-4
PulseTrain = {}

noiseAmpl = 2.0
noiseBandwidth = 500

def addGaussianNoise( data, noiseAmpl, bandwidth, samplingFreq ):
    # Generate white Gaussian noise
    noise = np.random.randn(len(data))
    # Design a low-pass filter
    nyquistFreq = 0.5 * samplingFreq
    cutoffFreq = bandwidth / nyquistFreq
    b, a = signal.butter(4, cutoffFreq, 'lowpass')
    # Filter the noise
    filtered_noise = signal.filtfilt(b, a, noise)
    # Scale the noise
    scaled_noise = filtered_noise * noiseAmpl
    # Add noise to the signal
    noisy_signal = data + scaled_noise

    return noisy_signal

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

def findSimPeaks( pkDelay, Vm, pulses, width = 0.016 ):
    NUM_SIM_SAMPLES = NUM_SAMPLES//10
    SIM_SAMPLE_FREQ = SAMPLE_FREQ//10
    if len( Vm ) < NUM_SIM_SAMPLES or len( pulses ) != 33:
        return []
    widthSamples = int( np.round( width * SIM_SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in pulses:
        idx1 = idx
        idx2 = idx + pkDelay - widthSamples
        #print( idx, idx2, len( Vm ) )
        pk = np.quantile( Vm[idx2:idx2 + widthSamples*2], 0.95 )
        val = np.quantile( Vm[idx2 - widthSamples*2: idx2], 0.05 )
        #val = np.quantile( Vm[idx - widthSamples:idx + widthSamples], 0.01 )
        #pk = np.quantile( Vm[idx2:idx2 + widthSamples*2], 0.99 )
        pks.append( pk )
        #print( "{}  {:.5f}  {:.5f}".format( idx, val, pk ) )
    return pks

def findPeaksEI( pkDelay, Vm, pulses, width = 0.008 ):
    if len( Vm ) < NUM_SAMPLES or len( pulses ) != 33:
        return []
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    pks = []
    for idx in pulses:
        idx1 = idx
        idx2 = idx + pkDelay - widthSamples
        #print( idx, idx2, len( Vm ) )
        #val = np.quantile( Vm[idx - widthSamples:idx + widthSamples], 0.01 )
        pk = np.max( Vm[idx2:idx2 + widthSamples*2] )
        pks.append( pk*1000 )
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
    return pks 

def parseSimRow( df, cell, ff, rowCache = None ):
    # Like parseRow but with the compact sim form with SAMPLE_FREQ=2000
    global PulseTrain
    NUM_SIM_SAMPLES = NUM_SAMPLES//SAMPLE_RATIO
    SIM_SAMPLE_FREQ = SAMPLE_FREQ//SAMPLE_RATIO
    cacheKey = (df.index[0], cell, ff)
    if rowCache is not None and cacheKey in rowCache:
        epsp_raw, pulseIndices, pkDelay = rowCache[cacheKey]
    else:
        alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
        pkDelay = pkDelay // SAMPLE_RATIO
        epsp_raw = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SIM_SAMPLES ])
        edges = np.array(df.iloc[0, SAMPLE_START+NUM_SIM_SAMPLES: ] )
        assert( len( edges ) == 33 )
        pulseIndices = np.array( SIM_SAMPLE_FREQ * edges, dtype = int)
        if rowCache is not None:
            rowCache[cacheKey] = (epsp_raw, pulseIndices, pkDelay)
    PulseTrain[ff] = pulseIndices
    epsp = addGaussianNoise( epsp_raw.copy(), noiseAmpl=noiseAmpl, bandwidth=noiseBandwidth, samplingFreq = SIM_SAMPLE_FREQ )
    baseline = min(epsp )
    epsp -= baseline # hack to handle traces with large ipsps.
    pks = findSimPeaks( pkDelay, epsp, PulseTrain[ff] )
    temp = sum( pks[3:8] + pks[11:16] + pks[19:24] + pks[27:32] ) / 16.0
    pks = np.array( pks ) / temp
    return pks

def parseRowEI( df, cell, ff ):
    global PulseTrain
    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    iGlu = np.array(df.iloc[0, SAMPLE_START+1*NUM_SAMPLES:SAMPLE_START+2*NUM_SAMPLES ] )
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    iGABA = np.array(df.iloc[0, SAMPLE_START+3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )
    tepsp = np.linspace( 0, 11, len(iGlu) )
    padt = np.pad( pulseTrig, 1 )
    edges = ((padt[:-1]< 0.02) & (padt[1:]>0.02) )
    if sum( edges ) != 33:
        return []
    PulseTrain[ff] = np.arange( 0, NUM_SAMPLES, 1, dtype = int )[edges[:NUM_SAMPLES]]
    epks = findPeaksEI( pkDelay, iGlu, PulseTrain[ff] )
    #print( "EPKS = ", epks )
    ipks = findPeaksEI( pkDelay, -iGABA, PulseTrain[ff] )
    return epks, ipks 

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

def panelC_extra( ax, pk5, ff, ymin = -8, ymax = 18, doLegend = False ):
    pk5 = np.array( pk5 )
    mean5 = -np.mean( pk5, axis = 0 )
    med5 = np.median( pk5, axis = 0 )
    #print( "IPSP SHAPE =  ", pk5.shape, med5.shape )
    ax.plot( PulseTrain[ff] / SAMPLE_FREQ, mean5, color="red", markersize=10, label = "IPSC" )
    if doLegend:
        ax.legend( loc='center left', frameon=False, ncol=2 )
    ax.set_ylim( ymin, ymax )

def panelC_epspVsTime( ax, pk5, ff, label, ylabel, ymax = 1.8, hideTime = False, isCompactSim = False ):
    pk5 = np.array( pk5 )
    mean5 = np.mean( pk5, axis = 0 )
    med5 = np.median( pk5, axis = 0 )
    #print( "SHAPE =  ", pk5.shape, med5.shape )
    SF = SAMPLE_FREQ // SAMPLE_RATIO if isCompactSim else SAMPLE_FREQ
    t0 = PulseTrain[ff][:9]/SF
    t1 = PulseTrain[ff][9:17]/SF
    t2 = PulseTrain[ff][17:25]/SF
    t3 = PulseTrain[ff][25:]/SF
    wid = 0.5/ff
    for pp in pk5:
        ax.scatter( t0 + np.random.rand( len( t0 ) ) * wid, pp[:9], color="cyan", s=1 )
        ax.scatter( t1 + np.random.rand( len( t1 ) ) * wid, pp[9:17], color="yellow", s=1 )
        ax.scatter( t2 + np.random.rand( len( t2 ) ) * wid, pp[17:25], color="pink", s=1 )
        ax.scatter( t3 + np.random.rand( len( t3 ) ) * wid, pp[25:], color="palegreen", s=1 )
    #ax.scatter( PulseTrain[ff / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq" )
    #ax.scatter( PulseTrain[ff] / SAMPLE_FREQ, mean5, color="blue", s=10, label = "Means" )
    ax.plot( PulseTrain[ff] / SF, med5, color="blue", markersize=10, label = "EPSC" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if hideTime:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( ylabel )
    ax.set_ylim( 0, ymax )
    transitions = [ (PulseTrain[ff][ii]+PulseTrain[ff][ii+1])/(2*SF) for ii in [8,16,24]]
    h = -11 if label[0] == 'A' else 0.05
    ax.scatter( transitions,[h,h,h], marker = '^', color='red' )
    ax.text( 0.05, 0.95, str(ff)+" Hz", fontsize = 14, transform=ax.transAxes )
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

def scanData( df, isCompactSim = False, rowCache = None ):
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
                    if isCompactSim:
                        epks = parseSimRow( dseq, cell, ff, rowCache = rowCache )
                    else:
                        epks = parseRow( dseq, cell, ff )
                    #OK = (len( epks ) > 0 and epks[0] > np.mean( epks )*0.0)
                    OK = (len( epks ) > 0)
                    #print( "{}, cell{}, freq{}, sweep{}, seq{}, {}".format( idx, cell, ff, ss, seq, "" if OK else "bad" ) )
                    if OK:
                        #norm = np.array(epks)/max(epks)
                        norm = np.array(epks)
                        pk5.append( norm )
                        freq5[ff].append( norm )

    return pk5, freq5

def scanDataEI( df ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    temp = df['stimFreq'].unique()
    freq5 = { ff:[] for ff in sorted(temp) }
    efreq = { ff:[] for ff in sorted(temp) }
    ifreq = { ff:[] for ff in sorted(temp) }
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
                    pkGlu, pkGABA = parseRowEI( dseq, cell, ff )
                    #OK = (len( epks ) > 0 and epks[0] > np.mean( epks )*0.0)
                    OK = (len( epks ) > 0)
                    print( "EI{}, cell{}, freq{}, sweep{}, seq{}, {}".format( 
                        idx, cell, ff, ss, seq, "" if OK else "bad" ) )
                    if OK:
                        #norm = np.array(epks)/max(epks)
                        norm = np.array(epks)
                        pk5.append( norm )
                        freq5[ff].append( norm )
                        efreq[ff].append( pkGlu )
                        ifreq[ff].append( pkGABA )

    return pk5, freq5, efreq, ifreq

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
    #print( "setFittingParams for cell ", cell )
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
    elif cell <= 6:
        alphaTau1 = 0.010 * SAMPLE_FREQ
        alphaTau2 = 0.012 * SAMPLE_FREQ
        alphaDelay = int( 0.0015 * SAMPLE_FREQ )
        pkDelay = int( 0.018 * SAMPLE_FREQ )
    else:
        alphaTau1 = 0.002 * SAMPLE_FREQ
        alphaTau2 = 0.030 * SAMPLE_FREQ
        alphaDelay = int( 0.007 * SAMPLE_FREQ )
        pkDelay = int( 0.015 * SAMPLE_FREQ )

    alphaTab = np.array( [ dualAlphaFunc(t, alphaTau1, alphaTau2) for t in range( ALPHAWINDOW ) ] )
    alphaTab = alphaTab / max( alphaTab )

    return alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2


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

def plotTransitionResponseHisto( df, ax, panel, numSq = 15 ):
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

def plotTransitionHeatmap( tdfList, fig, panel, text, ax ):
    data = []
    sdevData = []
    for numSq in [5, 15]:
        row = []
        sdevRow = []
        for ff in [8, 20, 50]:
            numSigVals = []
            for tdf in tdfList:
                pdf = tdf.loc[(tdf['numSq'] == numSq) & (tdf['freq'] == ff)]
                numSigVals.append( sum( pdf['pval'] < 0.05 ) )
            print( "HEAT = ", numSq, ff, np.mean(numSigVals), len(pdf), len(tdfList), "iters" )
            pctVals = 100 * np.array( numSigVals ) / len( pdf )
            row.append( np.mean( pctVals ) )
            sdevRow.append( np.std( pctVals ) )
        data.append( row )
        sdevData.append( sdevRow )
    data = np.array(data)
    sdevData = np.array(sdevData)

    x_labels = [8, 20, 50]
    y_labels = ["5 Sq", "15 Sq"]

    # Create the heatmap within ax
    #vmax = max( data )
    cax = ax.imshow(data, cmap='viridis', aspect='auto', vmin = 0, vmax = 40)
    # Add color bar
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label("% selective")  # Optional label for color bar

    # Add text annotations (percentage values)
    for i in range(data.shape[0]):  # Iterate over rows (numSq)
        for j in range(data.shape[1]):  # Iterate over columns (freq)
            val = int( round( data[i,j] ) )
            #color = "black" if (i == 1 and j in [1,2]) else "white"
            color = "black" if val>25 else "white"
            if len( tdfList ) > 1:
                sdev = int( round( sdevData[i,j] ) )
                ax.text(j, i, f"{val}±{sdev}%", ha='center', va='center', color=color, fontsize=12, fontweight='bold')
            else:
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
    ax.text( 0.25, 1.10, text, fontsize = 14, transform=ax.transAxes )
    ax.text( -0.25, 1.10, panel, fontsize = 22, weight = "bold", transform=ax.transAxes )

def obtainTransitions( df, isCompactSim = False, nIter = 16 ):
    doPrint = False
    if isCompactSim:
        np.random.seed( 1234 )
    cells = df['cellID'].unique()
    prepostidx = {8:[7,9,9,11], 20:[6,9,9,12], 50:[6,9,10,13]}
    headers = ["cell", "numSq", "freq", "transition", "pval", "pre1", "pre2", "post1", "post2", "delta"]
    # Pre-filter ndf for each (cell, numSq) pair and build the row cache on the
    # first iteration so DataFrame extraction only happens once across all iterations.
    ndfCache = { (cc, numSq): df.loc[(df['cellID'] == cc) & (df['numSq'] == numSq)]
                 for cc in cells for numSq in [5, 15] }
    rowCache = {} if isCompactSim else None
    allDfs = []
    for iterIdx in range( nIter ):
        data = []
        for cc in cells:
            for numSq in [5, 15]:
                ndf = ndfCache[(cc, numSq)]
                pks, freqs = scanData( ndf, isCompactSim, rowCache = rowCache )
                flist = ndf['stimFreq'].unique()
                for ff in flist:
                    for tt in [0, 1, 2]:
                        pre1 = np.array(freqs[ff])[:,7+tt*8].flatten()
                        pre2 = np.array(freqs[ff])[:,8+tt*8].flatten()
                        post1 = np.array(freqs[ff])[:,9+tt*8].flatten()
                        post2 = np.array(freqs[ff])[:,10+tt*8].flatten()
                        delta = (2*(post1 + post2) - (pre1 + pre2)) / (post1+post2+pre1+pre2)

                        p,q,r,s = prepostidx[ff]
                        pre = np.array(freqs[ff])[:,p+tt*8:q+tt*8].mean(axis=1)
                        post =np.array(freqs[ff])[:,r+tt*8:s+tt*8].mean(axis=1)
                        if (pre == post).all():
                            print( f"\nOOPS, all equal: cell {cc}.{numSq}.{ff}.{tt}: {pre}\n" )
                            wp = 1
                        else:
                            wp = wilcoxon( pre, post, alternative="less" ).pvalue
                        data.append( [cc, numSq, ff, tt, wp ,pre1, pre2, post1, post2, delta ] )
                        if doPrint and (wp < 0.05):
                            print( f"cell {cc}.{numSq}.{ff}.{tt}: p = {wp}, means = {np.mean(pre)},{np.mean(post)}" )
                            for qq, (ii, jj) in enumerate(zip ( pre, post )):
                                print( "{:-4d} {:.4f}, {:.4f}, {}".format(
                                    qq, ii, jj, ii < jj ) )
                            print( "     MeanPre MeanPost SumRising")
                            print( "      {:.4f}, {:.4f}, {}".format(
                                np.mean(pre), np.mean(post),
                                sum( pre < post ) ) )
                            print( "---------------------------------------" )
        if iterIdx == 0:
            print( "LENS = ", len( data[0] ), len( pre ), len( post ), len( headers) )
        iterDf = pandas.DataFrame(data, columns=headers)
        if iterIdx == 0:
            print( "transition data frame shape = ", iterDf.shape )
        allDfs.append( iterDf )
    return allDfs


def main():
    global pulseTrig, noiseAmpl
    parser = argparse.ArgumentParser()
    parser.add_argument( "--noiseAmpl", type=float, default=noiseAmpl, help="Amplitude of Gaussian noise added to sim traces (default: %(default)s)" )
    args = parser.parse_args()
    noiseAmpl = args.noiseAmpl
    plt.rcParams.update( {"font.size": 16} )
    fig = plt.figure( figsize = (15,15) )
    #fig.suptitle( "Fig8_v21", fontsize = 16 )
    gs = fig.add_gridspec( 5, 3 ) # 5 rows, 3 cols
    dfdet = pandas.read_hdf( simfile_ei )
    pk5, freq5, efreq, ifreq = scanDataEI( dfdet )
    ax = fig.add_subplot(gs[0,0])
    panelC_epspVsTime( ax, efreq[8], 8, "Ai", "E/IPSC (pA)", hideTime = True )
    panelC_extra( ax, ifreq[8], 8, ymin=-12, ymax=19, doLegend = True )

    ax = fig.add_subplot(gs[0,1])
    panelC_epspVsTime( ax, efreq[20], 20, "Aii", "E/IPSC (pA)", hideTime = True )
    panelC_extra( ax, ifreq[20], 20, ymin=-12, ymax=19 )

    ax = fig.add_subplot(gs[0,2])
    panelC_epspVsTime( ax, efreq[50], 50, "Aiii", "E/IPSC (pA)", hideTime=True )
    panelC_extra( ax, ifreq[50], 50, ymin=-12, ymax=19 )

    pk5, freq5 = scanData( dfdet )
    panelC_epspVsTime( fig.add_subplot(gs[1,0]), freq5[8], 8, "Bi", "Norm EPSP", hideTime = True)
    panelC_epspVsTime( fig.add_subplot(gs[1,1]), freq5[20], 20, "Bii", "Norm EPSP", hideTime = True)
    panelC_epspVsTime( fig.add_subplot(gs[1,2]), freq5[50], 50, "Biii", "Norm EPSP", hideTime = True)

    dfdet = pandas.read_hdf( simfile_37deg_runs )
    pk5, freq5 = scanData( dfdet )
    panelC_epspVsTime( fig.add_subplot(gs[2,0]), freq5[8], 8, "Ci", "37°C norm EPSP")
    panelC_epspVsTime( fig.add_subplot(gs[2,1]), freq5[20], 20, "Cii", "37°C norm EPSP")
    panelC_epspVsTime( fig.add_subplot(gs[2,2]), freq5[50], 50, "Ciii", "37°C norm EPSP")

    print( "READING real data" )

    df = pandas.read_hdf( datafile )
    exampleDf = df[df['cellID'] == exampleCell]
    pk5, freq5 = scanData( exampleDf )
    panelC_epspVsTime( fig.add_subplot(gs[3,0]), freq5[8], 8, "Di", "Expt norm EPSP")
    panelC_epspVsTime( fig.add_subplot(gs[3,1]), freq5[20], 20, "Dii", "Expt norm EPSP")
    ax = fig.add_subplot(gs[3,2])
    panelC_epspVsTime( ax, freq5[50], 50, "Diii", "Expt norm EPSP")
    ax.text( PulseTrain[50][24]/SAMPLE_FREQ, 1.50, "*", fontsize = 18 )

    print( "Plotting heatmap for real data" )
    tdf = obtainTransitions( df, isCompactSim = False, nIter = 1 )
    cell50df = tdf[0][(tdf[0]["cell"] == exampleCell) & (tdf[0]['freq']==50)]
    print( cell50df[["numSq","freq","transition","delta","pval"]])
    plotTransitionHeatmap( tdf, fig, 'Ei', 'Experiment', ax=fig.add_subplot( gs[4,0] ) )

    print( "Plotting heatmap for sim data, 0 ms jitter" )
    jdfs = [ pandas.read_hdf( f"{jitterDir}/sim{ii}_jitterMean_0.0.h5") for ii in range(1, 7 )]
    combinedJdf = pandas.concat(jdfs, ignore_index=True)
    tjdf = obtainTransitions( combinedJdf, isCompactSim = True )
    plotTransitionHeatmap( tjdf, fig, 'Eii', 'Sim: Zero jitter', ax=fig.add_subplot( gs[4,1] ) )

    print( "Plotting heatmap for sim data, jitter = 6ms" )
    jdfs = [ pandas.read_hdf( f"{jitterDir}/sim{ii}_jitterMean_0.006.h5") for ii in range(1, 7 )]
    combinedJdf = pandas.concat(jdfs, ignore_index=True)
    tjdf = obtainTransitions( combinedJdf, isCompactSim = True )
    plotTransitionHeatmap( tjdf, fig, 'Eiii', 'Sim: Jitter = 6 ms', ax=fig.add_subplot( gs[4,2] ) )


    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
