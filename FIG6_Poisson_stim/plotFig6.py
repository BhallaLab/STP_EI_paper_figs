import pandas
import pylab
import numpy as np
import argparse
import math
import scipy.stats as stats
import scipy.signal as signal
import scipy.optimize as optimize
import statsmodels.api as sm

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
repeatPatterns = False
inputs = []
stimList = []
pulseTrig = []
exptData = "../../2022/VC_DATA/all_cells_SpikeTrain_CC_long.h5"
simData = "fig6poisson_orig_0.h5"
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 11
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
EPSPTHRESH = 0.4 # Threshold for a distinct EPSP pk. In mV
ALPHAWINDOW = int( 0.32 * SAMPLE_FREQ )
MINIMUM_PULSE_THRESH = 4e-4
NUMPULSE = 32   # Number of initial pulses used for slope analysis

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


def normalizePks( pk ):
    """Normalize each trial's EPSP peaks to that trial's mean. Trials with
    zero mean (no responses) are dropped."""
    return [ np.array(row) / np.mean(row) for row in pk if np.mean(row) > 0 ]


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
    return np.array(data - filtered_data)

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

def parseRow( df, cell, args ):
    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    epsp = addGaussianNoise( epsp, args.noise*1e-3, args.noiseFreq, 20000 )
    if cell == 0:
        epsp *= 1000     # Scale to mV.
    baseline = min( np.percentile(epsp, 25 ), np.mean( epsp[100:1000 ]) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pks = findPeaks( pkDelay, epsp, threshold=args.threshold*1e-3 )

    if cell == 4041:
        # 1. Get the raw data slice from the DataFrame.
        field_data_slice = df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES]

        # 2. Convert the slice to a numeric type.
        #    'errors='coerce'' turns non-numeric values (like '') into NaN.
        #    '.fillna(0)' then replaces those NaN values with 0.
        #    '.values' gets the clean NumPy array.
        field = pandas.to_numeric(field_data_slice, errors='coerce').fillna(0).values

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

def panelD_probVsTime( ax, column, pk5, pk15, perCell=None ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    pk5 = (pk5 > EPSPTHRESH)
    prob5 = np.sum( pk5, axis = 0 ) * 100 / len( pk5 )
    pk15 = (pk15 > EPSPTHRESH)
    prob15 = np.sum( pk15, axis = 0 ) * 100 / len( pk15 )
    if perCell is not None:
        for cell, data in perCell.items():
            if len(data['pk5']) > 1:
                cpk5 = (np.array(data['pk5']) > EPSPTHRESH)
                cp5 = np.sum(cpk5, axis=0) * 100 / len(cpk5)
                ax.plot(PulseTrain/SAMPLE_FREQ, cp5, color='blue', alpha=0.2, linewidth=0.8, zorder=1)
            if len(data['pk15']) > 1:
                cpk15 = (np.array(data['pk15']) > EPSPTHRESH)
                cp15 = np.sum(cpk15, axis=0) * 100 / len(cpk15)
                ax.plot(PulseTrain/SAMPLE_FREQ, cp15, color='orange', alpha=0.2, linewidth=0.8, zorder=1)
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob5, color="blue", s = 10, label = "5 Sq", zorder=3 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, prob15, color="orange", s = 10, label = "15 Sq", zorder=3 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "Probability (%)" )
    ax.set_ylim( -5, 105 )
    label = "D" + ["i", "ii"][column]
    ax.text( -0.20, 1.10, label, fontsize = 22, weight="bold", transform=ax.transAxes )
    return prob5, prob15

def panelE_epspVsTime( ax, column, pk5, pk15, perCell=None ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    mean5 = np.mean( pk5, axis = 0 )
    mean15 = np.mean( pk15, axis = 0 )
    if perCell is not None:
        for cell, data in perCell.items():
            if len(data['pk5']) > 1:
                cm5 = np.mean(np.array(data['pk5']), axis=0)
                ax.plot(PulseTrain/SAMPLE_FREQ, cm5, color='blue', alpha=0.2, linewidth=0.8, zorder=1)
            if len(data['pk15']) > 1:
                cm15 = np.mean(np.array(data['pk15']), axis=0)
                ax.plot(PulseTrain/SAMPLE_FREQ, cm15, color='orange', alpha=0.2, linewidth=0.8, zorder=1)
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq", zorder=3 )
    ax.scatter( PulseTrain / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq", zorder=3 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim( top=5.0 )
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    label = "E" + ["i", "ii"][column]
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return mean5, mean15

def panelF_epspVsISI( ax, column, pk5, pk15, perCell=None ):
    pk5 = np.array( pk5 )
    pk15 = np.array( pk15 )
    mean5 = np.mean( pk5, axis = 0 )
    mean15 = np.mean( pk15, axis = 0 )
    padt = np.pad( PulseTrain, 1)
    isi = PulseTrain - padt[:len( PulseTrain )]
    if perCell is not None:
        for cell, data in perCell.items():
            if len(data['pk5']) > 1:
                cm5 = np.mean(np.array(data['pk5']), axis=0)
                ax.scatter(isi/SAMPLE_FREQ, cm5, color='blue', alpha=0.2, s=4, zorder=1)
            if len(data['pk15']) > 1:
                cm15 = np.mean(np.array(data['pk15']), axis=0)
                ax.scatter(isi/SAMPLE_FREQ, cm15, color='orange', alpha=0.2, s=4, zorder=1)
    ax.scatter( isi / SAMPLE_FREQ, mean5, color="blue", s=10, label = "5 Sq", zorder=3 )
    ax.scatter( isi / SAMPLE_FREQ, mean15, color="orange", s=10, label = "15 Sq", zorder=3 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim( top=5.0 )
    ax.set_xlabel( "ISI (s)" )
    ax.set_ylabel( "EPSP (mV )" )
    ax.set_xlim( -0.01, 0.2 )
    label = "F" + ["i", "ii"][column]
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    return mean5, mean15, isi / SAMPLE_FREQ

def panelG_epkHisto( ax, column, pk5, pk15, perCell=None ):
    pk5 = np.array( pk5 ).flatten()
    pk5 = pk5[pk5 > EPSPTHRESH]
    pk15 = np.array( pk15 ).flatten()
    pk15 = pk15[pk15 > EPSPTHRESH]
    if perCell is not None:
        for cell, data in perCell.items():
            if len(data['pk5']) > 1:
                cv5 = np.array(data['pk5']).flatten()
                cv5 = cv5[cv5 > EPSPTHRESH]
                if len(cv5) > 0:
                    ax.hist(cv5, bins=20, alpha=0.15, histtype="step", linewidth=1,
                            edgecolor="blue", range=(0, 12), zorder=1)
            if len(data['pk15']) > 1:
                cv15 = np.array(data['pk15']).flatten()
                cv15 = cv15[cv15 > EPSPTHRESH]
                if len(cv15) > 0:
                    ax.hist(cv15, bins=20, alpha=0.15, histtype="step", linewidth=1,
                            edgecolor="orange", range=(0, 12), zorder=1)
    ax.hist( pk5, bins = 20, alpha = 0.5, label = "5 sq", histtype = "step", linewidth = 2, edgecolor = "blue", zorder=3 )
    ax.hist( pk15, bins = 20, alpha = 0.5, label = "15 sq", histtype = "step", linewidth = 2, edgecolor = "orange", zorder=3 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim( -0.5, 12 )
    ax.set_xlabel( "EPSP (mV )" )
    ax.set_ylabel( "#" )
    label = "G" + ["i", "ii"][column]
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelC_FFT( ax, dcell, column, args ):
    MAX_FREQ = 60
    ipat = dcell["patternList"].astype(int).unique()
    recordings5 = []   # list of (cellID, epsp_array)
    recordings15 = []
    for pp in ipat:
        df = dcell.loc[dcell["patternList"] == pp]
        print( f"FFT for pat {pp}, column {column}" )
        for i in range(len(df)):
            row = df.iloc[[i]]  # Get a single-row DataFrame
            cell = int(row['cellID'].iloc[0])
            pulseTrig = np.array(row.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
            pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
            if pulseThresh < MINIMUM_PULSE_THRESH:
                continue
            epsp = np.array(row.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
            if cell == 0:
                epsp *= 1000
            if pp < 51:
                recordings5.append( (cell, epsp) )
            else:
                recordings15.append( (cell, epsp) )

    def compute_psds( recordings ):
        """Returns (avg_psd, per_cell_avg_psd, frequencies, freq_idx)."""
        cell_psd_lists = {}
        all_psds = []
        frequencies = None
        freq_idx = None
        for cell, voltage in recordings:
            f, psd = signal.welch(voltage, fs=SAMPLE_FREQ, nperseg=65536)
            if freq_idx is None:
                freq_idx = np.where(f <= MAX_FREQ)[0][-1]
                frequencies = f
            psd_trunc = np.abs(psd[:freq_idx+1])
            all_psds.append(psd_trunc)
            cell_psd_lists.setdefault(cell, []).append(psd_trunc)
        avg_psd = np.mean(all_psds, axis=0) if all_psds else np.array([])
        per_cell = {c: np.mean(v, axis=0) for c, v in cell_psd_lists.items()}
        return avg_psd, per_cell, frequencies, freq_idx

    avg_psd5, perCell_psd5, frequencies, freq_idx = compute_psds(recordings5)
    avg_psd15, perCell_psd15, _, _ = compute_psds(recordings15)

    if frequencies is not None:
        ax.plot(frequencies[:freq_idx+1], avg_psd5, color = "blue", linewidth=2)
        ax.plot(frequencies[:freq_idx+1], avg_psd15, color = "orange", linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (V^2/Hz)')
    ax.set_ylim( top=0.7 )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    label = "C" + ['i', 'ii'][column]
    ax.text( -0.20, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    if column == 0:
        ax.text( 0.30, 1.00, "Experiment", fontsize = 20, transform=ax.transAxes )
    else:
        ax.text( 0.30, 1.00, "Simulation", fontsize = 20, transform=ax.transAxes )
    return avg_psd5.astype( float ), avg_psd15.astype( float ), perCell_psd5, perCell_psd15

def panelA_SampleTrace( ax, dcell, column, args ):
    ipat = dcell["patternList"].astype(int)
    df = dcell.loc[(ipat == 46)]
    cell = df['cellID'].unique()[0]
    alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
    longAlpha = np.zeros(NUM_SAMPLES)
    longAlpha[:ALPHAWINDOW] += alphaTab
    PLOTLEN = 6.0
    pulseTrig = np.array(df.iloc[0, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    if pulseThresh < MINIMUM_PULSE_THRESH:
        return [], [], 0
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    epsp = addGaussianNoise( epsp, args.noise*1e-3, args.noiseFreq, 20000 )
    if cell == 0:
        epsp *= 1000
    else:
        pulseTrig *= 100
    field = np.array(df.iloc[0, SAMPLE_START + 3*NUM_SAMPLES:SAMPLE_START+4*NUM_SAMPLES ] )

    baseline = min( np.percentile(epsp, 25 ), np.mean( epsp[100:1000 ]) )
    epsp -= baseline # hack to handle traces with large ipsps.
    tepsp = np.linspace( 0, 11, len(epsp) )
    pks = findPeaks( pkDelay, epsp, threshold=args.threshold*1e-3 )
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
    ax.plot( tepsp, epsp[:len(tepsp)], "b", label = "Data EPSP   " )
    ax.plot( tepsp, fitEPSP[:len(tepsp)], "r", label = "Fit EPSP" )
    ax.plot( tepsp, pt, "g", label = "Trigger" )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel( "EPSP (mV)" )
    ax.set_xlabel( "Time (s)" )
    ax.legend( loc = "upper right", bbox_to_anchor=(1.0, 1.05), ncol=3, frameon=False, fontsize = 14 )
    ax.set_xlim( 0, PLOTLEN )
    ax.set_ylim( top=10 )
    label = chr( ord("A") + column )
    ax.text( -0.10, 1.05, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
    #ax.set_xlabel("Time (s)")

def scanData( df, args ):
    idx = 0
    cellList = df['cellID'].unique()
    pk5 = []
    pk15 = []
    fpk5 = []
    fpk15 = []
    patDict = { pp:[] for pp in [46,47,48,49,50,52,53,55] } # 8 patterns.
    perCell = {}    # per-cell peak lists: {cellID: {'pk5': [...], 'pk15': [...]}}
    for cellIdx, cell in enumerate( cellList ):
        if cell == 4001:
            continue
        alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2 = setFittingParams( cell )
        dcell = df.loc[df["cellID"] == cell]
        ipat = dcell["patternList"].astype(int)
        patList = ipat.unique()
        perCell[cell] = {'pk5': [], 'pk15': []}
        for pp in patList:
            dpat = dcell.loc[ipat == pp]
            sweepList = dpat['sweep'].unique()
            print( f"Row{idx}, cell{cell}, pat{pp}", flush=True )
            for ss in sweepList:
                dsweep = dpat.loc[ dpat['sweep'] == ss ]
                seqList = dsweep['exptSeq'].unique()
                for seq in seqList:
                    dseq = dsweep.loc[dsweep['exptSeq'] == seq]
                    idx += 1
                    fpks, epks = parseRow( dseq, cell, args )
                    patDict[pp].append( [fpks, epks] )
                    if pp in [46,47,48,49,50]:
                        pk5.append(epks)
                        perCell[cell]['pk5'].append(epks)
                        if cell == 4041: # The only one with field data
                            fpk5.append( [fpks, epks] )
                    else:
                        pk15.append(epks)
                        perCell[cell]['pk15'].append(epks)
                        if cell == 4041:
                            fpk15.append( [fpks, epks] )

    return pk5, pk15, fpk5, fpk15, patDict, perCell


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
        alphaTau1 = 0.008 * SAMPLE_FREQ
        alphaTau2 = 0.016 * SAMPLE_FREQ
        alphaDelay = int( 0.002 * SAMPLE_FREQ )
        pkDelay = int( 0.008 * SAMPLE_FREQ )
    else: # Use defaults from above
        pass

    alphaTab = np.array( [ dualAlphaFunc(t, alphaTau1, alphaTau2 ) for t in range( ALPHAWINDOW ) ] )
    alphaTab = alphaTab / max( alphaTab )
    return alphaTab, alphaDelay, pkDelay, alphaTau1, alphaTau2


def doLinFit( time, Y, name, silent=False ):
    """Fits a straight line to the data and returns parameter estimates,
    error estimates, R-squared, and p-value.

    Args:
        x: The independent variable data (1D array).
        y: The dependent variable data (1D array).
        silent: If True, suppress printed output.
    """
    #Y2 = Y / Y.mean(axis=1, keepdims=True)
    Y2 = Y
    Y_all = Y2.flatten()  # Flatten the (75, 40) array to a single vector
    numSamples = Y2.shape[0]
    time_all = np.tile( time, numSamples )
    try:
        # Use statsmodels for linear regression (more stats info)
        X = sm.add_constant(time_all)  # Add a constant for the intercept
        model = sm.OLS(Y_all, X).fit()
        slope = model.params[1]
        intercept = model.params[0]
        slope_err = model.bse[1]
        intercept_err = model.bse[0]
        r_squared = model.rsquared
        resid_std = model.resid.std()
        p_value = model.pvalues[1]  # p-value for the slope
        fitted_y = model.fittedvalues

        #Confidence Intervals
        alpha = 0.05
        n = len(Y_all)
        p = 2 #Number of parameters
        dof = n - p
        t_critical = stats.t.ppf(1 - alpha/2, dof)
        slope_ci = model.params[1] + np.array([-1, 1]) * t_critical * model.bse[1]
        intercept_ci = model.params[0] + np.array([-1, 1]) * t_critical * model.bse[0]

        if not silent:
            print( "Panel E Linear Fitting of first 1.8 sec: ", name )
            print("slope: {:.3g}±{:.3g}  (95% CI: {:.3g}, {:.3g})".format(
                slope, slope_err, float(slope_ci[0]), float(slope_ci[1])) )
            print("intcpt: {:.3g}±{:.3g}  (95% CI: {:.3g}, {:.3g})".format( intercept, intercept_err, intercept_ci[0], intercept_ci[1]) )
            print("r_sq: {:.3g}  p_value: {:.3g}    resid_std: {:.3g}".format(
                r_squared, p_value, resid_std ) )
            print()
        results = {
            'slope': slope,
            'intercept': intercept,
            'slope_err': slope_err,
            'intercept_err': intercept_err,
            'slope_ci': slope_ci,
            'intercept_ci': intercept_ci,
            'r_squared': r_squared,
            'p_value': p_value,
            'fitted_y': fitted_y
        }
        return results

    except Exception as e:  # Catch any potential errors
        print(f"Error fitting line: {e}")
        return None  # or raise the exception if you prefer


def doExpFitISI( Y, name, silent=False ):
    padt = np.pad( PulseTrain, 1)
    isi = (PulseTrain - padt[:len( PulseTrain )]) / SAMPLE_FREQ
    n_samples = Y.shape[0]
    n_timepoints = len( isi )
    Y2 = Y / Y.mean(axis=1, keepdims=True)
    Y2 = Y

    # Reshape Y to a 1D array
    Y_all = Y2.flatten()  # Flatten the (75, 40) array to a single vector
    time_all = np.tile(isi, n_samples) #tile time so it matches the length of Y_all

    # Exponential decay function
    def exponential_decay(t, y0, tau, y1):
        return y0 * np.exp(-t / tau) + y1

    # Initial guesses (adjust these as needed)
    y0_guess = np.max(Y_all) - np.min(Y_all)
    tau_guess = 0.03
    y1_guess = np.min(Y_all)

    try:
        popt, pcov = optimize.curve_fit(exponential_decay, time_all, Y_all, p0=[y0_guess, tau_guess, y1_guess], maxfev=5000)
        y0, tau, y1 = popt
        perr = np.sqrt(np.diag(pcov))  # Standard deviations of the params
        alpha = 0.05  # Significance level (1 - confidence level)
        n = len(Y_all)
        p = len(popt)
        dof = n - p  # Degrees of freedom
        t_critical = stats.t.ppf(1 - alpha/2, dof) #t-critical value for 2-sided test
        ci_lower = popt - t_critical * perr
        ci_upper = popt + t_critical * perr

        if not silent:
            print( "Panel F: ISI decline Exp Fitting: ", name )
            print("y0: {:.3g} ± {:.3g}  (95% CI: {:.3g}, {:.3g})".format( y0, perr[0], ci_lower[0], ci_upper[0]))
            print("tau: {:.3g} ± {:.3g}  (95% CI: {:.3g}, {:.3g})".format(
                tau, perr[1], ci_lower[1], ci_upper[1]))
            print("y1: {:.3g} ± {:.3g}  (95% CI: {:.3g}, {:.3g})".format(
                y1, perr[2], ci_lower[2], ci_upper[2]))

        residuals = Y_all - exponential_decay(time_all, y0, tau, y1)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((Y_all - np.mean(Y_all))**2)
        r_squared = 1 - (ss_res / ss_tot)

        if not silent:
            print(f"R-squared: {r_squared}")

            if hasattr(Y2, 'shape') and len(Y2.shape) == 2:
                n_replicates = Y2.shape[0]
                n_timepoints = Y2.shape[1]
                n_total = n_replicates * n_timepoints

                ss_pure_error = 0
                for t_idx in range(n_timepoints):
                    y_at_t = Y2[:, t_idx]
                    ss_pure_error += np.sum((y_at_t - np.mean(y_at_t))**2)

                df_pure_error = n_total - n_timepoints
                df_residual = n_total - len(popt)
                ms_residual = ss_res / df_residual
                ms_pure_error = ss_pure_error / df_pure_error

                f_statistic = ms_residual / ms_pure_error
                p_value = 1 - stats.f.cdf(f_statistic, df_residual, df_pure_error)

                print(f"Lack-of-fit F-statistic: {f_statistic}")
                print(f"Lack-of-fit p-value: {p_value}\n")
            else:
                print("Lack-of-fit test requires replicated data. Skipping.")
            print( flush = True )
        return { 'y0': y0, 'y0_err': perr[0], 'tau': tau, 'tau_err': perr[1], 'y1': y1, 'y1_err': perr[2] }

    except RuntimeError as e:
        print(f"Fit failed: {e}")
        return None


def _perCellMetrics( perCell ):
    """Compute per-cell scalar summaries used by the stats functions.

    Returns a dict keyed by cellID, each containing:
      prob5, prob15   : mean response probability (t>2s), scalar
      slope5, slope15 : linear-fit slope on absolute peaks (t<=1.8s)
      tau5, tau15     : ISI exp-decay tau on absolute peaks
      mu5, mu15       : lognormal mu of absolute non-zero peaks
      sig5, sig15     : lognormal sigma of absolute non-zero peaks
    """
    t2_mask  = PulseTrain / SAMPLE_FREQ > 2.0
    t18_mask = PulseTrain / SAMPLE_FREQ <= 1.8
    pt = PulseTrain / SAMPLE_FREQ

    out = {}
    for cell, data in perCell.items():
        d = {}
        for sq, pk_key in [('5', 'pk5'), ('15', 'pk15')]:
            cpk = data[pk_key]
            if len(cpk) < 2:
                d[f'prob{sq}']  = None
                d[f'slope{sq}'] = None
                d[f'tau{sq}']   = None
                d[f'mu{sq}']    = None
                d[f'sig{sq}']   = None
                continue

            cpk_arr = np.array(cpk)

            # probability (t>2s)
            cprob = np.sum(cpk_arr > EPSPTHRESH, axis=0) * 100.0 / len(cpk_arr)
            d[f'prob{sq}'] = float(np.mean(cprob[t2_mask]))

            # slope (t<=1.8s) — absolute peaks
            lf = doLinFit(pt[:NUMPULSE], cpk_arr[:, :NUMPULSE],
                          f"cell{cell} {sq}sq", silent=True)
            d[f'slope{sq}'] = lf['slope'] if lf else None

            # ISI tau — absolute peaks
            ef = doExpFitISI(cpk_arr, f"cell{cell} {sq}sq", silent=True)
            d[f'tau{sq}'] = ef['tau'] if ef else None

            # lognormal fit on non-zero absolute peaks
            ev = cpk_arr.flatten()
            ev = ev[ev > 0]
            if len(ev) >= 5:
                try:
                    s_v, _, sc_v = stats.lognorm.fit(ev, floc=0)
                    d[f'mu{sq}']  = float(np.log(sc_v))
                    d[f'sig{sq}'] = float(s_v)
                except Exception:
                    d[f'mu{sq}']  = None
                    d[f'sig{sq}'] = None
            else:
                d[f'mu{sq}']  = None
                d[f'sig{sq}'] = None

        out[cell] = d
    return out


def _reportRange( metric_name, cell_vals, sim_val ):
    """Print per-cell values and whether sim_val falls within their range."""
    valid = [(c, v) for c, v in cell_vals if v is not None]
    if not valid:
        print(f"    {metric_name}: no valid expt cells")
        return
    vals = [v for _, v in valid]
    lo, hi = min(vals), max(vals)
    mn, sd = np.mean(vals), np.std(vals)
    cell_str = "  ".join(f"cell{c}={v:.3g}" for c, v in valid)
    print(f"    {metric_name} per cell: {cell_str}")
    print(f"    expt range [{lo:.3g}, {hi:.3g}], mean={mn:.3g}±{sd:.3g}")
    if sim_val is not None:
        within = lo <= sim_val <= hi
        print(f"    sim={sim_val:.3g};  within expt range: {'YES' if within else 'NO'}")
    print()


def statsC( edata, sdata ):
    """Panel C: for each expt cell compute cosine similarity of its PSD to the
    sim PSD; report whether the sim-vs-cell similarity is in the range of
    cell-vs-pooled-mean similarities (i.e. the sim is as 'close' as any cell)."""
    print( "Panel C: PSD — per-cell vs simulation comparison" )
    for sq, e_pc, s_psd, e_avg in [
        ("5sq",  edata.get('perCellPsd5',  {}), sdata['psd5'],  edata['psd5']),
        ("15sq", edata.get('perCellPsd15', {}), sdata['psd15'], edata['psd15'])
    ]:
        if not e_pc:
            # sim dataset — no per-cell breakdown, fall back to global metrics
            corr, pval = stats.pearsonr( e_avg, s_psd )
            print( f"  {sq}: Pearson corr={corr:.3f}, p={pval:.4e} (no per-cell data)" )
            continue

        cos_sim_to_sim  = []   # cosine(cell_psd, sim_psd)
        cos_sim_to_mean = []   # cosine(cell_psd, expt_mean_psd)  — "self" similarity
        for cell, cpsd in e_pc.items():
            n = min(len(cpsd), len(s_psd))
            if n == 0:
                continue
            c2sim  = np.dot(cpsd[:n], s_psd[:n])  / (np.linalg.norm(cpsd[:n])  * np.linalg.norm(s_psd[:n]))
            c2mean = np.dot(cpsd[:n], e_avg[:n]) / (np.linalg.norm(cpsd[:n]) * np.linalg.norm(e_avg[:n]))
            cos_sim_to_sim.append( (cell, float(c2sim)) )
            cos_sim_to_mean.append( (cell, float(c2mean)) )

        sim_to_mean_n = min(len(s_psd), len(e_avg))
        sim_cos = float(np.dot(s_psd[:sim_to_mean_n], e_avg[:sim_to_mean_n]) /
                        (np.linalg.norm(s_psd[:sim_to_mean_n]) * np.linalg.norm(e_avg[:sim_to_mean_n])))

        self_vals = [v for _, v in cos_sim_to_mean]
        print(f"  {sq}: cosine(sim PSD, expt mean PSD) = {sim_cos:.4f}")
        print(f"  {sq}: cosine(each expt cell PSD, expt mean PSD) — range of 'typical' similarity:")
        for cell, v in cos_sim_to_mean:
            print(f"    cell {cell}: {v:.4f}")
        lo, hi = min(self_vals), max(self_vals)
        print(f"  {sq}: expt cell range [{lo:.4f}, {hi:.4f}];  "
              f"sim within range: {'YES' if lo <= sim_cos <= hi else 'NO'}")
    print( "  (cosine similarity: 1 = identical shape, lower = diverging spectral profile)" )
    print( flush = True )


def statsD( edata, sdata ):
    """Panel D: mean response probability (t>2s) per expt cell vs sim."""
    t2_mask = PulseTrain / SAMPLE_FREQ > 2.0
    print( "Panel D: Response probability (t>2s) — per-cell comparison" )
    pcm = edata.get('perCellMetrics', {})
    for sq, s_prob in [("5sq", sdata['prob5']), ("15sq", sdata['prob15'])]:
        pk_key = sq[:-2]   # '5' or '15'
        cell_vals = [(c, d[f'prob{pk_key}']) for c, d in pcm.items()]
        sim_val   = float(np.mean(s_prob[t2_mask]))
        print(f"  {sq}:")
        _reportRange("mean prob t>2s (%)", cell_vals, sim_val)
    print( flush = True )


def statsE( edata, sdata ):
    """Panel E: initial-dip slope and steady-state amplitude per expt cell vs sim."""
    t2_mask = PulseTrain / SAMPLE_FREQ > 2.0
    print( "Panel E: EPSP vs time — per-cell comparison" )
    pcm = edata.get('perCellMetrics', {})
    for sq, s_linfit, s_pk in [
        ("5sq",  sdata['linfit5'],  sdata['pk5']),
        ("15sq", sdata['linfit15'], sdata['pk15'])
    ]:
        pk_key = sq[:-2]
        print(f"  {sq} — initial dip slope (t<=1.8s, absolute):")
        cell_slopes = [(c, d[f'slope{pk_key}']) for c, d in pcm.items()]
        sim_slope   = s_linfit['slope'] if s_linfit else None
        _reportRange("slope (/s)", cell_slopes, sim_slope)

        print(f"  {sq} — mean absolute amplitude (t>2s):")
        cell_amps = []
        for c, d in pcm.items():
            cpk = edata['perCell'][c][f'pk{pk_key}']
            if len(cpk) < 2:
                cell_amps.append((c, None))
                continue
            cell_amps.append((c, float(np.mean(np.array(cpk)[:, t2_mask]))))
        sim_amp = float(np.mean(np.array(s_pk)[:, t2_mask])) if len(s_pk) > 0 else None
        _reportRange("abs amp (mV, t>2s)", cell_amps, sim_amp)
    print( flush = True )


def statsF( edata, sdata ):
    """Panel F: ISI exp-decay tau per expt cell vs sim."""
    print( "Panel F: ISI decay tau — per-cell comparison" )
    pcm = edata.get('perCellMetrics', {})
    for sq, s_expfit in [("5sq", sdata['expfit5']), ("15sq", sdata['expfit15'])]:
        pk_key = sq[:-2]
        cell_taus = [(c, d[f'tau{pk_key}']) for c, d in pcm.items()]
        sim_tau   = s_expfit['tau'] if s_expfit else None
        print(f"  {sq}:")
        _reportRange("tau (s)", cell_taus, sim_tau)
    print( flush = True )


def statsG( edata, sdata ):
    """Panel G: lognormal amplitude distribution per expt cell vs sim."""
    print( "Panel G: EPSP amplitude distribution — per-cell lognormal comparison" )
    pcm = edata.get('perCellMetrics', {})
    for sq, s_pk in [("5sq", sdata['pk5']), ("15sq", sdata['pk15'])]:
        pk_key = sq[:-2]

        # sim lognormal — absolute peaks
        sv = np.array(s_pk).flatten(); sv = sv[sv > 0]
        if len(sv) >= 5:
            s_s, _, sc_s = stats.lognorm.fit(sv, floc=0)
            sim_mu, sim_sig = float(np.log(sc_s)), float(s_s)
        else:
            sim_mu, sim_sig = None, None

        cell_mus  = [(c, d[f'mu{pk_key}'])  for c, d in pcm.items()]
        cell_sigs = [(c, d[f'sig{pk_key}']) for c, d in pcm.items()]

        print(f"  {sq} — lognormal mu (log-mean):")
        _reportRange("mu", cell_mus, sim_mu)
        print(f"  {sq} — lognormal sigma (log-std):")
        _reportRange("sigma", cell_sigs, sim_sig)
    print( flush = True )


def main():
    parser = argparse.ArgumentParser( description = "Read and plot sim data" )
    parser.add_argument( "-f", "--fname", type = str, help = "Optional: Name of hdf5 pandas file with data.", default = exptData )
    parser.add_argument( "-f2", "--fname2", type = str, help = "Optional: Name of another hdf5 pandas file with data.", default = simData )
    parser.add_argument( "-n", "--noise", type = float, help = "Optional: Noise to add to EPSP trace, in mV. Default = 0.5.", default = 0.5 )
    parser.add_argument( "-freq", "--noiseFreq", type = float, help = "Optional: Bandwidth of noise to add to EPSP trace, in Hz. Default = 200.", default = 200 )
    parser.add_argument( "-t", "--threshold", type = float, help = "Optional: Threshold for classifying an EPSP trace as an event. In mV. Default = 0.5.", default = 0.5 )
    args = parser.parse_args()
    np.random.seed( 12345 )
    df = pandas.read_hdf( args.fname )
    plt.rcParams.update( {"font.size": 20} )
    fig = plt.figure( figsize = (10,24) )
    gs = fig.add_gridspec( 7, 2 ) # 7 rows, 2 cols
    edata = plotFrame( gs, fig, args, df, 0, 0 )
    if args.fname2:
        df2 = pandas.read_hdf( args.fname2 )
        # Here we check if it is real or synth data
        sdata = plotFrame( gs, fig, args, df2, 1, 0)
    statsG( edata, sdata )
    statsD( edata, sdata )
    statsE( edata, sdata )
    statsF( edata, sdata )
    statsC( edata, sdata )
    fig.tight_layout()
    plt.show()


def plotFrame(gs, fig, args, df, column = 0, cell = 0):
    global pulseTrig
    global pulseThresh

    # Set up the stimulus timings
    ax = fig.add_subplot( gs[column,:] )
    panelA_SampleTrace( ax, df, column, args )
    pk5, pk15, fpk5, fpk15, patDict, perCell = scanData( df, args )
    npk5, npk15 = normalizePks( pk5 ), normalizePks( pk15 )
    print( "PULSE TRAIN 32 = ", PulseTrain[NUMPULSE]/ SAMPLE_FREQ )
    dataname = "Expt" if column == 0 else "Sim"
    pt = np.array( PulseTrain ) / SAMPLE_FREQ
    linfit5  = doLinFit( pt[:NUMPULSE], np.array(pk5)[:, :NUMPULSE], dataname +  " 5 sq" )
    linfit15 = doLinFit( pt[:NUMPULSE], np.array(pk15)[:, :NUMPULSE], dataname + " 15 sq" )
    expfit5  = doExpFitISI( np.array(pk5), dataname +  " 5 sq" )
    expfit15 = doExpFitISI( np.array(pk15), dataname + " 15 sq" )
    print( flush = True )

    psd5, psd15, perCellPsd5, perCellPsd15 = panelC_FFT(fig.add_subplot(gs[2,column]), df, column, args)
    prob5, prob15 = panelD_probVsTime( fig.add_subplot(gs[3,column]), column, pk5, pk15 )
    mean5, mean15 = panelE_epspVsTime( fig.add_subplot(gs[4,column]), column, pk5, pk15 )
    fmean5, fmean15, isi = panelF_epspVsISI( fig.add_subplot(gs[5,column]), column, pk5, pk15 )
    panelG_epkHisto( fig.add_subplot(gs[6,column]), column, pk5, pk15 )

    # Compute per-cell scalar metrics (used by stats functions)
    perCellMetrics = _perCellMetrics( perCell ) if column == 0 else {}

    return { 'psd5': psd5, 'psd15': psd15,
             'perCellPsd5': perCellPsd5, 'perCellPsd15': perCellPsd15,
             'pk5': pk5, 'pk15': pk15,
             'npk5': npk5, 'npk15': npk15,
             'prob5': prob5, 'prob15': prob15,
             'mean5': mean5, 'mean15': mean15,
             'linfit5': linfit5, 'linfit15': linfit15,
             'fmean5': fmean5, 'fmean15': fmean15,
             'isi': isi,
             'expfit5': expfit5, 'expfit15': expfit15,
             'perCell': perCell,
             'perCellMetrics': perCellMetrics }


if __name__ == "__main__":
    main()
