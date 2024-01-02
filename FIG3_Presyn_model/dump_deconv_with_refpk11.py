import pandas
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math
import scipy.optimize as sci

datadir = "../../../2022/VC_DATA"
#datadir = "/home1/bhalla/adityaa/Lab/Projects/EI_Dynamics/Analysis/parsed_data"
sampleRate = 20000.0
runtime = 1.0
doFigs = False
OnsetDelay = 0.009  # Delay from the light stimulus to onset of syn current
startT = 0.2 + OnsetDelay
endT = 0.5 + OnsetDelay

def tauFit( kernel, baseline):
    y = kernel[ int( round(0.05*sampleRate )): ]
    #print( "KERN = ", y.shape, len( y ), min(y), max(y) )
    pk = y[0]
    x = np.linspace( 0, len(y)/sampleRate, len(y), endpoint = False )
    '''
    plt.figure()
    plt.plot( x, y )
    plt.show()
    '''
    ret, cov = sci.curve_fit(lambda t,a,tau: a*np.exp(-t/tau), x, y, p0=(pk-baseline,0.02) )
    return ret

def calcKernel( dat ):

    startIdx = int( round(startT*sampleRate ) )
    endIdx = int( round(endT*sampleRate ) )
    baseline  = np.mean( dat.iloc[startIdx - int( 0.005*sampleRate ):startIdx] )
    
    rawKernel = np.array( dat.iloc[startIdx:endIdx] )
    kmax = max( rawKernel )
    kmin = min( rawKernel )
    if math.isnan( kmax ) or math.isnan( kmin ):
        raise FloatingPointError( "calcKernel: kmax or kmin is a nan" )

    if abs( kmax ) > abs( kmin ): # Inhib trace has a positive peak.
        return kmax, rawKernel, baseline, tauFit( rawKernel, baseline )
    else:
        return kmin, rawKernel, baseline, tauFit( rawKernel, baseline )

def findStpScale( kernel, kpk, ret, si, stimWidth, tau ):
    if ret[0] < endT and si < (endT * sampleRate):
        return 1.0
    if kpk < 0 : # Exc
        kpkIdx = np.argmin( kernel[:-stimWidth] )
    else:
        kpkIdx = np.argmax( kernel[:-stimWidth] )
    pkIdx = int( round ( (ret[2]) * sampleRate )) - si
    riseIdx = int( round ( (ret[2] - ret[0]) * sampleRate ))
    riseDelta1 = kernel[kpkIdx + stimWidth - riseIdx ] - kernel[kpkIdx + stimWidth]
    riseDelta = ret[1] - ret[1]*np.exp( -(ret[2] - ret[0]) / tau[1] )
    if ret[0] < endT + 0.01:   # First pulse after ref.
        riseTotal = ret[3] - ret[1]
    else:
        riseTotal = riseDelta + ret[3] - ret[1]
    return riseTotal / kpk

def findPkVal( dat, freq, startIdx, isExc ):
    # Returns times and absolute peak and val preceding peak of specified stim.
    stimWidth = int( round( 0.7 * sampleRate / freq ) )
    d2 = np.array(dat.iloc[startIdx:startIdx + stimWidth])
    if isExc: # Sign flipped in current. imin is peak.
        imin = np.argmin(d2)
        # Look for valley preceding this.
        d3 = np.array(dat[startIdx + imin-stimWidth:imin + startIdx])
        imax = np.argmax(d3)
        if imax + imin - stimWidth > imin:
            print( "WARNING: reversal" )
        return [(imax + startIdx + imin - stimWidth)/sampleRate, d3[imax], (startIdx+imin)/sampleRate, d2[imin]]
    else:   # This is the inhibitory input, which has positive currents.
        imax = np.argmax(d2)
        d3 = np.array(dat.iloc[startIdx + imax-stimWidth:imax + startIdx])
        imin = np.argmin(d3)
        if imax < imin + imax - stimWidth:
            print( "WARNING: reversal" )
        return [(imin + startIdx+imax - stimWidth )/sampleRate, d3[imin], (startIdx + imax)/sampleRate, d2[imax]]

def deconv( dat, freq ):
    startIdx = int( round( 0.5 * sampleRate ) )
    stimWidth = int( round( sampleRate/freq ) )
    # Don't reuse stimWidth for stimIdx because of cumulative error.
    stimIdx = [int(0.2 * sampleRate)] + [ int( round( sampleRate* (0.5 + i/freq ) ) ) for i in range( 8 ) ]

    kpk, kernel, baseline, tau = calcKernel( dat )
    kpkidx = np.argmax( kernel ) if kpk > 0 else np.argmin( kernel )

    scaleList = []
    absPk = [ kpk ]
    absVal = [ baseline ]
    pv = []
    correctedStimIdx = []

    for si in stimIdx:
        ret = findPkVal( dat, freq, si + kpkidx//2, (kpk < 0) )
        pv.append( ret )
        scale = findStpScale( kernel, kpk, ret, si, stimWidth, tau )
        scaleList.append( scale )
        if kpk > 0:
            label = "Inh"
        else:
            label = "Exc"

    npv = np.array( pv ).transpose()
    synthPlot = plotFromKernel(scaleList, stimIdx, kernel, freq, npv, label)
    return np.array(scaleList), synthPlot, npv

def plotFromKernel( scaleList, stimIdx, kernel, freq, npv, label ):
    ret = np.zeros( int( round( sampleRate * 1.5 ) ) ) 
    ret[int(round(sampleRate*endT)):] += npv[1][1]
    #print( "LK = ", len( kernel ), len( ret ), stimIdx, scaleList )
    for ss, idx in zip( scaleList, stimIdx ):
        #idx -= int( round( 0.5 * sampleRate) )
        #print( "SS, IDX = ", ss, idx )
        if idx > 0 :
            ret[idx:len(kernel)+idx] += kernel * ss

    t = np.arange( 0.0, 1.0 - 1e-6, 1.0/sampleRate )
    #t2 = np.insert( t, 0, 0.2 )
    #t2 = np.array( [0.2]+[0.5 + ii/freq for ii in range(8)] )
    #t2 = np.array( [0.2]+[0.5 + ii/sampleRate for ii in stimIdx[1:]] )
    plt.plot( t, ret[0:len(t)], ":", label = label + "_est" )
    #print( "npv shape = ", npv.shape, ", len scaleList = ", len( scaleList ) )
    plt.plot( npv[0], npv[1], "c*-" )
    plt.plot( npv[2], npv[3], "y*-" )
    return ret[:len(t)]

def innerAnalysis( sqDat, freq, pattern, cellNum ):
    numSamples = int( round( sampleRate * runtime ) )
    #dataStartColumn = len( sqDat.columns ) - 80000 # Was 29, new is 39.
    dataStartColumn = len( sqDat.columns ) - 80000 # Was 29, new is 39.
    #print( " dataStartColumn = ", dataStartColumn, len( sqDat.columns ), numSamples )
    inh = sqDat.loc[ (sqDat['stimFreq'] == freq) & (sqDat['clampPotential'] > -0.05 ) & (sqDat['patternList'] == pattern ) ]
    exc = sqDat.loc[ (sqDat['stimFreq'] == freq) & (sqDat['clampPotential'] < -0.05 ) & (sqDat['patternList'] == pattern ) ]
    #print( mypat )
    #yexc = exc.iloc[:,numSamples + dataStartColumn: 2*numSamples+dataStartColumn]
    #yinh = inh.iloc[:,numSamples + dataStartColumn: 2*numSamples+dataStartColumn]
    yexc = exc.iloc[:,dataStartColumn: numSamples+dataStartColumn]
    yinh = inh.iloc[:,dataStartColumn: numSamples+dataStartColumn]
    emean = yexc.mean()
    imean = yinh.mean()
    #print( " LENS = ", len( emean ), len( imean ) )

    if doFigs:
        plt.figure()
    try:
        ideconv, iSynthPlot, ipv = deconv( imean, freq )
        edeconv, eSynthPlot, epv = deconv( emean, freq )
    except FloatingPointError as error:
        print( "innerAnalysis: freq = {}, pattern = {}, cellNum = {}".format( freq, pattern, cellNum ) )
        raise
    assert( len( eSynthPlot ) == len( emean ) )
    assert( len( iSynthPlot ) == len( imean ) )

    if doFigs:
        t = np.arange( 0.0, runtime - 1e-6, 1.0/sampleRate )
        plt.plot( t, imean, label="Inh" )
        plt.plot( t, emean, label="Exc" )

        plt.xlabel( "Time (s)" )
        plt.ylabel( "Synaptic current (pA)" )
        plt.legend()
        plt.title( "cell {}, Freq = {}".format( cellNum, freq ) )

        plt.figure()
        tr = np.arange( 0.0, len( imean ) - 1e-6, 1.0 ) / sampleRate
        plt.plot( tr, iSynthPlot - imean, label = "Inh" )
        tr = np.arange( 0.0, len( emean ) - 1e-6, 1.0 ) / sampleRate
        plt.plot( tr, eSynthPlot - emean, label = "Exc" )

        plt.xlabel( "Time (s)" )
        plt.ylabel( "Residual (pA)" )
        plt.legend()
        plt.title( "cell {}, Freq = {}: Residual ".format( cellNum, freq ) )

        plt.figure()
        plt.plot( range( len( ideconv ) ), ideconv/ideconv[0], label="Inh" )
        plt.plot( range( len( edeconv ) ), edeconv/edeconv[0], label="Exc" )
        plt.xlabel( "Pulse # in burst" )
        plt.ylabel( "Synaptic response as ratio to first pulse." )
        plt.legend()
        plt.title("cell {}, Freq = {}: Deconv ".format( cellNum, freq) )
        plt.show()

    finh = np.append( ideconv, ipv ) 
    fexc = np.append( edeconv, epv )
    return finh, fexc

def badData( cell, freq, pattern ):
    '''
    bad = [(6301, 40, 48), (6301, 40, 49), (6301, 40, 47), (6301, 40, 46), 
            (6301, 40, 50), (6301, 40, 55), (6301, 40, 52), (6301, 40, 53),
            (6201, 40, 48), (6201, 40, 49), (6201, 40, 47), (6201, 40, 46),
            (6201, 40, 50), (6201, 40, 55), (6201, 40, 52), (6201,40,53), 
            (1523, 20, 48), (1523, 20, 49), (1523, 20, 47), (1523, 20, 46),
            (1523, 20, 50), (1523, 20, 55), (1523, 20, 52), (1523, 20, 53), 
            (1523, 50, 48), (1523, 50, 49), (1523, 50, 47), (1523, 50, 46),
            (1523, 50, 50), (1523, 50, 55), (1523, 50, 52), (1523, 50, 53)
            ]
    return (cell, freq, pattern) in bad
    '''
    #badFreq = {7491: 30, 7491:40, 6301:40, 6201:40, 1523:20, 1523:50, 1931:20, 1931:50, 1621:30}
    badFreq = { 7491: {20:[50], 30:[47,48,49], 50:[48,49,47]}, 
            6301: {40:-1, 50:-1, 30:[55,53,48,49,50] },
            6201: {40:-1, 30:[48,49,50]},
            1931: {20:-1, 50:-1, 30:[49,50,55]},
            1621: {30:[55,52,53], 40:[52,53,55]},
            1541: {50:[49,46]},
            1531: {50:[46]},
            1523: {20:-1, 50:-1},
            1522: {20:[48,47,52,53,55], 50:-1},
            111: {20:[49], 50:[46,47,48,49,55]},
            1491: {20:[48,47,55,52,53], 50:[48,49,47,50,55,52,53]}
            } 
    bb = badFreq.get(cell)
    if bb:
        ff = bb.get(freq)
        if ff:
            if ff == -1: # Indicates all patterns
                return True
            elif pattern in ff:
                return True
    return False


def main():
    global doFigs
    parser = argparse.ArgumentParser( description = "This program analyzes Aditya's voltage clamp data from a pandas file" )
    parser.add_argument( "-o", "--output", type = str, help = "Optional: output pandas hdf file for model params, default = STP_pks_and_refs.h5", default = "STP_pks_and_refs.h5")

    parser.add_argument( "-f", "--doFigs", action = "store_true", help = "Optional Flag: Decides if figures are to be displayed. Default False" )
    args = parser.parse_args()
    doFigs = args.doFigs

    frame = []
    alldat = pandas.read_hdf( "{}/all_cells_FreqSweep_VC_long.h5".format( datadir ) )
    alldat['patternList'] = alldat['patternList'].astype( int )
    #print( alldat.columns[:100] )
    #for cellNum in [7491, 6301, 6201, 1931, 1621, 1541, 1531, 1524, 1523, 1522, 1491, 111]:
    for cellNum in [ 7492, 7491, 6301, 6201, 1931, 1621, 1541, 1531, 1524, 1523, 1522, 1491, 111]:
        #dat = pandas.read_hdf( "{}/cell{}_trainingSet_longest.h5".format( datadir, cellNum ) )
        dat = alldat.loc[ alldat["cellID"] == cellNum ]
        pattern = 0
        for numSq in [5, 15]:
            sqDat = dat.loc[ dat['numSq'] == numSq ]
            freqList = sqDat['stimFreq'].unique()
            #sqDat['patternList'] = sqDat['patternList'].astype( int )
            #print (sqDat.columns[:51])
            #print (sqDat.columns)
            #print ("SHAPE = ", sqDat['patternList'].shape )
            '''
            print ( "PL", sqDat['patternList'][:4] )
            print ( "eexptid", sqDat['exptID'][:4] )
            print ( "numpatterns", sqDat['numPatterns'][:4] )
            print ( "exptSeq", sqDat['exptSeq'][:4] )
            print ( "PL0", sqDat['patternList'].iloc[0] )
            print ( "PL0.shape", sqDat['patternList'].iloc[0].shape )
            '''

            patternList = sqDat['patternList'].unique()
            print( 'patternList = ', patternList )
            for freq in freqList:
                for pattern in patternList:
                    pattern = int( pattern )
                    if badData( cellNum, freq, pattern ):
                        continue
                    print( cellNum, numSq, freq, pattern )
                    finh, fexc = innerAnalysis( sqDat, freq, pattern, cellNum )
                    line = [ cellNum, 0, numSq, freq, pattern ]
                    frame.append( np.append( line, finh ) )
                    line = [ cellNum, 1, numSq, freq, pattern ]
                    frame.append( np.append( line, fexc ) )
                    if doFigs:
                        plt.show()

    print( "Finished pre-analysis" )
    colNames =  ["cell", "exc", "numSq", "freq", "pattern"] + \
        [ "p"+str(ii) for ii in range(9) ] + \
        [ "tval"+str(ii) for ii in range(9) ] +  \
        [ "val"+str(ii) for ii in range(9) ] + \
        [ "tpk"+str(ii) for ii in range(9) ] + \
        [ "pk"+str(ii) for ii in range(9) ]
    assert( len( colNames ) == len( frame[0] ) )
    outputDF = pandas.DataFrame( frame, columns = colNames )
    outputDF.to_hdf( args.output, "w" )

    #Id = dat['exptID'].loc[ dat['numSq'] == numSq ]
    #sweep = dat['sweep'].loc[ dat['numSq'] == numSq ]
    #Rm = dat['InputRes'].loc[ dat['numSq'] == numSq ]
    #tau = dat['Tau'].loc[ dat['numSq'] == numSq ]


if __name__ == "__main__":
    main()
