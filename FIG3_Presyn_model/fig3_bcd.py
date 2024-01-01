import pandas
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math
import scipy.optimize as sci

datadir = "../../../2022/VC_DATA"
#datadir = "/home1/bhalla/adityaa/Lab/Projects/EI_Dynamics/Analysis/parsed_data"
sampleRate = 20000.0
OnsetDelay = 0.009  # Delay from the light stimulus to onset of syn current
startT = 0.2 + OnsetDelay
endT = 0.5 + OnsetDelay
runtime = 1.0
doFigs = True

def tauFit( kernel, baseline):
    y = kernel[ int( round(0.05*sampleRate )): ]
    pk = y[0]
    x = np.linspace( 0, len(y)/sampleRate, len(y), endpoint = False )
    ret, cov = sci.curve_fit(lambda t,a,tau: a*np.exp(-t/tau), x, y, p0=(pk-baseline,0.02) )
    #print( "len = {}, a = {:.4f}, tau = {:4f}, range = {:4f}, bl = {:.4f}".format( len( y ), ret[0], ret[1], pk - baseline, baseline ) )
    plt.plot( x + startT + 0.05, ret[0] * np.exp( -x/ret[1] ), ":" )
    return ret

def calcKernel( dat ):

    startIdx = int( round(startT*sampleRate ) )
    endIdx = int( round(endT*sampleRate ) )
    #baseline  = np.mean( dat[startIdx - int( 0.005*sampleRate ):startIdx] )
    baseline  = np.mean( dat.iloc[startIdx - int( 0.005*sampleRate ):startIdx] )
    
    rawKernel = np.array( dat.iloc[startIdx:endIdx] )
    #t = np.arange( 0.0, endT - startT - 1e-6, 1.0/sampleRate )
    #print( "SHAPE = ", dat.shape, rawKernel.shape, startIdx, endIdx )
    kmax = max( rawKernel )
    kmin = min( rawKernel )
    if math.isnan( kmax ) or math.isnan( kmin ):
        #print( "idx = ", startIdx, endIdx )
        #print( dat.iloc[startIdx:endIdx] )
        raise FloatingPointError( "calcKernel: kmax or kmin is a nan" )


    if abs( kmax ) > abs( kmin ): # Inhib trace has a positive peak.
        return kmax, rawKernel, baseline, tauFit(rawKernel, baseline)
    else:
        return kmin, rawKernel, baseline, tauFit(rawKernel, baseline)

def findStpScale( kernel, kpk, ret, si, stimWidth, tau ):
    # ret[0,1] = valley_t, y; ret[1,2] = pk_t, y
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
    label = "Min to Max" if (si < 11000 and kpk > 0) else None
    plt.plot( [ret[2], ret[2]], [-riseDelta + ret[1], ret[3]], "ro-", label = label )
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
        d3 = np.array(dat.iloc[startIdx + imin-stimWidth:imin + startIdx])
        #d3 = np.array(dat.iloc[startIdx:startIdx+stimWidth])
        #print( "Exc startIdx={}, imin = {}, stimWidth = {}, tot={}".format( startIdx, imin, stimWidth, startIdx + imin-stimWidth))
        imax = np.argmax(d3)
        if imax + imin - stimWidth > imin:
            print( "WARNING: reversal" )
        return [(imax + startIdx + imin - stimWidth)/sampleRate, d3[imax], (startIdx+imin)/sampleRate, d2[imin]]
    else:   # This is the inhibitory input, which has positive currents.
        imax = np.argmax(d2)
        d3 = np.array(dat.iloc[startIdx + imax-stimWidth:imax + startIdx])
        #print( "Inh startIdx={}, imax = {}, stimWidth = {}, tot={}".format( startIdx, imax, stimWidth, startIdx + imax-stimWidth))
        imin = np.argmin(d3)
        if imax < imin + imax - stimWidth:
            print( "WARNING: reversal" )
        return [(imin + startIdx+imax - stimWidth )/sampleRate, d3[imin], (startIdx + imax)/sampleRate, d2[imax]]

def deconv( dat, freq, ax ):
    startIdx = int( round( endT * sampleRate ) )
    stimWidth = int( round( sampleRate/freq ) )
    # Don't reuse stimWidth for stimIdx because of cumulative error.
    stimIdx = [int(startT * sampleRate)] + [ int( round( sampleRate* (endT + i/freq ) ) ) for i in range( 8 ) ]

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
    synthPlot = plotFromKernel(scaleList, stimIdx, kernel, freq, npv, label, ax)
    return np.array(scaleList), synthPlot, npv

def plotFromKernel( scaleList, stimIdx, kernel, freq, npv, label, ax ):
    ret = np.zeros( int( round( sampleRate * 1.5 ) ) ) 
    ret[int(round(sampleRate* endT )):] += npv[1][1]
    for ii in range( len( scaleList ) ):
        ss = scaleList[ii]
        idx = stimIdx[ii]
        if idx > 0 :
            ks = kernel * ss
            if label == "Inh":
                offset = npv[3,ii] - max( ks + ret[idx:len(kernel)+idx] )
            else: 
                offset = npv[3,ii] - min( ks + ret[idx:len(kernel)+idx] )
            ret[idx:len(kernel)+idx] += ks + offset

    t = np.arange( 0.0, 1.0 - 1e-6, 1.0/sampleRate )
    #plt.plot( t, ret[0:len(t)], ":", label = label + "_est" )
    el1 = None if label == "Inh" else "Troughs"
    el2 = None if label == "Inh" else "Peaks"
    ax.plot( npv[0], npv[1], "c*-", label = el1 )
    ax.plot( npv[2], npv[3], "y.-", label = el2 )
    return ret[:len(t)]

def innerAnalysis( sqDat, freq, pattern, cellNum ):
    numSamples = int( round( sampleRate * runtime ) )
    dataStartColumn = len( sqDat.columns ) - 80000 # Was 29, new is 39.
    inh = sqDat.loc[ (sqDat['stimFreq'] == freq) & (sqDat['clampPotential'] > -0.05 ) & (sqDat['patternList'] == pattern ) ]
    exc = sqDat.loc[ (sqDat['stimFreq'] == freq) & (sqDat['clampPotential'] < -0.05 ) & (sqDat['patternList'] == pattern ) ]
    yexc = exc.iloc[:,dataStartColumn: numSamples+dataStartColumn]
    yinh = inh.iloc[:,dataStartColumn: numSamples+dataStartColumn]
    emean = yexc.mean()
    imean = yinh.mean()

    if doFigs:
        fig = plt.figure( figsize = (10,12) )
        gs = fig.add_gridspec( 3, 2 ) # 3 rows, 2 cols
        plt.rcParams.update({'font.size': 14})
        ax = fig.add_subplot( gs[0,:] )
    try:
        ideconv, iSynthPlot, ipv = deconv( imean, freq, ax )
        edeconv, eSynthPlot, epv = deconv( emean, freq, ax )
    except FloatingPointError as error:
        print( "innerAnalysis: freq = {}, pattern = {}, cellNum = {}".format( freq, pattern, cellNum ) )
        raise
    assert( len( eSynthPlot ) == len( emean ) )
    assert( len( iSynthPlot ) == len( imean ) )

    if doFigs:
        t = np.arange( 0.0, runtime - 1e-6, 1.0/sampleRate )
        ax.plot( t, imean, "g-", label="IPSC" )
        ax.plot( t, emean, "b-", label="EPSC" )

        ax.set_xlabel( "Time (s)" )
        ax.set_ylabel( "Synaptic current (pA)" )
        ax.set_ylim( -300, 800 )
        ax.legend( loc="upper left", frameon = False )
        ax.text( -0.12, 1.05, "A", fontsize = 22, weight = "bold", transform=ax.    transAxes )
        #ax.title( "cell {}, Freq = {}".format( cellNum, freq ) )

        '''
        ax = fig.add_subplot( gs[1,0] )
        tr = np.arange( 0.0, len( imean ) - 1e-6, 1.0 ) / sampleRate
        ax.plot( tr, iSynthPlot - imean, label = "Inh" )
        tr = np.arange( 0.0, len( emean ) - 1e-6, 1.0 ) / sampleRate
        ax.plot( tr, eSynthPlot - emean, label = "Exc" )

        plt.xlabel( "Time (s)" )
        plt.ylabel( "Residual (pA)" )
        plt.legend( loc="upper left", frameon = False )
        plt.title( "cell {}, Freq = {}: Residual ".format( cellNum, freq ) )
        '''
        ax = fig.add_subplot( gs[1,0] )
        ax.plot( range( len( ideconv ) ), ideconv/ideconv[0], label="Inh" )
        ax.plot( range( len( edeconv ) ), edeconv/edeconv[0], label="Exc" )
        ax.set_xlabel( "Pulse # in burst" )
        ax.set_ylabel( "Min-to-Max ratio" )
        ax.legend( loc="upper right", frameon = False )
        ax.set_ylim( 0, 1.4)
        ax.text( -0.24, 1.05, "B", fontsize = 22, weight = "bold", transform=ax.    transAxes )

        ax = fig.add_subplot( gs[1,1] )
        y = ipv[0] - ipv[2]
        ax.plot( range( len( y ) ), y/y[0], label="Inh" )
        y = epv[0] - epv[2]
        ax.plot( range( len( y ) ), y/y[0], label="Exc" )
        ax.set_xlabel( "Pulse # in burst" )
        ax.set_ylabel( "Peak-to-Trough ratio" )
        ax.set_ylim( 0, 1.4)
        ax.legend( loc="upper right", frameon = False )
        ax.text( -0.24, 1.05, "C", fontsize = 22, weight = "bold", transform=ax.    transAxes )
        #ax.title("cell {}, Freq = {}: Deconv ".format( cellNum, freq) )
        plt.tight_layout()
        plt.show()

    finh = np.append( ideconv, ipv ) 
    fexc = np.append( edeconv, epv )
    return finh, fexc

def main():
    global doFigs
    parser = argparse.ArgumentParser( description = "This program analyzes Aditya's voltage clamp data from a pandas file" )
    parser.add_argument( "-o", "--output", type = str, help = "Optional: output pandas hdf file for model params, default = STP_pks_and_refs.h5", default = "STP_pks_and_refs.h5")
    args = parser.parse_args()

    frame = []
    alldat = pandas.read_hdf( "{}/all_cells_FreqSweep_VC_long.h5".format( datadir ) )
    alldat['patternList'] = alldat['patternList'].astype( int )
    #print( alldat.columns[:100] )
    cellNum = 7492
    numSq = 5
    dat = alldat.loc[ alldat["cellID"] == cellNum ]
    sqDat = dat.loc[ dat['numSq'] == numSq ]
    freqList = sqDat['stimFreq'].unique()
    patternList = sqDat['patternList'].unique()
    freq = freqList[0]
    pattern = int( patternList[1] )
    finh, fexc = innerAnalysis( sqDat, freq, pattern, cellNum )
    if doFigs:
        plt.show()

if __name__ == "__main__":
    main()
