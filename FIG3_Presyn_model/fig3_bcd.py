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

def alphaFunc( t, tp ):
    return (t/tp) * np.exp(1-t/tp)

def dualAlphaFunc( t, t1, t2 ):
    if t < 0:
        return 0.0
    if abs( t1 - t2 ) < 1e-6:
        return alphaFunc( t, t1 )
    return (1.0/(t1-t2)) * (np.exp(-t/t1) - np.exp(-t/t2))

class Alpha():
    def __init__( self, kernel ):
        self.kernel = kernel

    def score( self, params ):
        lk = len( self.kernel )
        dt = 1.0/sampleRate
        alpha = np.array([alphaFunc( i*dt, params[0] ) for i in range(lk)] )
        delta = alpha - self.kernel
        return np.dot( delta, delta )

    def fit( self ):
        initGuess = [0.02]
        result = sci.minimize( self.score, initGuess, method = "BFGS" )
        print( "Simple alpha func tp = {:.2f}, score = {:.2f}, initScore = {:.2f}".format( result.x[0], self.score( result.x ), self.score( initGuess ) ) )
        return result.x[0]

    def plotFit( self, freq ):
        plt.figure()
        tp = self.fit()
        t = np.arange( 0.0, 1.0- 1e-6, 1.0/sampleRate )
        plt.plot( t[:len(self.kernel)], self.kernel, label="data" )
        plt.plot( t, [alphaFunc( tt, tp ) for tt in t], label = "fit" )
        plt.xlabel( "Time (s)" )
        plt.ylabel( "frac open" )
        plt.title( "alpha func fit for " + str(freq) + " Hz" )
        plt.legend()

class DualAlpha():
    def __init__( self, kernel, delay, power ):
        self.kernel = kernel
        self.delay = delay
        self.power = power

    def scaledAlphaVec( self, t1, t2, length ):
        dt = 1.0/sampleRate
        alpha = np.array([dualAlphaFunc( i*dt-self.delay, t1, t2 ) for i in range(length)] )
        if abs( max( alpha ) ) > 1e-9:
            alpha = alpha / max( alpha )
        alpha = np.power( alpha, self.power )
        return alpha

    def score( self, params ):
        lk = len( self.kernel )
        dt = 1.0/sampleRate
        alpha = self.scaledAlphaVec( params[0], params[1], lk )
        delta = alpha - self.kernel
        return np.dot( delta, delta )

    def fit( self ):
        initGuess = [0.002, 0.02 ]
        result = sci.minimize( self.score, initGuess, method = "BFGS" )
        print( "delay={}, power={}, t1 = {:.3f}, t2={:.3f}, score = {:.2f}, initScore = {:.2f}".format(self.delay, self.power, result.x[0], result.x[1], self.score( result.x ), self.score( initGuess ) ) )
        return result.x

    def plotFit( self, freq ):
        plt.figure()
        [t1,t2] = self.fit()
        t = np.arange( 0.0, 1.0- 1e-6, 1.0/sampleRate )
        alpha = self.scaledAlphaVec( t1, t2, len(t) )
        plt.plot( t[:len(self.kernel)], self.kernel, label="data" )
        plt.plot( t, alpha, label = "fit" )
        plt.xlabel( "Time (s)" )
        plt.ylabel( "frac open" )
        plt.title( "dual alpha func fit for " + str(freq) + " Hz" )
        plt.title( "dual alpha func freq= {} Hz, delay = {} sec, power = {}".format(freq, self.delay, self.power) )
        plt.legend()


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

    #plt.plot( t, rawKernel )

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
    # This is the kernel change due to the time between valley and peak.
    #print( "riseDelta: ", kpkIdx, stimWidth, riseIdx, kpk )
    riseDelta1 = kernel[kpkIdx + stimWidth - riseIdx ] - kernel[kpkIdx + stimWidth]
    riseDelta = ret[1] - ret[1]*np.exp( -(ret[2] - ret[0]) / tau[1] )
    plt.plot( [ret[2], ret[2]], [-riseDelta + ret[1], ret[3]], "bx-" )
    #print( "origRiseDelta={:.4f}, new={:.4f}".format( riseDelta1, riseDelta) )
    if ret[0] < endT + 0.01:   # First pulse after ref.
        riseTotal = ret[3] - ret[1]
    else:
        riseTotal = riseDelta + ret[3] - ret[1]
        #riseTotal = ret[3] - ret[1]
    #print( "si, kpkidx, pkIdx, riseDelta, risetotal = ", int( sampleRate/stimWidth ), si, kpkIdx, pkIdx, riseDelta, riseTotal/kpk )
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

def deconv( dat, freq ):
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
        #ret = findPkVal( dat, freq, si, (kpk < 0) )
        pv.append( ret )
        #print( "FindSTPSCALE: ", kpk, ret, si, stimWidth )
        scale = findStpScale( kernel, kpk, ret, si, stimWidth, tau )
        #print( si, scale )
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
    ret[int(round(sampleRate* endT )):] += npv[1][1]
    #print( "LK = ", len( kernel ), len( ret ), stimIdx, scaleList )
    for ss, idx in zip( scaleList, stimIdx ):
        #idx -= int( round( endT * sampleRate) )
        #print( "SS, IDX = ", ss, idx )
        if idx > 0 :
            ret[idx:len(kernel)+idx] += kernel * ss

    t = np.arange( 0.0, 1.0 - 1e-6, 1.0/sampleRate )
    #t2 = np.insert( t, 0, startT )
    #t2 = np.array( [0.2]+[endT + ii/freq for ii in range(8)] )
    #t2 = np.array( [0.2]+[ endT + ii/sampleRate for ii in stimIdx[1:]] )
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
    pattern = int( patternList[0] )
    finh, fexc = innerAnalysis( sqDat, freq, pattern, cellNum )
    if doFigs:
        plt.show()

if __name__ == "__main__":
    main()
