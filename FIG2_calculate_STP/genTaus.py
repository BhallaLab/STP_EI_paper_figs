import pandas
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math
import scipy.optimize as sci

datadir = "../../../2022/VC_DATA"
#datadir = "/home1/bhalla/adityaa/Lab/Projects/EI_Dynamics/Analysis/parsed_data"
sampleRate = 20000.0
OnsetDelay = 0.025  # Delay from the light stimulus to peak of syn current
stimT = 0.2
startT = stimT + OnsetDelay
endT = 0.5
runtime = 1.0

def tauFit( kernel, baseline):
    y = np.array(kernel) - baseline
    pk = y[0]
    x = np.linspace( 0, len(y)/sampleRate, len(y), endpoint = False )
    #print( "LEN KERNEL = ", len( kernel ), len( y ), len( x ) )
    ret, cov = sci.curve_fit(lambda t,a,tau: a*np.exp(-t/tau), 
            x, y, p0=(pk,0.02) )
    #print( "len = {}, a = {:.4f}, tau = {:4f}, range = {:4f}, bl = {:.4f}".format( len( y ), ret[0], ret[1], pk - baseline, baseline ) )
    #plt.plot( x + startT + 0.05, ret[0] * np.exp( -x/ret[1] ), ":" )
    return ret

def calcKernel( cellID, dat ):
    dataStartColumn = len( dat.columns ) - 80000 # Was 29, new is 39.
    startIdx = int( round(startT*sampleRate ) ) + dataStartColumn
    stimIdx = int( round(stimT*sampleRate ) ) + dataStartColumn
    endIdx = int( round(endT*sampleRate ) ) + dataStartColumn
    if len( dat ) == 0:
        return None, None, None, [None,None]
    baseline  = np.mean( np.array(dat.iloc[:,stimIdx - int( 0.005*sampleRate ):stimIdx] ) )
    rawKernel  = np.mean( dat.iloc[:,startIdx:endIdx] )
    kmax = max( rawKernel )
    kmin = min( rawKernel )
    if math.isnan( kmax ) or math.isnan( kmin ):
        raise FloatingPointError( "calcKernel: kmax or kmin is a nan" )


    if abs( kmax ) > abs( kmin ): # Inhib trace has a positive peak.
        return kmax, rawKernel, baseline, tauFit(rawKernel, baseline)
    else:
        return kmin, rawKernel, baseline, tauFit(rawKernel, baseline)



def main():
    parser = argparse.ArgumentParser( description = "This program analyzes Aditya's voltage clamp data from a pandas file" )
    parser.add_argument( "-o", "--output", type = str, help = "Optional: output pandas hdf file for model params, default = STP_pks_and_refs.h5", default = "STP_pks_and_refs.h5")
    args = parser.parse_args()

    frame = []
    alldat = pandas.read_hdf( "{}/all_cells_FreqSweep_VC_long.h5".format( datadir ) )
    alldat['patternList'] = alldat['patternList'].astype( int )
    #print( alldat.columns[:100] )
    #cellNum = 7492
    cellList = alldat['cellID'].unique()

    ########## Set up graphics #########
    fig = plt.figure( figsize = (16,20) )
    gs = fig.add_gridspec( 10, 4 ) # 10 rows, 4 cols
    plt.rcParams.update({'font.size': 14})

    ######### Calculate taus ###########
    idx = 0
    for cellID in cellList:
        ax = fig.add_subplot( gs[ idx // 4, idx % 4 ] )
        dat = alldat.loc[ alldat["cellID"] == cellID ]
        excDat = dat.loc[ dat['clampPotential'] < -0.05 ]
        kmax, rawKernel, baseline, [ampl,tau] = calcKernel( cellID, excDat )
        if kmax != None:
            print( "{}  EPSC  a={:.3f}  t={:.1f}".format( cellID, ampl, tau*1000 ) )
            #excTau = analyzeTau( excDat, cellID, "exc" )
            x = np.arange( len(rawKernel) )/sampleRate
            ax.plot( x + startT, rawKernel, "-b" )
            if tau > 0.001:
                taustr = "  tau={:.1f}ms".format( tau * 1000 ) 
                ax.plot( x + startT, baseline + ampl * np.exp(-x/tau), ":b")
            ax.text( -0.12, 1.05, str( cellID) + taustr, fontsize = 16, weight = "bold", transform=ax.transAxes )
        idx += 1
        ax = fig.add_subplot( gs[ idx // 4, idx % 4 ] )
        inhDat = dat.loc[ dat['clampPotential'] > -0.05 ]
        kmax, rawKernel, baseline, [ampl,tau] = calcKernel( cellID, inhDat )
        if kmax != None:
            print( "{}  IPSC  a={:.3f}  t={:.1f}".format( cellID, ampl, tau*1000 ) )
            x = np.arange( len(rawKernel) )/sampleRate
            ax.plot( x + startT, rawKernel, "-r" )
            if tau > 0.001:
                taustr = " tau={:.1f}ms".format( tau * 1000 ) 
                ax.plot( x + startT, baseline + ampl * np.exp(-x/tau), ":r")
            ax.text( -0.12, 1.05, str( cellID) + taustr, fontsize = 16, weight = "bold", transform=ax.transAxes )
        idx += 1
        #inhTau = analyzeTau( inhDat, cellID, "inh" )

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
