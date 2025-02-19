import numpy as np
import matplotlib.pyplot as plt
import rdesigneur as rd
import pandas as pd
chemDt = 1e-4
chemPlotDt = 2e-4
elecDt = 50e-6
diffusionLength = 1e-3 # ensure single chem compt
stimTime = 0.002
eps = 1e-6
presettleTime = 0.35
settleTime = 0.65
stimCa = 50e-3  # 50 uM
baseCa = 80e-6  # 80 nM
doInnerPlot = False

runtime = 0.8

def calcSTP( cell, isExc, freqdf ):
    pkNames = ['pk{}'.format(ii) for ii in range(9)]
    tpkNames = ['tpk{}'.format(ii) for ii in range(9)]
    valNames = ['val{}'.format(ii) for ii in range(9)]
    tvalNames = ['tval{}'.format(ii) for ii in range(9)]
    edf = freqdf.loc[freqdf["exc"] == int( isExc )]
    meanVal = edf[valNames].mean(axis=0)
    meanTval = edf[tvalNames].mean(axis=0)
    meanPk = edf[pkNames].mean(axis=0)
    meanTpk = edf[tpkNames].mean(axis=0)

    plotvec = np.dstack((meanVal, meanPk)).flatten()
    tvec = np.dstack((meanTval, meanTpk)).flatten()
    #print( "SZ = ", plotvec.shape, tvec.shape )
    val1 = meanVal[1]
    pk1 = meanPk[1] - val1
    pk2 = meanPk[2] - val1
    val2 = meanVal[2] - val1
    if doInnerPlot:
        xp0 = meanTpk[0]
        xp1 = meanTpk[1]
        xv2 = meanTval[2]
        xp2 = meanTpk[2]
        colors = ["red", "blue", "green"]
        plt.figure()
    if isExc:
        #val1 = meanVal[1]
        #pk1 = -meanPk[1]
        #pk2 = -meanPk[2]
        #val2 = -meanVal[2]
        if cell == 7492:
            freq = freqdf['freq'].unique()
            print( "{}  p1={:.3f}  v2={:.3f}  p2={:.3f}  stp={:.3f}".format( freq, pk1/pk2, val2/pk2, 1, (pk2-val2)/pk1 ) )
        if doInnerPlot:
            plt.scatter( [xp1, xv2, xp2], [pk1, val2, pk2], color=colors, marker = "o")
            #plt.xlim( settleTime, min( runtime, settleTime + isi*2+0.005 ) )
            #plt.ylim( 0, (min(-pk2, -pk1) * 1.1 ) )
        #pk1 -= meanVal[1]
        #pk2 -= meanVal[2]
    else:
        #pk1 = meanPk[1]
        #pk2 = meanPk[2]
        #val2 = meanVal[2]
        if doInnerPlot:
            plt.scatter( [xp1, xv2, xp2], [pk1, val2, pk2], color=colors, marker = "o")
            #plt.xlim( settleTime, runtime )
            #plt.ylim( 0, (max(pk2, pk1) * 1.1 ) )
        #pk1 -= meanVal[1]
        #pk2 -= meanVal[2]

    #print( "{}  {}: {:.3g}  {:.3g}  {:.3g}  {:.4g}".format( cell, isExc, pk1, pk2, val1, baseline ) )
    if doInnerPlot:
        plt.title( "{} {}".format( cell, "Exc" if isExc else "Inh" ) )
        plt.plot( tvec, plotvec )
        plt.show()
    return pk1, (pk2-val2)

def main():
    dat = pd.read_hdf( "STP_pks_and_refs.h5")
    cellDict = {
        # key: (gluTau, gabaTau),
        1931: (8.9, 15.0),
        1491: (10.0, 10.0),
        1541: (10.0, 10.0),
        1524: (10.0, 15.0),
        #1523: (14.5, None),
        1522: (10.0, 16.0),
        1531: (9.6, 10.0),
        1621: (10.0, 10.0),
        111: (12.0, 12.0),
        7492: (15, 25.0),
        7491: (15.0, 16.0),
        6201: (9.0, 16.0),
        6301: (15.0, 20.0),
        #5501: (1701.9, None)

    }


    fig, axes = plt.subplots(4, 3, figsize=(10, 15.0))
    axes = axes.flatten()
    isiList = [ 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.08 ]
    idx = 0
    numSq = 5
    for cell in cellDict:
        celldf = dat.loc[dat['cell']==cell]
        freqList = sorted( celldf['freq'].unique() )
        isiList = [1.0/ff for ff in freqList]
        ax = axes[idx]
        eplot = []
        iplot = []
        for freq, isi in zip( freqList, isiList):
            freqdf = celldf.loc[(celldf['freq']==freq) & (celldf["numSq"]==numSq)]
            pk1, pk2 = calcSTP( cell, True, freqdf )
            eplot.append( pk2/pk1 )
            pk1, pk2 = calcSTP( cell, False, freqdf )
            iplot.append( pk2/pk1 )
            print( ".", flush=True, end ="" )
        ax.plot( isiList, eplot, "-or", label = "gluR" )
        ax.plot( isiList, iplot, "-ob", label = "GABAR" )
        ax.set_ylim( 0, 2.5 )
        ax.set_title( cell )
        ax.legend()
        print( cell )
        idx += 1
    plt.show()


if __name__ == "__main__":
    main()

