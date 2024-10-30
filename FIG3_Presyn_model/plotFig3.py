##
## plotFig3.py
## This file generates Figure 3 panels B to G.
##
## It is licenced under GPLv3.0 or later
## Copyright (c) Upinder S. Bhalla
##

import os
import findSim
import numpy as np
import pandas
import matplotlib.pyplot as plt

def doPlot( dataset, idx, panel, freq, doProb = False ):
    cell = dataset.split("_")[0]
    exc = "Exc" if dataset.split("_")[1] == "1" else "Inh"
    if not os.path.isfile("Expts/fs_{0}_{1}_pk.json".format(dataset,freq)):
        print("Not found: Expts/fs_{0}_{1}_pk.json".format(dataset,freq))
        return
    ax = plt.subplot( 3, 2, idx )
    ax.text( -0.25, 1.05, panel, fontsize = 22, weight = "bold", 
        transform=ax.transAxes )
    if doProb:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            "Expts/fs_{}_prob.json".format( dataset ),
            modelFile = "Models/{}_vclamp{}.py".format( exc, cell ), 
            chemFile = "Models/BothPresyn86.g",
            mapFile = "Maps/mapPresyn{}.json".format( exc ),
            bigFont = True, labelPos = "upper left", deferPlot = True )
        ex = np.array( diagnostics["exptX"] )
        ey = np.array( diagnostics["exptY"] )
        ax.plot( ex, ey, "bo-" )
    else:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            "Expts/fs_{0}_{1}_pk.json".format( dataset, freq), 
            modelFile = "Models/{}_vclamp{}.py".format( exc, cell ), 
            chemFile = "Models/BothPresyn86.g", 
            mapFile = "Maps/mapPresyn{}.json".format( exc ),
            bigFont = True, labelPos = "upper left", deferPlot = True )
        ex = np.array( diagnostics["exptX"] ) - 0.008
        ey = np.array( diagnostics["exptY"] )
        upx = ex[::2]
        dnx = ex[1::2]
        upy = ey[::2]
        dny = ey[1::2]
        ax.plot( upx, upy, "bo-" )
        ax.plot( dnx, dny, "bo-" )
        ax.text( 0.70, 0.85, "{} Hz".format( freq ), fontsize = 16, transform=ax.transAxes )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.lines[2].set_visible( False )
    ax.lines[1].set_visible( False )
    ax.lines[3].set_visible( False )
    if doProb:
        ax.set_xlim( 0.95, 1.05 )
        #ax.set_ylim( -1, 41 )
    else:
        ax.set_xlim( 0.95, 1.7 )
        ax.set_ylim( -0.05, 1.5 )
    plt.title("")

def doScoreHisto( ax ):
    bins = np.array( np.arange( 0, 0.7501, 0.05 ) )
    df = pandas.read_hdf( "model_params.h5" )
    #df.columns == ['cell', 'exc', 'numSq', 'pattern', 'initScore', 
    #'finalScore', 'Ca_bind_RR.Kd', 'Ca_bind_RR.tau', 'docking.Kf', 
    # 'vesicle_release.Kf', 'remove.Kf', 'replenish_vesicle.tau', 
    #'vesicle_pool.concInit', 'ligand_binding.tau', 'ligand_binding.Kd'],
    values = df.loc[(df['exc'] == 1) & (df['numSq'] == 5)]['finalScore']
    ax.hist(values, bins=bins, alpha=0.5, histtype = "step", 
            linewidth = 3, color = None, label = "glu-R" )
    values = df.loc[(df['exc'] == 0) & (df['numSq'] == 5)]['finalScore']
    ax.hist(values, bins=bins, alpha=0.5, histtype = "step", 
            linewidth = 3, color = None, label = "GABA-R" )
    ax.set_xlabel('Optimized model score', fontsize = 16)
    ax.set_ylabel('Frequency', fontsize = 16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend( fontsize = 14, frameon = False )
    ax.text( -0.25, 1.05, "G", fontsize = 22, weight = "bold",
        transform=ax.transAxes )


def main():
    cell = "7492"
    exc = "1"
    pat = "5"
    dataset = "{}_{}_{}".format( cell, exc, pat )
    fig = plt.figure( figsize = (10, 12), facecolor = "white" )
    fig.suptitle( dataset )
    doPlot( dataset, 1, "B", 20 )
    doPlot( dataset, 2, "C", 30 )
    doPlot( dataset, 3, "D", 40 )
    doPlot( dataset, 4, "E", 50 )
    doPlot( dataset, 5, "F", 50, doProb = True )
    ax = plt.subplot( 3,2,6 )
    doScoreHisto(ax)
    plt.show()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
