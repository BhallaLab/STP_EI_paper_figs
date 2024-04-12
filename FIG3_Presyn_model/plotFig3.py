import os
import findSim
import numpy as np
import pandas
import matplotlib.pyplot as plt


# Get rid of the frames
# Get rid of time axis except for bottom plot.
# Convert solid blue line into just the dots. Better yet, 
# link above and below
def doPlot( dataset, idx, panel, freq, doProb = False ):
    if not os.path.isfile("Expts/fs_{0}_{1}_pk.json".format(dataset,freq)):
        return
    ax = plt.subplot( 3, 2, idx )
    ax.text( -0.20, 1.05, panel, fontsize = 22, weight = "bold", 
        transform=ax.transAxes )
    if doProb:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            #"ExptFilesForFigure/fs_7492_1_5_49_{}_pk.json".format( freq ), 
            #modelFile = "ResultsForFigure/opt7492_1_5_49.g", 
            "Expts/fs_{}_prob.json".format( dataset ),
            modelFile = "Results/opt{}.g".format( dataset ), 
            mapFile = "Maps/mapPresyn.json",
            bigFont = True, labelPos = "upper left", deferPlot = True )
    else:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            #"ExptFilesForFigure/fs_7492_1_5_49_{}_pk.json".format( freq ), 
            #modelFile = "ResultsForFigure/opt7492_1_5_49.g", 
            "Expts/fs_{0}_{1}_pk.json".format( dataset, freq), 
            modelFile = "Results/opt{}.g".format( dataset ), 
            mapFile = "Maps/mapPresyn.json",
            bigFont = True, labelPos = "upper left", deferPlot = True )
    #ax.lines[3].set_visible( False )
    ex = np.array( diagnostics["exptX"] )
    ey = np.array( diagnostics["exptY"] )
    upx = ex[::2]
    dnx = ex[1::2]
    upy = ey[::2]
    dny = ey[1::2]
    ax.plot( upx, upy, "bo-" )
    ax.plot( dnx, dny, "bo-" )
    ax.lines[2].set_visible( False )
    ax.lines[1].set_visible( False )
    #ax.lines[0].set_visible( False )
    ax.lines[3].set_visible( False )
    #ax.lines[1].set_linestyle( None )
    if doProb:
        ax.set_xlim( 0.95, 1.05 )
        #ax.set_ylim( -1, 41 )
    else:
        ax.set_xlim( 0.95, 1.7 )
        ax.set_ylim( -0.05, 1.6 )
    plt.title("")

def doScoreHisto( ax ):
    df = pandas.read_hdf( "model_params.h5" )
    #df.columns == ['cell', 'exc', 'numSq', 'pattern', 'initScore', 
    #'finalScore', 'Ca_bind_RR.Kd', 'Ca_bind_RR.tau', 'docking.Kf', 
    # 'vesicle_release.Kf', 'remove.Kf', 'replenish_vesicle.tau', 
    #'vesicle_pool.concInit', 'ligand_binding.tau', 'ligand_binding.Kd'],
    values = df.loc[(df['exc'] == 1) & (df['numSq'] == 5)]['finalScore']
    ax.hist(values, bins=10, alpha=0.5, histtype = "step", 
            linewidth = 3, color = None, label = "Exc" )
    values = df.loc[(df['exc'] == 0) & (df['numSq'] == 5)]['finalScore']
    ax.hist(values, bins=10, alpha=0.5, histtype = "step", 
            linewidth = 3, color = None, label = "Inh" )
    ax.set_xlabel('Optimized model score', fontsize = 16)
    ax.set_ylabel('Frequency', fontsize = 16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.legend()
    ax.text( -0.10, 1.05, "F", fontsize = 22, weight = "bold",
        transform=ax.transAxes )


def main():
    #for cell in ["7492", "7491", "1621"]:
    for cell in ["7492", "1491", "6301", "6201", "1541", "1531"]:
        #for exc in ["0", "1"]:
        for exc in ["0",]:
            for pat in [ "5_46", "5_47", "5_48", "5_49", "5_50", "15_52", "15_53", "15_55"]: 

                dataset = "{}_{}_{}".format( cell, exc, pat )
                fig = plt.figure( figsize = (10, 10), facecolor = "white" )
                fig.suptitle( dataset )
                doPlot( dataset, 1, "B", 20 )
                if cell not in ["7491", "1931"]:
                    doPlot( dataset, 2, "C", 30 )
                    doPlot( dataset, 3, "D", 40 )
                doPlot( dataset, 4, "E", 50 )
                doPlot( dataset, 5, "5", 50, doProb = True )
                ax = plt.subplot( 3,2,6 )
                doScoreHisto(ax)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
