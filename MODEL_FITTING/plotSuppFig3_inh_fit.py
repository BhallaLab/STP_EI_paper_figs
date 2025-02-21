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
    cell = dataset.split("_")[0]
    exc = "Exc" if dataset.split("_")[1] == "1" else "Inh"
    if not os.path.isfile("Expts/fs_{0}_{1}_pk.json".format(dataset,freq)):
        print( "Could not find ", "Expts/fs_{0}_{1}_pk.json".format(dataset,freq))
        return
    ax = plt.subplot( 3, 2, idx )
    ax.text( -0.10, 1.05, panel, fontsize = 22, weight = "bold", 
        transform=ax.transAxes )
    if doProb:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            #"ExptFilesForFigure/fs_7492_1_5_49_{}_pk.json".format( freq ), 
            #modelFile = "ResultsForFigure/opt7492_1_5_49.g", 
            "Expts/fs_{}_prob.json".format( dataset, freq ),
            modelFile = "Models/{}_vclamp{}.py".format( exc, cell ), 
            chemFile = "Models/BothPresyn86.g",
            mapFile = "Maps/mapPresyn{}.json".format( exc ),
            bigFont = True, labelPos = "upper left", deferPlot = True )
        ex = np.array( diagnostics["exptX"] )
        ey = np.array( diagnostics["exptY"] )
        ax.plot( ex, ey, "bo-" )
        #upx = ex[::2]
        #dnx = ex[1::2]
        #upy = ey[::2]
        #dny = ey[1::2]
    else:
        score, elapsedTime, diagnostics = findSim.innerMain( 
            #"ExptFilesForFigure/fs_7492_1_5_49_{}_pk.json".format( freq ), 
            #modelFile = "ResultsForFigure/opt7492_1_5_49.g", 
            "Expts/fs_{0}_{1}_pk.json".format( dataset, freq), 
            modelFile = "Models/{}_vclamp{}.py".format( exc, cell ), 
            chemFile = "Models/BothPresyn86.g",
            mapFile = "Maps/mapPresyn{}.json".format( exc ),
            bigFont = True, labelPos = "upper left", deferPlot = True )
        ex = np.array( diagnostics["exptX"] )
        ey = np.array( diagnostics["exptY"] )
        upx = ex[::2]
        dnx = ex[1::2]
        upy = ey[::2]
        dny = ey[1::2]
        ax.plot( upx, upy, "bo-" )
        ax.plot( dnx, dny, "bo-" )
    #ax.lines[3].set_visible( False )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.lines[2].set_visible( False )
    ax.lines[1].set_visible( False )
    ######ax.lines[0].set_visible( False )
    #ax.lines[3].set_visible( False )
    #ax.lines[1].set_linestyle( None )
    if doProb:
        ax.set_xlim( 0.95, 1.05 )
        #ax.set_ylim( -1, 41 )
    else:
        ax.set_xlim( 0.95, 1.7 )
        ax.set_ylim( -0.05, 1.7 )
    plt.title("")

def doScoreHisto( ax ):
    #bins = "auto"
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
    ax.text( -0.10, 1.05, "F", fontsize = 22, weight = "bold",
        transform=ax.transAxes )


def main():
    for cell in ["7492"]:
    #for cell in ["1931", "1524", "1522", "1621", "111", "7491"]:
    #for cell in ["7492", "1621", "1491", "6301", "6201", "1541", "1531"]:
        #for exc in ["0", "1"]:
        for exc in ["0",]:
            for pat in [ "15"]: 

                dataset = "{}_{}_{}".format( cell, exc, pat )
                fig = plt.figure( figsize = (10, 12), facecolor = "white" )
                fig.suptitle( dataset )
                doPlot( dataset, 1, "A", 20 )
                if cell not in ["7491", "1931"]:
                    doPlot( dataset, 2, "B", 30 )
                    doPlot( dataset, 3, "C", 40 )
                doPlot( dataset, 4, "D", 50 )
                doPlot( dataset, 5, "E", 20, doProb = True )
                ax = plt.subplot( 3,2,6 )
                doScoreHisto(ax)
                plt.show()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
