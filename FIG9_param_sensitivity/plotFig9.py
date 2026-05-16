import pandas as pd
import pylab
import numpy as np
import math
from scipy.stats import linregress
from scipy.stats import wilcoxon
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

SAMPLE_FREQ = 20000
SAMPLE_RATIO = 10 # Ratio between SAMPLE_FREQ and SIM_SAMPLE_FREQ
SPIKE_RATIO = 2 # Ratio between SAMPLE_FREQ and SPK_SAMPLE_FREQ
SAMPLE_TIME = 5
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
NUM_SIM_SAMPLES = NUM_SAMPLES // SAMPLE_RATIO
SIM_SAMPLE_FREQ = SAMPLE_FREQ//SAMPLE_RATIO
NUM_SPK_SAMPLES = NUM_SAMPLES // SPIKE_RATIO
SPK_SAMPLE_FREQ = SAMPLE_FREQ//SPIKE_RATIO
pkDelay = int(0.003 * (SAMPLE_FREQ // SAMPLE_RATIO))
SAMPLE_START = 49
BAR_DATA_DIR = "./BARDATA/"
SPK_DATA_DIR = "./"

paramSweepFile = "param_sweep_digest.h5"
referenceVals = { "Pattern Overlap (%)": 8.87, "wtGlu": 5, "wtGABA": 10,
        "pCA3_CA1": 0.02 ,"pCA3_Inter": 0.01 ,"pInter_CA1":0.01 }

freq = 80.0 # Hz
params = {
    "zeroIndices": [ 32, 64, 96, 128, 192, 224, 240],
    "wtGlu": [0.2, 0.5, 1, 2, 3, 5, 10, 20],
    "wtGABA": [1, 2, 5, 10, 20, 50, 100],
    "pCA3_CA1": [0.005, 0.01, 0.02, 0.05, 0.10],
    "pCA3_Inter": [0.002, 0.005, 0.01, 0.02, 0.05],
    "pInter_CA1": [0.002, 0.005, 0.01, 0.02, 0.05]
}
FREQ_SWEEP_DIR = "./FREQ_SWEEP"

RomNums = ["i", "ii", "iii", "iv", "v", "vi"]
Overlap = {32:33.72, 64:28.968, 96:24.32, 128:19.56, 192:8.87, 224:4.88, 240:3.096}

def panelB( dfs, fig, gs, row, paramList ):
    for idx, pp in enumerate( paramList ):
        dp = dfs.loc[dfs["param"]==pp]
        flist = dp["freq"].unique()
        xlab = pp
        if pp == "zeroIndices": # flip to overlap
            xlab = "Pattern Overlap (%)"
        ax = fig.add_subplot( gs[ row, idx ] )
        for ii, ff in enumerate( reversed( flist ) ):
            dfreq = dp.loc[dp["freq"]==ff].copy()
            numsig_cols = ['n_wp1', 'n_wp2', 'n_wp3']
            dfreq[numsig_cols] = dfreq[['wp1', 'wp2', 'wp3']]<0.01
            dfreq['tot_numsig'] = dfreq[numsig_cols].sum(axis=1)/3.0
            # Group by 'val' to average across the multiple cells/samples
            final_plot_data = dfreq.groupby('val')['tot_numsig'].sum().reset_index()

            x = final_plot_data['val']
            if pp == "zeroIndices": # flip to overlap
                x = [ Overlap[zz] for zz in x ]
            y = 100 * final_plot_data['tot_numsig'] / len( final_plot_data )
            colors = ["magenta", "seagreen", "blue"]
            ax.plot( x,y, color=colors[ii], markersize=5, label = str(ff)+ " Hz" )
        ax.scatter( [referenceVals[xlab]], [2.5], 
            marker = '^', color = 'red', s=100 )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel( xlab )
        if idx == 0:
            ax.set_ylabel( "% selective" )
        ax.set_ylim( 0, 64 )
        if idx == 0 and row == 1:
            #ax.legend( fontsize=16, frameon=False, facecolor='none', ncol=3)
            ax.legend( fontsize=14, frameon=False, facecolor='none', loc="upper right", ncol=2)
            ax.text( -0.25, 1.10, "B", fontsize = 22, weight = "bold", transform=ax.transAxes )
        label = RomNums[row*3+idx-3]
        ax.text( -0.15, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )

# This is tuned to higher freqs but gives very similar output as reference.
def findPeaksMedHigh( Vm, pulses, freq ):
    pkDelaySamples = int(SIM_SAMPLE_FREQ * 0.016)
    widthSamples = int( np.round( 0.012 * SIM_SAMPLE_FREQ ) )
    pks = []
    for pp in pulses:
        idx = int(pp * SIM_SAMPLE_FREQ)
        #print( "IDX=", idx, pkDelay, widthSamples, freq)
        idx2 = idx + pkDelaySamples - widthSamples
        pk = np.max( Vm[idx2:idx2 + widthSamples*2] )
        pks.append( pk )
    return pks

# This is tuned to the reference 8, 20, 50Hz stim patterns.
def findPeaks( Vm, pulses, freq ):
    pkDelaySamples = int(SIM_SAMPLE_FREQ * 0.018)
    widthSamples = int( np.round( 0.016 * SIM_SAMPLE_FREQ ) )
    pks = []
    for pp in pulses:
        idx = int(pp * SIM_SAMPLE_FREQ)
        #print( "IDX=", idx, pkDelay, widthSamples, freq)
        idx2 = idx + pkDelaySamples - widthSamples
        pk = np.max( Vm[idx2:idx2 + widthSamples*2] )
        pks.append( pk )
    return pks

# This uses the tuning parameters from Sep 2025. Doesn't work well.
def findPeaks2025( Vm, pulses, freq, width = 0.001 ):
    if len( Vm ) < NUM_SIM_SAMPLES or len( pulses ) != 33:
        print( "BAD" )
        return []
    widthSamples = int( np.round( width * SIM_SAMPLE_FREQ ) )
    pd = int( 0.0056 * SIM_SAMPLE_FREQ )
    pks = []
    for pp in pulses:
        idx = int(pp * SIM_SAMPLE_FREQ)
        idx1 = idx + pd - widthSamples
        idx2 = idx + pd + widthSamples
        #print( idx, idx2, len( Vm ) )
        pk = max( Vm[idx1:idx2])
        pks.append( pk )
        #print( "{}  {:.3f}  {:.3f}".format( idx, val, pk ) )
    return pks


def findBlockValPk( pks ):
    x = [ 8, 16, 24 ]
    pre = [ np.mean( pks[ii-3:ii] ) for ii in x ]
    post =[ np.mean( pks[ii:ii+3] ) for ii in x ]
    return np.array(x), np.array(pre), np.array(post)

def parseRow( df, cell, ff ):
    PulseTrain = {}

    # Finds the field and epsp peaks for each pulse.
    # If any are too small, it puts in a zero.
    epsp = np.array(df.iloc[0, SAMPLE_START:SAMPLE_START+NUM_SIM_SAMPLES ])
    pulseTrig = np.array(df.iloc[0, SAMPLE_START+NUM_SIM_SAMPLES: ] )
    epsp -= min(epsp) # hack to handle traces with large ipsps.
    '''
    padt = np.pad( pulseTrig, 1 )
    edges = pulseTrig
    assert( len( edges ) == 33 )
    PulseTrain[ff] = np.array( SIM_SAMPLE_FREQ * edges, dtype = int)
    pks = findPeaks( epsp, PulseTrain[ff], ff )
    '''
    pks = findPeaks( epsp, pulseTrig, ff )
    x, miny, maxy = findBlockValPk( pks )
    return maxy / miny

def scanData( df ):
    idx = 0
    cellStats = {}
    cellList = df['cellID'].unique()
    temp = df['stimFreq'].unique()
    freq5 = { ff:[] for ff in temp }
    pk5 = []
    for cellIdx, cell in enumerate( cellList ):
        dcell = df.loc[(df["cellID"] == cell)]
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
                    pkRatio = parseRow( dseq, cell, ff )
                    pk5.append( pkRatio )
                    freq5[ff].append( pkRatio )

    return pk5, freq5

def panelC(fig, gs):
    pidx = 0
    for parmname, pvals in params.items():
        maxy = []
        for vv in pvals:
            df = pd.read_hdf( f"{FREQ_SWEEP_DIR}/fparm_{parmname}_{vv}.h5" )
            for nidx, nn in enumerate( df['numSq'].unique() ):
                ndf = df[df['numSq']==nn]
                pk5, freq5 = scanData( ndf )
                freqlist = sorted(list(freq5.keys()))
                y = [np.mean(freq5[ff]) for ff in freq5 ]
                maxy.append(freqlist[np.argmax(np.array(y))])
        x = list( pvals )
        xlab = parmname
        if parmname == "zeroIndices": # flip to overlap
            xlab = "Pattern Overlap (%)"
            x = [ Overlap[pp] for pp in pvals ]
        ax = fig.add_subplot(gs[3 + pidx//3, pidx % 3])
        ax.plot( x, maxy )
        ax.set_ylim( 0, 210 )
        #ax.set_xscale( 'log' )
        ax.scatter( [referenceVals[xlab]], [10], 
            marker = '^', color = 'red', s=100 )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel( xlab, fontsize = 16 )
        if pidx in [0,3]:
            ax.set_ylabel( "Peak Freq (Hz)", fontsize = 16 )
        if pidx == 0:
            ax.text( -0.25, 1.10, "C", fontsize = 22, weight = "bold", transform=ax.transAxes )
        label = RomNums[pidx]
        ax.text( -0.15, 1.10, label, fontsize = 22, weight = "bold", transform=ax.transAxes )
        pidx += 1

def processPanelD(fname):
    """Loads file and calculates the 3-block peak ratios."""
    file_path = BAR_DATA_DIR + fname
    ret = []
    #plt.figure()
    try:
        df = pd.read_hdf(file_path)
        for numSq in [5,15]:
            sf = df[df['numSq'] == numSq]
            for freq in [8, 20, 50]:
                fdf = sf[sf['stimFreq'] == freq]
                for cellID in range( 1, 7 ):
                    cdf = fdf[fdf['cellID']==cellID]
                    epsp = np.array(cdf.iloc[:, SAMPLE_START : SAMPLE_START + NUM_SIM_SAMPLES]).mean( axis = 0 )
                    pulses = np.array(cdf.iloc[0, SAMPLE_START + NUM_SIM_SAMPLES : ])
                    #plt.plot( epsp )
                    #plt.show()
                    #quit()
                    epsp -= min(epsp)
                    pks = findPeaks(epsp, pulses, freq)
                    block, miny, maxy = findBlockValPk(pks)
                    bvp = (maxy - miny)/(maxy+miny)
                    for bb in bvp:
                        ret.append({
                            "fname": fname,
                            "cellID": cellID,
                            "numSq": numSq,
                            "freq": freq,
                            "val": bb,
                        })
        return pd.DataFrame(ret)
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return []

def spkProcessPanelD(fname):
    """Loads file and calculates the 3-block peak ratios."""
    file_path = BAR_DATA_DIR + fname
    ret = []
    try:
        df = pd.read_hdf(file_path)
        #for numSq in [5,15]:
        print( "CELLID = ", df['cellID'].unique() )
        print( "freq = ", df['stimFreq'].unique() )
        print( "sweep = ", df['sweep'].unique() )
        for numSq in [5,15]:
            sf = df[df['numSq'] == numSq]
            for freq in [8, 20, 50]:
                valPreIdx = int(3*SPK_SAMPLE_FREQ/freq)
                pkPreIdx = 0
                pkPostIdx = int(3*SPK_SAMPLE_FREQ/freq)
                fdf = sf[sf['stimFreq'] == freq]
                for cellID in range( 1, 7 ):
                    miny = []
                    maxy = []
                    cdf = fdf[fdf['cellID']==cellID]
                    epsp = np.array(cdf.iloc[:, SAMPLE_START : SAMPLE_START + NUM_SPK_SAMPLES])
                    padEpsp = np.pad( epsp, pad_width = ( (0,0), (0,1)), mode='constant', constant_values = -60 )
                    spks = (epsp < -30) & (padEpsp[:,1:] > -30)
                    pulses = np.array(cdf.iloc[0, SAMPLE_START + NUM_SPK_SAMPLES : ])
                    pulse_idx = ( pulses * SPK_SAMPLE_FREQ ).astype(int)
                    for ii in [pulse_idx[8], pulse_idx[16], pulse_idx[24]]:
                        miny.append( sum( spks[:,ii - valPreIdx:ii].flatten() ) )
                        maxy.append(sum( spks[:,ii + pkPreIdx:ii+pkPostIdx].flatten() ) )
                    print( f"{file_path}.{freq}: miny = {[int(yy) for yy in miny]}, maxy = {[int(yy) for yy in maxy]}" )

                    bvp = (np.array(maxy) - np.array(miny))/(np.array(maxy)+np.array(miny))
                    for bb in bvp:
                        ret.append({
                            "fname": fname,
                            "cellID": cellID,
                            "numSq": numSq,
                            "freq": freq,
                            "val": bb,
                        })
        return pd.DataFrame(ret)
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return []

def statsPanelD(df, fileList):
    """Generates stats comparing each group with the baseline (determ_orig_0.h5)."""
    ref_name = "determ_orig_0.h5"
    ref_data = df[df["fname"] == ref_name]["val"]
    groups = [f for f in fileList if f != ref_name]
    p_values = []
    
    print("\n--- Statistical Results (Welch's T-Test) ---")
    for g in groups:
        group_data = df[df["fname"] == g]["val"]
        # Welch's t-test (equal_var=False)
        t_stat, p_val = ttest_ind(group_data, ref_data, equal_var=False)
        p_values.append(p_val)
        print(f"Group: {g:20s} | T-stat: {t_stat:.4f} | P-value: {p_val:.4e}")
        
    return p_values

def get_sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

def plotPanelD(fig, gs, df, fileList, p_values): 
    ax = fig.add_subplot(gs[5:7,:])
    # --- Reference ---
    ref_name = "determ_orig_0.h5"
    ref = df[df["fname"] == ref_name].set_index("freq")
    
    # --- Other groups ---
    labelList = ["Stochasticity", "no STP", "Physiol Temp", "Spiking" ]
    groups = [f for f in fileList if f != ref_name]
    freqs = [8, 20, 50]
    colors = {50:'#C000C0', 20:'lightseagreen', 8:'blue'}
    
    # --- Plot setup ---
    bar_width = 0.25
    x = np.arange(len(groups))
    
    # To keep track of max height for asterisk placement
    all_heights_with_err = []
    
    for i, freq in enumerate(freqs):
        heights = []
        errors = []
        
        for g in groups:
            val = df[(df["fname"]==g) & (df["freq"]==freq)]["val"]
            mean_val = val.mean( axis = 0 )
            sem_val = val.std( axis = 0 ) / np.sqrt( len(val ) )
            
            ref_mean = ref.loc[freq]["val"].mean( axis = 0 )
            ref_sem = ref.loc[freq]["val"].std( axis = 0 ) / np.sqrt( len( val ) )
            
            # Percent difference
            pct = 100*(mean_val - ref_mean)/ref_mean
            
            # Error propagation for (a-b)/b: avoids dividing by mean_val
            print( f"{g:20s}: mean_val={mean_val:.3f}, ref_mean={ref_mean:.3f}")
            pct_sem = 100 * np.sqrt((sem_val/ref_mean)**2 + (mean_val*ref_sem/ref_mean**2)**2)
            
            heights.append(pct)
            errors.append(pct_sem)
            all_heights_with_err.append(pct + pct_sem)
        
        ax.bar(x + i*bar_width,
                heights,
                bar_width,
                yerr=errors,
                color = colors[freq],
                capsize=4,
                label=f"{freq} Hz")
    
    # Plot Significance Asterisks at the same level
    y_level = max(all_heights_with_err) + 5 if all_heights_with_err else 10
    for i, p in enumerate(p_values):
        ax.text(x[i] + bar_width, y_level, get_sig_label(p), 
                ha='center', va='bottom', fontsize=14, fontweight='bold')

    # Formatting
    ax.axhline(0, color='black', linewidth=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(labelList, fontsize=16)
    ax.set_ylabel("% change", fontsize = 16)
    ax.set_ylim( -200, 250 )
    #ax.set_title("Analysis of Peak Differences across Blocks")
    ax.legend( title="Frequency", fontsize=14, title_fontsize=14,
        frameon=False, loc="upper left" )
    ax.text( -0.072, 1.10, "D", fontsize = 22, weight = "bold", transform=ax.transAxes )

def panelD( fig, gs ):
    fileList = ["determ_orig_0.h5", "stoch_orig_0.h5","noSTP_orig_0.h5", "temp37deg_orig_0.h5", "spk_orig_0.h5" ]
    
    results = []
    for fname in fileList:
        #if fname == "spk_orig_0.h5":
        if "spk" in fname:
            tempdf = spkProcessPanelD( fname )
        else:
            tempdf = processPanelD( fname )
        results.append( tempdf )
    df = pd.concat(results, ignore_index=True)

    df_filtered = df[df['numSq']==15]
    
    # 1. Generate stats and get p-values
    p_values = statsPanelD(df_filtered, fileList)
    
    # 2. Plot using the p-values for annotation
    plotPanelD(fig, gs, df_filtered, fileList, p_values)

def main():
    plt.rcParams.update( {"font.size": 16} )
    fig = plt.figure( figsize = (15,21) )
    #fig.suptitle( "Fig8_v18", fontsize = 16 )
    gs = fig.add_gridspec( 7, 3 ) # 7 rows, 3 cols, but top row is empty
    # Row 1: Leave blank, for the schematics.
    dfs = pd.read_hdf( paramSweepFile )
    # Row 2,3: Param variations: zeroIndices, wtGlu, wtGABA
    panelB( dfs, fig, gs, 1, ["zeroIndices", "wtGlu", "wtGABA"] )
    panelB( dfs, fig, gs, 2, ["pCA3_CA1","pCA3_Inter","pInter_CA1"])
    panelC( fig, gs ) # Row 3,4: Freq peak param variations
    panelD( fig, gs ) # Row 5: Barcharts of effs on response.

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
