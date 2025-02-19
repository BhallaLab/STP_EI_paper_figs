import sys
import os
import numpy as np
import pandas
import argparse
import multiprocessing
import hoss

sampleRate = 20000.0
fracVesicleRelease = 1.0  # frac of vesicle pool released in first pulse.
presynModelNumber = 79
hossParams = '"Buffer.concInit", "Ca_bind_buffer.tau", "remove_Ca.Kd", "remove_Ca.tau", "Ca_bind_RR.Kd", "Ca_bind_RR.tau", "docking.Kf", "vesicle_release.Kf", "remove.Kf", "replenish_vesicle.tau", "vesicle_pool.concInit", "undocking.Kf"'
numHossParams = len( hossParams.split(",") ) + 1
excProb = 0.25
inhProb = 0.5

class EvalFunc:
    def __init__( self, cell, celldf, exc, tau, numSq, freqs, args ):
        self.cell = int( round( cell ) )
        self.celldf = celldf
        self.exc = int( round( exc ) )
        self.tau = tau  # Time course of vclamp for EPSC in seconds.
        self.numSq = int( round( numSq ) )
        self.freqs = [int( round ( ff ) ) for ff in freqs]
        self.pkDelay = args.peakDelay
        self.presynModelNumber = args.presynModelNumber
        #self.matchMinima = args.matchMinima
        self.algorithm = args.algorithm
        #self.doPairedPulse = args.pairedPulse
        self.doPairedPulse = True

    def fitFreq( self, freq ):
        pkNames = ['pk{}'.format(ii) for ii in range(9)]
        tpkNames = ['tpk{}'.format(ii) for ii in range(9)]
        valNames = ['val{}'.format(ii) for ii in range(9)]
        tvalNames = ['tval{}'.format(ii) for ii in range(9)]
        freqdf = self.celldf.loc[
                (self.celldf['exc']==self.exc) & 
                (self.celldf['numSq']==self.numSq) &
                (self.celldf['freq']==freq)
            ]
        if freqdf.empty:
            return None, None
        # Average over all patterns and all repeats
        meanVal = freqdf[valNames].mean(axis=0)
        meanTval = freqdf[tvalNames].mean(axis=0)
        meanPk = freqdf[pkNames].mean(axis=0)
        meanTpk = freqdf[tpkNames].mean(axis=0)
        means = np.array( [meanTval, meanVal, meanTpk, meanPk ] ).transpose()
        self.dumpFindSim( means, freq, False ) # for the pk/val
        if self.doPairedPulse:
            self.dumpFindSim( means, freq, True ) # for the paired pulse

    def fitCell( self ):
        if len( self.freqs) == 0:
            return None, None
        curr_proc=multiprocessing.current_process()
        # uncomment following line to get this to work
        curr_proc.daemon=False
        #stpNames = ['p0', 'p1', 'p2','p3', 'p4','p5', 'p6','p7']
        for ff in self.freqs:
            self.fitFreq( ff )
        self.dumpProbFindSim()
        self.dumpHoss()
        excLabel = "Exc" if self.exc else "Inh"
        self.dumpVclampModelDefFile("Models/" + excLabel + "_vclamp.py")
        fileSig = "{}_{}_{}".format( self.cell, self.exc, self.numSq )
        fname = "Configs/run{}.json".format( fileSig )
        hoss.main( [fname, "--algorithm", self.algorithm ] )
        df = pandas.read_table( "Results/_opt{}result.txt".format( fileSig ), sep='\s+', skiprows=3, nrows=numHossParams, engine="python", header = None )
        #assert(os.path.exists( "Expts/fs_{}_20_pk.json".format(fileSig) ) )
        paramNames = ["cell", "exc", "numSq", "freq", "initScore", "finalScore"] + list(df.iloc[:,0])
        initScore = 0.0
        finalScore = 0.0
        with open( "Results/_opt{}result.txt".format(fileSig) , "r") as fp:
            last_line = fp.readlines()[-1]
            sp = last_line.replace( ",", " ").split( " " )
            initScore = float(sp[3])
            finalScore = float( sp[-1])
        return paramNames, [self.cell, self.exc, self.numSq, 20, initScore, finalScore] + list( df.iloc[:, 2] )
        '''
        df = pandas.read_table( "Results/_opt{}results.txt".format( fileSig ), sep='\s+', skiprows=3, nrows=8, engine="python", header = None )
        paramNames = ["cell", "exc", "numSq", "freq"] + list(df.iloc[:,0])
        return paramNames, [self.cell, self.exc, self.numSq, freq] + list( df.iloc[:, 2] )
        '''

    def ticker( self, arg ):
        # arg should be the return value from the fitcell func
        #print( "{}_{}_{}    ".format( self.cell, self.exc, self.numSq ), end = "" )
        print( ".", end = '' )
        sys.stdout.flush()

    def dumpVclampModelDefFile(self, origFile ):
        # Specify the path to your file here. This is just an example path.
        base, ext = os.path.splitext(origFile)

        # Generate the new filename by inserting the string
        newFileName = base + str(self.cell) + ext
    
        # Open the original file to read
        with open(origFile, 'r') as file:
            fileContents = file.read()
    
        # Replace "TAU" with the string tau.
        modifiedContents = fileContents.replace("TAU", str(self.tau) )
    
        # Write the modified contents to a new file
        with open(newFileName, 'w') as newFile:
            newFile.write(modifiedContents)


    def dumpFindSim( self, means, freq, isPP ):
        fileSig = "{}_{}_{}_{}".format( self.cell, self.exc, self.numSq, freq)
        stimAmpl = 50.0
        Ca_basal = 0.08
        stimWidth = 0.002
        pkDelay = self.pkDelay * 0.001
        #windowStartt = -0.010
        #windowEndt = 0.015
        windowStartt = -0.015 if self.exc else -0.010
        windowEndt = 0.010 if self.exc else 0.015
        pd = [pkDelay] * 9
        label = "Exc" if self.exc else "Inh"
        clampPotl = -70.0 if self.exc else 0.0
        pkLabel = "pp" if isPP else "pk"
        if isPP:
            opLabel = "oscPk" if self.exc else "oscVal"
            baselineOp = "max" if self.exc else "min"
            sampling = "min" if self.exc else "max"
        else:
            opLabel = "oscVal" if self.exc else "oscPk"
            baselineOp = "max" if self.exc else "min"
            sampling = "min" if self.exc else "max"

        transmitter = "glu" if self.exc else "GABA"
        fsStartStr  = '''{{
    "FileType": "FindSim",
    "Version": "2.0",
    "Metadata": {{
        "transcriber": "Upi Bhalla",
        "organization": "NCBS",
        "email": "bhalla@ncbs.res.in",
        "source": {{
            "sourceType": "simulation",
            "authors": "U.S. Bhalla",
            "descriptor": "Generated FindSim file from tab_presyn_patterns_vclamp2.py"
        }},
        "testModel": "./Models/{0}_vclamp{1}.py",
        "testMap": "./Maps/mapPresyn{0}.json"
    }},  
    "Experiment": {{
        "design": "TimeSeries",
        "species": "mouse",
        "cellType": "CA1 pyramidal neuron",
        "notes": "Data is from recording of IPSC to 5 or 15 sq optical stim by Aditya Asopa"
    }},  
    "Stimuli": [
        {{
            "timeUnits": "sec",
            "quantityUnits": "mV",
            "entity": {{"name": "soma"}},
            "field": "Vclamp",
            "data": [
                [ 0.001, {2:.3f}]
            ]
        }},
        {{
            "timeUnits": "sec",
            "quantityUnits": "uM",
            "entity": {{"name": "CaInput"}},
            "field": "conc",
            "data": [
                [ 0, 0.08]'''

        fsMidStr = '''
            ]
        }}
    ],  
    "Readouts": {{
        "timeUnits": "sec",
        "quantityUnits": "ratio",
        "entity": {{"name": "vclamp"}},
        "field": "current",
        "window": {{"startt": {0}, "endt": {1}, "dt": 0.001,
            "operation": "{2}", "baseline": "{4}"
        }},
        "normalization": {{
            "entity": {{"name": "vclamp"}},
            "sampling": "{3}"
        }},
        "data": ['''.format( windowStartt, windowEndt, opLabel, sampling, baselineOp )

        fsEndStr = '\n        ]\n  }\n}\n'

        fname = "Expts/fs_{}_{}.json".format( fileSig, pkLabel )
        settleTime = 1.0
        with open( fname, "w" ) as fp:
            ############# Write the Ca stim sequence ################
            fp.write( fsStartStr.format( label, self.cell, clampPotl ) )
            for idx in range( len( means ) ):
                if isPP and idx > 2 :   # Only do first 3 points: ref, p0,p1
                    continue
                t = settleTime if idx == 0 else settleTime + 0.3 + (idx-1)/freq
                fp.write( ",\n                [{:.4f}, {:.3f}],\n".format( t, stimAmpl ) )
                fp.write( "                [{:.4f}, {:.3f}]".format( stimWidth + t, Ca_basal ) )
            fp.write( fsMidStr )

            ############# Now on to the expt/sim output ################
            comma = ""
            #for idx, dd in enumerate( deconv ):
            #print( "LENS = ",  len( meanSTP ), len( meanTpk ) )
            #print( meanSTP )
            val1 = means[1][1]
            if isPP:
                refPk = means[2][3]-val1
                for idx in [1,2]:
                    [tval, val, tpk, pk] = means[idx]
                    fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".
                        format ( comma, tval + 0.8, 
                        fracVesicleRelease * (val-val1)/refPk ) )
                    comma = ","
                    fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".
                        format( comma, tpk + 0.8, 
                        fracVesicleRelease * (pk-val1)/refPk ) )
            else:
                #refPk = means[0][3] # Earlier used first pk. Now overall max
                if sampling == "max":
                    refPk = np.max( means[:,3] ) - val1
                else:
                    refPk = np.min( means[:,3] ) - val1
                for idx, [tval, val, tpk, pk] in enumerate( means ):
                    ## Use max as normalization reference in the simulations.
                    if idx > 0:
                        fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".
                            format( comma, tval + 0.8, 
                            fracVesicleRelease * (val-val1)/refPk ) )
                        fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".
                            format( comma, tpk + 0.8, 
                            fracVesicleRelease * (pk-val1)/refPk ) )
                    else: # First pk is always one.
                        fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".
                            format( comma, tpk + 0.8, 
                            fracVesicleRelease * (pk-val1)/refPk ) )
                    comma = ","
            fp.write( fsEndStr )


    def dumpProbFindSim( self ):
        fileSig = "{}_{}_{}".format( self.cell, self.exc, self.numSq )
        stimAmpl = 50.0
        Ca_basal = 0.08
        stimWidth = 0.002
        pkDelay = self.pkDelay * 0.001
        windowStartt = 0.000
        windowEndt = 0.010
        pd = [pkDelay] * 9
        if self.exc:
            label = "Exc" if self.exc else "Inh"
            value = -np.log(1-excProb)/(windowEndt - windowStartt) 
        else:
            label = "Inh"
            value = -np.log(1-inhProb)/(windowEndt - windowStartt)

        transmitter = "glu" if self.exc else "GABA"
        fsStr  = '''{{
    "FileType": "FindSim",
    "Version": "1.0",
    "Metadata": {{
        "transcriber": "Upi Bhalla",
        "organization": "NCBS",
        "email": "bhalla@ncbs.res.in",
        "source": {{
            "sourceType": "simulation",
            "authors": "U.S. Bhalla",
            "descriptor": "Generated FindSim file from tab_presyn_patterns_vclamp2.py"
        }},
        "testModel": "./Models/{0}_vclamp{1}.py",
        "testMap": "./Maps/mapPresyn{0}.json"
    }},  
    "Experiment": {{
        "design": "TimeSeries",
        "species": "mouse",
        "cellType": "CA1 pyramidal neuron",
        "notes": "Data is from recording of IPSC to 5 or 15 sq optical stim by Aditya Asopa"
    }},  
    "Stimuli": [
        {{
            "timeUnits": "sec",
            "quantityUnits": "uM",
            "entity": {{"name": "CaInput"}},
            "field": "conc",
            "data": [
                [ 0, 0.08],
                [1.0100, 50.000],
                [1.0120, 0.080],
                [1.0400, 0.080]
            ]
        }}
    ],  
    "Readouts": {{
        "timeUnits": "sec",
        "quantityUnits": "#",
        "entity": {{"name": "{2}"}},
        "field": "n",
        "window": {{"startt": {3}, "endt": {4}, "dt": 0.001,
            "operation": "mean"
        }},
        "data": [
            [1.0, 0, 0],
            [1.011, {5}, 0],
            [1.020, {5}, 0],
            [1.030, {5}, 0]
        ]
    }},
    "Modifications": {{
        "itemsToDelete": [ {{"name": "remove"}}]
    }}
}}'''.format( label, self.cell, transmitter, windowStartt, windowEndt, value )

        fname = "Expts/fs_{}_prob.json".format( fileSig )
        with open( fname, "w" ) as fp:
            ############# Write the Ca stim sequence ################
            fp.write( fsStr )

    def dumpHoss( self ):
        excLabel = "Exc" if self.exc else "Inh"
        fileSig = "{}_{}_{}".format( self.cell, self.exc, self.numSq )
        fname = "Configs/run{}.json".format( fileSig )
        exptStr = ''
        for iff in self.freqs:
            exptStr += '"fs_{}_{}_pk.json": {{"weight": {} }},\n'.format(
                fileSig, iff, 100 )
            if self.doPairedPulse:
                exptStr+='"fs_{}_{}_pp.json": {{"weight": {} }},\n'.format(
                    fileSig, iff, 200 )
        exptStr += '"fs_{}_prob.json": {{"weight": 400 }}\n'.format(fileSig)
    
        startStr = '''{{
    "FileType":"Hoss",
    "Version":"2.0",
    "author": "TabPresyn program",
    "model": "Models/{0}_vclamp{1}.py",
    "map":"Maps/mapPresyn{0}.json",
    "exptDir": "./Expts",
    "outputDir": "./Results",
    "scoreFunc": "NRMS",
    "tolerance": 0.0001,
    "algorithm": "COBYLA",
    "comment": "Program-generated HOSS config for presynaptic release.",
    "hossMethod": {{
        "method": "hoss",
        "numProcesses": 1
    }},
        "HOSS": [
            {{
            "name": "IndividualPathways", 
            "hierarchyLevel": 1,
            "presyn": {{
                "comment": "This is the only block at present",
                "expt": {{{2}
                }},
                "params": [{3}],
                "paramBounds": {{ "vesicle_pool.concInit":[0.3e-4,10.0e-3,0],
                "docking.Kd":[1e-14, 0.00001,0]
                }},

                "resultFile": "_opt{4}result.txt",
                "optModelFile": "opt{4}.g"
            }}
        }}
    ]\n}}'''.format( excLabel, self.cell, exptStr, hossParams, fileSig)

        with open( fname, "w" ) as fp:
            fp.write( startStr )

def main():
    global fracVesicleRelease
    parser = argparse.ArgumentParser( description = "Fits presynaptic release model to STP curves from pandas file. Required directories: Configs, Models, Expts, Results, Maps" )
    parser.add_argument( "-o", "--output", type = str, help = "Optional: output pandas hdf file for model params, default = model_params.h5", default = "model_params.h5")
    parser.add_argument( "-f", "--file", type = str, 
            help = "Optional: Name of tabulated STP file in pandas hdf5 format. Default = 'STP_pks_and_refs.h5'", 
            default = "STP_pks_and_refs.h5" )
    parser.add_argument( "-pd", "--peakDelay", type = float, 
            help = "Optional: Delay for time of glu peak in ms. Default = 10", default = 10.0)
    parser.add_argument( '-n', '--numProcesses', type = int, default = 0, help='Optional: Number of processes to spawn' )
    parser.add_argument( '-pmn', '--presynModelNumber', type = int, default = 76, help='Optional: Presynaptic model number' )
    parser.add_argument( '-fvr', '--fracVesicleRelease', type = float, default = 1.0, help='Optional: Fraction of total vesicle count released on first pulse. default 1.0' )
    parser.add_argument( '-pp', '--pairedPulse', action = "store_true", help='Flag: Turns on dumping of experiments to do paired pulse matching' )
    parser.add_argument( '-a', '--algorithm', type = str, help='Optional: Algorithm name to use, from the set available to scipy.optimize.minimize. Options are CG, Nelder-Mead, Powell, BFGS, COBYLA, SLSQP, TNC. The library has other algorithms but they either require Jacobians or they fail outright. There is also L-BFGS-B which handles bounded solutions, but this is not needed here because we already take care of bounds. COBYLA works well and is the default.', default = "COBYLA" )

    args = parser.parse_args()
    dat = pandas.read_hdf( args.file )
    #cellList = dat['cell'].unique()
    cellDict = {
        # key: (gluTau, gabaTau),
        1931: (8.9, 20.0),
        1491: (10.0, 16.0),
        1541: (12.0, 16.0),
        1524: (15.0, 24.0),
        #1523: (14.5, None),
        1522: (15.0, 24.0),
        1531: (11.6, 16.0),
        1621: (15.0, 16.0),
        111: (20.0, 16.0),
        7492: (20, 24.0),
        7491: (20.0, 24.0),
        6201: (10.0, 20.0),
        6301: (20.0, 24.0),
        #5501: (1701.9, None)

    }
    #print( "Cell list = ", cellList, len( cellList ) )

    numProcesses = args.numProcesses

    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count() // 4
    if numProcesses == 0:
        numProcesses = 1
    numProcesses = min( numProcesses, len( cellDict ) * 4 )
    fracVesicleRelease = args.fracVesicleRelease

    print( "NUM PROCESSES = ", numProcesses )
    pool = multiprocessing.Pool( processes = numProcesses )
    params = []
    ret = []
    for cell in cellDict:
        for exc in [ 0, 1 ]:
            for numSq in [ 5, 15 ]:
                print( "cell={} exc={} numSq={}".format( int(cell), exc, numSq ))
                celldf = dat.loc[dat['cell']==cell]
                freqs = celldf.loc[celldf['numSq']==numSq]['freq'].unique()
                ev = EvalFunc( cell, celldf, exc, cellDict[cell][1-exc]*1e-3, numSq, freqs, args )
                ret.append( pool.apply_async( ev.fitCell, callback=ev.ticker ))
    ans = [ rr.get() for rr in ret ]

    for aa in ans:
        [paramNames, paramVals] = aa
        if paramNames:
            params.append( paramVals )
    df = pandas.DataFrame( params, columns = paramNames )
    df.info()
    df.to_hdf( args.output, "w" )

    #pool.join()
    pool.close()

if __name__ == "__main__":
    main()
