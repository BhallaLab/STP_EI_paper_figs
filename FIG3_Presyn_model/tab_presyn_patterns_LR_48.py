import sys
import os
import numpy as np
import pandas
import argparse
import multiprocessing
import hoss

sampleRate = 20000.0
fracVesicleRelease = 1.0  # frac of vesicle pool released in first pulse.
presynModelNumber = 75
transmitterDecayTime = 5e-3 # 5 ms for sure.
transmitterBaseline = 0.0   # Approaches it
hossParams = '"Ca_bind_RR.Kd", "Ca_bind_RR.tau", "docking.Kf", "vesicle_release.Kf", "remove.Kf", "replenish_vesicle.tau", "vesicle_pool.concInit", "ligand_binding.tau", "ligand_binding.Kd"'
numHossParams = len( hossParams.split(",") ) + 1

class EvalFunc:
    def __init__( self, cell, rowList, exc, numSq, pattern, idx, args ):
        # rowList has an entry for each of the freqs, in order. idx is trial idx.
        self.cell = int( round( cell ) )
        self.rowList = rowList
        self.exc = exc
        self.numSq = numSq
        self.pattern = pattern
        self.idx = idx
        self.pkDelay = args.peakDelay
        self.presynModelNumber = args.presynModelNumber
        #self.matchMinima = args.matchMinima
        self.algorithm = args.algorithm
        #self.doPairedPulse = args.pairedPulse
        self.doPairedPulse = True

    def fitCell( self ):
        curr_proc=multiprocessing.current_process()
        # uncomment following line to get this to work
        curr_proc.daemon=False
        #stpNames = ['p0', 'p1', 'p2','p3', 'p4','p5', 'p6','p7']
        pkNames = ['pk{}'.format(ii) for ii in range(9)]
        tpkNames = ['tpk{}'.format(ii) for ii in range(9)]
        valNames = ['val{}'.format(ii) for ii in range(9)]
        tvalNames = ['tval{}'.format(ii) for ii in range(9)]
        for rr in self.rowList:
            iff = int( round( rr['freq'] ) )
            '''
            if self.cell == 7492 and self.exc == 0 and self.numSq == 5 and self.pattern == 48:
                print( "cell=", self.cell, " IFF=", iff, " idx=", self.idx )
            '''
            meanVal = rr[valNames]
            meanTval = rr[tvalNames]
            meanPk = rr[pkNames]
            meanTpk = rr[tpkNames]
            means = np.array( [meanTval, meanVal, meanTpk, meanPk ] ).transpose()
            self.dumpFindSim( means, iff, False ) # for the pk/val
            if self.doPairedPulse:
                self.dumpFindSim( means, iff, True ) # for the paired pulse
        self.dumpHoss()
        excLabel = "Exc" if self.exc else "Inh"
        fileSig = "{}_{}_{}_{}_{}".format( self.cell, self.exc, self.numSq, self.pattern, self.idx )
        fname = "TConfigs/run{}.json".format( fileSig )
        #os.sync()
        hoss.main( [fname, "--algorithm", self.algorithm ] )
        df = pandas.read_table( "TResults/_opt{}result.txt".format( fileSig ), sep='\s+', skiprows=3, nrows=numHossParams, engine="python", header = None )
        assert(os.path.exists( "TExpts/fs_{}_20_pk.json".format(fileSig) ) )
        for rr in self.rowList:
            iff = int( round( rr['freq'] ) )
            #print( "remove TExpts/fs_{}_{}_pk.json".format( fileSig, iff ) )
            os.remove( "TExpts/fs_{}_{}_pk.json".format( fileSig, iff ) )
            os.remove( "TExpts/fs_{}_{}_pp.json".format( fileSig, iff ) )
        paramNames = ["cell", "exc", "numSq", "pattern", "idx", "initScore", "finalScore"] + list(df.iloc[:,0])
        initScore = 0.0
        finalScore = 0.0
        with open( "TResults/_opt{}result.txt".format(fileSig) , "r") as fp:
            last_line = fp.readlines()[-1]
            sp = last_line.replace( ",", " ").split( " " )
            initScore = float(sp[3])
            finalScore = float( sp[-1])
        return paramNames, [self.cell, self.exc, self.numSq, self.pattern, self.idx, initScore, finalScore] + list( df.iloc[:, 2] )
        '''
        df = pandas.read_table( "TResults/_opt{}results.txt".format( fileSig ), sep='\s+', skiprows=3, nrows=8, engine="python", header = None )
        paramNames = ["cell", "exc", "numSq", "pattern"] + list(df.iloc[:,0])
        return paramNames, [self.cell, self.exc, self.numSq, pattern] + list( df.iloc[:, 2] )
        '''

    def ticker( self, arg ):
        # arg should be the return value from the fitcell func
        #print( "{}_{}_{}    ".format( self.cell, self.exc, self.numSq ), end = "" )
        print( ".", end = '' )
        sys.stdout.flush()



    def dumpFindSim( self, means, freq, isPP ):
        fileSig = "{}_{}_{}_{}_{}".format( self.cell, self.exc, self.numSq, self.pattern, self.idx )
        stimAmpl = 50.0
        Ca_basal = 0.08
        stimWidth = 0.002
        pkDelay = self.pkDelay * 0.001
        windowStartt = -0.010
        windowEndt = 0.005
        pd = [pkDelay] * 9
        label = "Exc" if self.exc else "Inh"
        pkLabel = "pp" if isPP else "pk"

        transmitter = "glu" if self.exc else "GABA"
        fsStartStr  = '''{{
    "FileType": "FindSim",
    "Version": "1.0",
    "Metadata": {{
        "transcriber": "Upi Bhalla",
        "organization": "NCBS",
        "email": "bhalla@ncbs.res.in",
        "source": {{
            "sourceType": "other",
            "PMID": 0,
            "year": 2022,
            "figure": ""
        }},
        "testModel": "./Models/{}Presyn.g",
        "testMap": "./Maps/mapPresyn.json"
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
            "entity": "CaInput",
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
        "entities": [
            "{0}"
        ],
        "field": "conc",
        "window": {{"startt": {1}, "endt": {2}, "dt": 0.001,
            "operation": "{3}"
        }},
        "normalization": {{
            "entities": [
                "{5}"
            ],
            "sampling": "{4}"
        }},
        "data": ['''.format( "L_R", windowStartt, windowEndt, "oscPk", "start", "L_R")

        fsEndStr = '\n        ]\n  }\n}\n'

        fname = "TExpts/fs_{}_{}_{}.json".format( fileSig, freq, pkLabel )
        settleTime = 1.0
        with open( fname, "w" ) as fp:
            ############# Write the Ca stim sequence ################
            fp.write( fsStartStr.format( label ) )
            for idx in range( len( means ) ):
                if isPP and idx > 2 :   # Only do first 3 points: ref, p0,p1
                    continue
                t = settleTime if idx == 0 else settleTime + 0.3 + (idx-1)/freq
                fp.write( ",\n                [{:.4f}, {:.3f}],\n".format( t + pd[idx], stimAmpl ) )
                fp.write( "                [{:.4f}, {:.3f}]".format( stimWidth + t + pd[idx], Ca_basal ) )
            fp.write( fsMidStr )

            ############# Now on to the expt/sim output ################
            comma = ""
            #for idx, dd in enumerate( deconv ):
            #print( "LENS = ",  len( meanSTP ), len( meanTpk ) )
            #print( meanSTP )
            refPk = means[0][3]
            for idx, [tval, val, tpk, pk] in enumerate( means ):
                if isPP and idx > 2:   # Only do first 2 points.
                    continue
                ## Nasty. Skip the reference valley, and instead use the
                ## reference peak as the first point so it can be the
                ## normalization reference in the simulations.
                if idx > 0:
                    fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".format( 
                        comma, tval + 0.8, 
                        fracVesicleRelease * val/refPk ) )
                fp.write( "{}\n            [{:.4f}, {:.4f}, 0]".format( 
                    comma, tpk + 0.8, 
                    fracVesicleRelease * pk/refPk ) )
                comma = ","
            fp.write( fsEndStr )

    def dumpHoss( self ):
        excLabel = "Exc" if self.exc else "Inh"
        fileSig = "{}_{}_{}_{}_{}".format( self.cell, self.exc, self.numSq, self.pattern, self.idx )
        fname = "TConfigs/run{}.json".format( fileSig )
        exptStr = ''
        comma = "\n                    "
        for rr in self.rowList:
            iff = int( round( rr['freq'] ) )
            exptStr += '{}"fs_{}_{}_pk.json": {{"weight": {} }}'.format(
                comma, fileSig, iff, 200 if iff in [20, 50] else 100 )
            comma = ",\n                    "
            if self.doPairedPulse:
                exptStr += '{}"fs_{}_{}_pp.json": {{"weight": {} }}'.format(
                    comma, fileSig, iff, 200 if iff in [20, 50] else 100 )
    
        startStr = '''{{
    "FileType":"Hoss",
    "Version":"2.0",
    "author": "TabPresyn program",
    "model": "Models/{0}Presyn{1}.g",
    "map":"Maps/mapPresyn.json",
    "exptDir": "./TExpts",
    "outputDir": "./TResults",
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
    ]\n}}'''.format( excLabel, self.presynModelNumber, exptStr, hossParams, fileSig)

        with open( fname, "w" ) as fp:
            fp.write( startStr )

def main():
    global fracVesicleRelease
    parser = argparse.ArgumentParser( description = "Fits presynaptic release model to STP curves from pandas file. Required directories: Configs, Models, Expts, Results, Maps" )
    parser.add_argument( "-o", "--output", type = str, help = "Optional: output pandas hdf file for model params, default = model_params_trials.h5", default = "model_params_trials.h5")
    parser.add_argument( "-f", "--file", type = str, 
            help = "Optional: Name of tabulated STP file in pandas hdf5 format. Default = 'STP_trial_pks_and_refs.h5'", 
            default = "STP_trial_pks_and_refs.h5" )
    parser.add_argument( "-pd", "--peakDelay", type = float, 
            help = "Optional: Delay for time of glu peak in ms. Default = 10", default = 10.0)
    parser.add_argument( '-n', '--numProcesses', type = int, default = 0, help='Optional: Number of processes to spawn' )
    parser.add_argument( '-pmn', '--presynModelNumber', type = int, default = 75, help='Optional: Presynaptic model number' )
    parser.add_argument( '-fvr', '--fracVesicleRelease', type = float, default = 1.0, help='Optional: Fraction of total vesicle count released on first pulse. default 1.0' )
    parser.add_argument( '-pp', '--pairedPulse', action = "store_true", help='Flag: Turns on dumping of experiments to do paired pulse matching' )
    parser.add_argument( '-a', '--algorithm', type = str, help='Optional: Algorithm name to use, from the set available to scipy.optimize.minimize. Options are CG, Nelder-Mead, Powell, BFGS, COBYLA, SLSQP, TNC. The library has other algorithms but they either require Jacobians or they fail outright. There is also L-BFGS-B which handles bounded solutions, but this is not needed here because we already take care of bounds. COBYLA works well and is the default.', default = "COBYLA" )

    args = parser.parse_args()
    dat = pandas.read_hdf( args.file )
    cellList = dat['cell'].unique()
    print( "Cell list = ", cellList, len( cellList ) )

    numProcesses = args.numProcesses

    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count() // 4
    if numProcesses == 0:
        numProcesses = 1
    numProcesses = min( numProcesses, len( cellList ) * 4 )
    fracVesicleRelease = args.fracVesicleRelease

    print( "NUM PROCESSES = ", numProcesses )
    pool = multiprocessing.Pool( processes = numProcesses )
    params = []
    ret = []
    for cell in cellList:
        for exc in [ 0, 1 ]:
            for numSq in [ 5, 15 ]:
                celldf = dat.loc[dat['cell']==cell]
                freqList = celldf['freq'].unique()
                patterns = celldf.loc[celldf['numSq']==numSq]['pattern'].unique()
                for pp in patterns:
                    pattern = int( round( pp ) )
                    #print( "cell={} exc={} numSq={} pattern={}".format( int(cell), exc, numSq, pattern ))
                    rows = {}
                    numTrials = []
                    for ff in freqList:
                        iff = int( round( ff ) )
                        rows[iff] = []
                        patdf = celldf.loc[
                            (celldf['freq']==iff) & 
                            (celldf['exc']==exc) & 
                            (celldf['numSq']==numSq) &
                            (celldf['pattern']==pattern)
                        ]
                        numTrials.append( len( patdf ) )
                        for idx in range( len( patdf ) ):
                            rows[iff].append( patdf.iloc[idx,:] )
                    # Check that all the cases are repeated same # of times.
                    if numTrials.count(numTrials[0]) != len(numTrials):
                        print( "Warning: Inconsistent numTrials = ", numTrials,
                                "\nIn freqs: ", rows.keys() )
                        print( "Discarding smallest entry" )
                        mt = min( numTrials )
                        if numTrials.count( mt ) == 1:
                            for key, value in rows.items():
                                if len( value ) == mt:
                                    rows.pop( key )
                                    print( "popping entry for ", key )
                                    break
                        else:
                            quit()

                    for idx in range( numTrials[0] ):
                        r2 = []
                        for iff in rows:
                            #iff = int( round( ff ) )
                            r2.append( rows[iff][idx] )
                            #print( "Appended ", rows[iff][idx]["cell"], iff, idx )
                        #print( "DONE APPPPPPPPPP", flush = True )
                        ev = EvalFunc( cell, r2, exc, numSq, pattern, idx, args )
                        ret.append( pool.apply_async( ev.fitCell, callback=ev.ticker ))
    ans = [ rr.get() for rr in ret ]

    for aa in ans:
        [paramNames, paramVals] = aa
        params.append( paramVals )
    df = pandas.DataFrame( params, columns = paramNames )
    df.info()
    df.to_hdf( args.output, "w" )

    #pool.join()
    pool.close()

if __name__ == "__main__":
    main()
