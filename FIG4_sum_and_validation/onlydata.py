import moose
import pandas
import pylab
import numpy as np
import math
import matplotlib.pyplot as plt
import rdesigneur as rd
import rdesigneur.rmoogli as rmoogli
import argparse

freq = 80.0 # Hz
settleTime = 0.1    # seconds
stimDuration = 0.002   # seconds
numPulses = 16
#stimEnd = settleTime + numPulses/(freq+1)
postStim = 0.4
stimAmpl = 5e-2     # mM
basalCa = 0.08e-3   # mM
GABAdelay = 5.0e-3  # seconds
width = 0.002
firstPlotEntry = ['soma', '1', '.', 'Vm', 'Membrane potential']
#runtime = settleTime + numPulses / freq + postStim
numSq = 15

gluStimStr = "8e-5"
GABAStimStr = "8e-5"
gluR_clamp_potl = "-0.07"
GABAR_clamp_potl = "0.0"
GABAR_clamp_offset = 0.1    # nA
gluConductanceScale = 0.5   # Relative to default value in the spine proto
gluTau2Scale = 4   # Relative to default value in the spine proto

numCA1Exc = 100
numCA1Inh = 200
pCA3_CA1 = 0.0002
pCA3_Inter = 0.0008
pInter_CA1 = 1.0/256.0
interState = 0
thresh_CA3_Inter = 0.9999   # Avoid doing exact float comparisons to cross thresh.
thresh_CA3_CA1 = 0.9999
thresh_Inter_CA1 = 0.9999
repeatPatterns = False
inputs = []
stimList = []
pulseTrig = []
pulseThresh = 0.001
patternData = "../../../2022/VC_DATA/all_cells_SpikeTrain_CC_long.h5"
SAMPLE_FREQ = 20000
chemDt = 0.0005
SAMPLE_TIME = 11
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
PKDELAY = int( 0.020 * SAMPLE_FREQ )
PKTHRESH = 0.0005 # Threshold for a distinct EPSP pk. In Volts
ALPHATAU = 0.005 * SAMPLE_FREQ
ALPHAWINDOW = int( ALPHATAU * 4 )
ALPHA = np.array( [ (t/ALPHATAU)*np.exp(1-t/ALPHATAU) for t in range( ALPHAWINDOW ) ] )

moogList = [
    ['#', '1', '.', 'Vm', 'Membrane potential', -70.0, 0.0],
    ['#', '1', 'glu/Ca', 'n', 'Glu Ca n', 0.0, 1000.0],
    ['#', '1', 'GABA/Ca', 'n', 'GABA Ca n', 0.0, 1000.0]
]

def patternDict():
    patternZeros = [0]*64
    patternOnes = [1]*64

    patternA =  [
             0,0,0,1,1,0,0,0,
             0,0,1,0,0,1,0,0,
             0,0,1,0,0,1,0,0,
             0,1,0,0,0,0,1,0,
             0,1,1,1,1,1,1,0,
             1,1,0,0,0,0,1,1,
             1,0,0,0,0,0,0,1,
             1,0,0,0,0,0,0,1,]

    patternB =  [
             0,1,1,1,1,1,0,0,
             0,1,0,0,0,1,0,0,
             0,1,0,0,0,1,0,0,
             0,1,1,1,1,0,0,0,
             0,1,0,0,0,1,0,0,
             0,1,0,0,0,0,1,0,
             0,1,1,1,1,1,0,0,
             0,0,0,0,0,0,0,0,]

    patternC =  [
             0,1,1,1,1,1,1,0,
             1,1,0,0,0,0,0,1,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,1,
             0,1,1,1,1,1,1,0,]

    patternD =  [
             0,0,0,0,0,0,0,0,
             1,1,1,1,1,1,0,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,1,1,1,1,1,0,0,]

    patternE =  [
             0,1,1,1,1,1,1,0,
             0,1,1,0,0,0,0,0,
             0,1,1,0,0,0,0,0,
             0,1,1,1,0,0,0,0,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,0,0,0,
             0,1,1,1,1,1,1,0,]

    patternF =  [ 
             0,0,0,0,0,0,0,0,
             1,1,1,1,1,1,1,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,1,1,1,1,1,0,0,
             1,1,1,1,1,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,]

    patternG =  [
             0,0,1,1,1,1,1,0,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,1,1,0,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,0,1,1,1,1,1,0,]

    patternH =  [
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,1,1,1,1,1,1,0,
             1,0,0,0,0,1,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,
             1,0,0,0,0,0,1,0,]

    patternI =  [
             0,0,1,1,1,1,0,0,
             0,0,0,1,1,0,0,0,
             0,0,0,1,1,0,0,0,
             0,0,0,1,1,0,0,0,
             0,0,0,1,1,0,0,0,
             0,0,0,1,1,0,0,0,
             0,0,0,1,1,0,0,0,
             0,1,1,1,1,1,1,0,]

    return {
            '0': np.array( patternZeros, dtype = float),
            '1': np.array( patternOnes, dtype = float),
            'A': np.array( patternA, dtype = float),
            'B': np.array( patternB, dtype = float),
            'C': np.array( patternC, dtype = float),
            'D': np.array( patternD, dtype = float),
            'E': np.array( patternE, dtype = float),
            'F': np.array( patternF, dtype = float),
            'G': np.array( patternG, dtype = float),
            'H': np.array( patternH, dtype = float),
            'I': np.array( patternI, dtype = float)
            }

def generatePatterns( seq, seed ):
    global inputs
    global CA3_Inter
    global CA3_CA1
    global Inter_CA1
    pd = patternDict()
    inputs = []
    for cc in seq:
        inputs.append( pd[cc] )

    # Assumes all are random projections.
    np.random.seed( seed )
    CA3_Inter = (np.random.rand( 256, 256 ) < pCA3_Inter) * 1.0
    CA3_CA1 = (np.random.rand( numCA1Exc, 256 ) < pCA3_CA1) * 1.0
    Inter_CA1 = (np.random.rand( numCA1Inh, 256 ) < pInter_CA1) * 1.0



class MooArg:
    def __init__( self, title, field ):
        self.diaScale = 1.0
        self.ymin = 0
        self.ymax = 1000
        self.title = title
        if self.title == "Inter":
            self.ymin = -1000
        self.field = field
        self.relpath = "."


def makeInputs( name, xOffset ):
    spacing = 2e-6
    size = spacing / 1.5
    yOffset = 15e-6
    zOffset = 0.0
    CA3 = moose.Neutral( "/model/elec/" + name )
    CA3cells = moose.Compartment( "/model/elec/{}/soma".format(name), 256 )
    for idx, compt in enumerate( CA3cells.vec ):
        ii = idx % 16
        jj = 16 - idx // 16
        compt.x0 = compt.x = ii * spacing + xOffset
        compt.y0 = compt.y = jj * spacing + yOffset
        compt.z = zOffset
        compt.z0 = zOffset + size
        compt.diameter = size
        compt.tick = -1
    #moose.showfield( compt )
    #print( "LEN = ", len( moose.wildcardFind( CA3.path + "/#[ISA=CompartmentBase]" ) ) )


    ma = MooArg( name, "Vm" )

    return CA3cells.vec, ma 


def buildModel( presynModelName, seed, useGssa, vGlu, vGABA, spiking ):
    rGlu = pow( vGlu, 1.0/3.0)
    rGABA = pow( vGABA, 1.0/3.0)
    NaGbar = 400.0 if spiking else 6.0
    KGbar = 450.0 if spiking else 3.5

    rdes = rd.rdesigneur(
        elecDt = 50e-6,
        chemDt = chemDt,
        funcDt = 0.0005,
        diffDt = 0.0005,
        chemPlotDt = 0.0005,
        elecPlotDt = 0.0005,
        diffusionLength = 1e-3, # ensure single chem compt
        useGssa = useGssa,
        # cellProto syntax: ['ballAndStick', 'name', somaDia, somaLength, dendDia, dendLength, numDendSeg]
        cellProto = [['ballAndStick', 'soma',  10e-6, 10e-6, 2e-6, 200e-6, 1]],
        spineProto = [['makeActiveSpine()', 'spine']],
        spineDistrib = [
            # Put 100 spines, 2 um apart
            ['spine', 'dend#', '2e-6', '-1e-6']  
            ],
        chemProto = [
            ['Models/{}'.format( presynModelName ), 'chem'],
        ],
        chanProto = [
            ['make_Na()', 'Na'],
            ['make_K_DR()', 'K'],
            ['make_GABA()', 'GABA']
        ],
        # Some changes to the default passive properties of the cell.
        #passiveDistrib = [['soma', 'CM', '0.03', 'Em', '-0.065']],
        passiveDistrib = [['soma', 'CM', '0.01', 'Em', '-0.065']],
        chemDistrib = [
            # Args: chem_model, elec_compts, mesh_type, spatial_distrib, r_scalefactor, radius_sdev
            ['glu', 'head#', 'presyn_spine', '1', rGlu, 0 ],
            # Args: chem_model, elec_compts, mesh_type, spatial_distrib, r_absolute, r_sdev, spacing
            ['GABA', 'dend#', 'presyn_dend', '1', 0.5e-6 * rGABA, 0, 1e-6 ],
        ],
        chanDistrib = [
            ['Na', 'soma', 'Gbar', str(NaGbar) ],
            ['K', 'soma', 'Gbar', str(KGbar) ],
            ['GABA', 'dend#', 'Gbar', '3.2' ]
        ],
        adaptorList = [
            [ 'glu/glu', 'n', 'glu', 'activation', 0.0, 8000 ],
            [ 'GABA/GABA', 'n', 'GABA', 'activation', 0.0, 8000 ],
        ],
        stimList = stimList,
        plotList = [
            firstPlotEntry,
            ['head#', '1', 'glu', 'Ik', 'glu current'],
            ['dend#', '1', 'GABA', 'Ik', 'GABA current'],
            ['head#', '1', 'glu/Ca', 'n', 'glu_presyn_Ca'],
            ['dend#', '1', 'GABA/Ca', 'n', 'GABA_presyn_Ca'],
            ['head#', '1', 'glu/Docked', 'n', 'glu_Docked'],
            ['dend#', '1', 'GABA/Docked', 'n', 'GABA_Docked'],
            ['head#', '1', 'glu/RR_pool', 'n', 'glu_RR'],
            ['dend#', '1', 'GABA/RR_pool', 'n', 'GABA_RR'],
            ['head#', '1', 'glu/glu', 'n', 'glu_released'],
            ['dend#', '1', 'GABA/GABA', 'n', 'GABA_released'],
        ],
        moogList = moogList
    )
    moose.seed( seed ) 
    gluReceptor = moose.element( '/library/spine/head/glu' )
    gluReceptor.Gbar *= gluConductanceScale # Tweak conductance
    gluReceptor.tau2 *= gluTau2Scale # Tweak closing time
    moose.element( '/library/GABA' ).Ek = -0.07 # Tweak Erev.
    rdes.buildModel()
    if useGssa:
        #moose.showfield( "/model/chem/glu/ksolve" )
        moose.element( "/model/chem/glu/ksolve" ).useClockedUpdate = 1
        moose.element( "/model/chem/GABA/ksolve" ).useClockedUpdate = 1
    # Here we can add a few messages from the periodic synaptic input to the
    # Ca.
    
    #rs = moose.element( '/model/elec/soma/GABAR/sh/synapse/synInput_rs' )
    #moose.connect( rs, "spikeOut", "/model/chem/dend/glu/Ca/tab", "input" )
    #moose.reinit()
    return rdes


'''
def getAUC( plot, freq ):
    p1 = moose.element( plot )
    vec = p1.vector

    numPtsInPk = int( 1.0 / ( p1.dt * freq ) )
    print( "numPtsInPk = ", numPtsInPk )
    startIdx = int ( 1.0 / p1.dt )
    auc = []
    for i in range( numPulses ):
        nextIdx = startIdx + numPtsInPk
        auc.append( sum( vec[startIdx:nextIdx] ) )
        startIdx = nextIdx
    return -1e9 * np.array( auc )
'''

def addSubplot( subplotIdx, rdesIdx, field, freq, sequence ):
    plt.subplot( 6, 1, subplotIdx )
    p1 = "/model/graphs/plot{}".format( rdesIdx )
    p2 = "/model/graphs/plot{}".format( rdesIdx + 1 )
    e1 = moose.element( p1 )
    e2 = moose.element( p2 )
    dt = e1.dt
    d1 = e1.vector
    d2 = e2.vector
    ylabel = "# of molecules" if subplotIdx != 1 else "Current (pA)"
    scaling = 1.0e12 if subplotIdx == 1 else 1.0
    t = np.arange( 0, len( d1 ) * dt - 1e-6, dt )
    plt.plot( t, d1 * scaling, "r", label = "glu."+field )
    plt.plot( t, d2 * scaling, "b", label = "GABA."+field )
    plt.xlabel("Time")
    plt.ylabel( ylabel )
    plt.legend()
    if subplotIdx == 1:
        plt.title( "Freq = {} Hz, Seq = {}".format( freq, sequence ) )

def stimFunc():
    global interState
    t = moose.element( '/clock' ).currentTime
    stimEnd = SAMPLE_TIME + settleTime
    if t < settleTime or t > stimEnd:
        return
    ts = t - settleTime
    patternNum = 0
    #patternNum = math.floor( ts * freq )
    #print( "Pattern Num = {}, t = {}".format( patternNum,  t ) )
    #CA3isActive = patternNum > math.floor( (ts-stimDuration)*freq )
    #InterIsActive = math.floor( (ts-GABAdelay)*freq ) > math.floor( (ts-GABAdelay-stimDuration)*freq )
    pulseIdx = int( round( ts / chemDt ) )
    if pulseIdx >= len( pulseTrig ) :
        #print( "ERROR: Pulse idx = {} > {} at ts = {}".format( pulseIdx, len( pulseTrig ), ts ) )
        CA3isActive = False
        InterIsActive = False
    else:
        CA3isActive = (pulseTrig[pulseIdx] > pulseThresh )
        pulseIdx = int( round( (ts - GABAdelay) / chemDt ) )
        InterIsActive = (pulseTrig[pulseIdx] > pulseThresh)
    if ( pulseIdx % 500 ) == 0:
        print( "PL = ", len( pulseTrig), " PI = ", pulseIdx, " CA3:", CA3isActive, InterIsActive )

    gluInput = moose.vec( "/model/chem/glu/Ca_ext" )
    gabaInput = moose.vec( "/model/chem/GABA/Ca_ext" )
    thresh = 2.0
    if CA3isActive:
        ca3cells = moose.vec( "/model/elec/CA3/soma" )
        pat = inputs[patternNum % len( inputs ) ]
        pat = pat.reshape( 8, 8 )
        if numSq == 5:
            ca3cells.Vm = np.append(pat.repeat( 2, axis=0 ).reshape(-1), np.zeros( 128 ) )
        else:
            ca3cells.Vm = pat.repeat( 4, axis=0 ).reshape(-1)
        #ca3cells.Vm = np.repeat( inputs[ patternNum % len( inputs ) ], 4 )
        Inter = moose.vec( "/model/elec/Inter/soma" )
        Inter.Vm = (np.matmul( CA3_Inter, ca3cells.Vm ) >= thresh_CA3_Inter ) * 1.0
        gluInput.concInit = (np.matmul( CA3_CA1, ca3cells.Vm ) >= thresh_CA3_CA1 ) * stimAmpl
        '''
        print( "MEAN CA3_Inter = ", np.mean( np.matmul( CA3_Inter, ca3cells.Vm ) ),
            " CA3_CA1 = ", np.mean( np.matmul( CA3_CA1, ca3cells.Vm ) ),
            " Inter_CA1 = ", np.mean( np.matmul( Inter_CA1, Inter.Vm ) ),
            )
        '''
    else:
        moose.vec( "/model/elec/CA3/soma" ).Vm = 0.0
        #moose.vec( "/model/elec/Inter/soma" ).Vm = 0.0
        gluInput.concInit = basalCa

    if InterIsActive:
        Inter = moose.vec( "/model/elec/Inter/soma" )
        gabaInput.concInit = (np.matmul( Inter_CA1, Inter.Vm ) >= thresh_Inter_CA1 ) * stimAmpl
        interState = 1
    else:
        gabaInput.concInit = basalCa
        if interState == 1:
            interState = 0
            moose.vec( "/model/elec/Inter/soma" ).Vm = 0
    #print( "{:.4f}  CA3={}, Inter={}, pattern={}".format( t, CA3isActive, InterIsActive, patternNum))

def makeNetwork( rdes ):
    origNeuronId = rdes.elecid
    CA3cells, CA3args = makeInputs( "CA3", 20e-6 )
    interneurons, interneuronArgs = makeInputs( "Inter", 70e-6 )
    if len( moogList ) > 0:
        rdes.elecid = moose.element( "/model/elec/CA3" )
        rdes.moogNames.append( rmoogli.makeMoogli( rdes, CA3cells, CA3args, rd.knownFieldsDefault['Vm'] ) )
        rdes.elecid = moose.element( "/model/elec/Inter" )
        rdes.moogNames.append( rmoogli.makeMoogli( rdes, interneurons, interneuronArgs, rd.knownFieldsDefault['Vm'] ) )
        rdes.elecId = origNeuronId

def findPeaks( pulseTrig, pulseThresh, Vm, width = 0.005 ):
    widthSamples = int( np.round( width * SAMPLE_FREQ ) )
    half = widthSamples // 2
    trigIdx = []
    lastPP = 0.0
    for idx, pp in enumerate( pulseTrig ):
        if pp > pulseThresh and lastPP < pulseThresh:
            trigIdx.append( idx )
            #print( idx )
        lastPP = pp

    pks = []
    pkIdx = []
    for idx in trigIdx:
        idx2 = idx + PKDELAY
        # NOTE that maxIdx is relative to idx2, ie, within a window.
        maxIdx = np.argmax( Vm[idx2-widthSamples:idx2+widthSamples] )
        #print( "{:.3f}".format( idx2 / SAMPLE_FREQ ) )
        if maxIdx > half and maxIdx < (3*half):
            thisMin = min(Vm[ idx2-widthSamples-half : idx2-widthSamples ])
            thisMax = Vm[idx2 + maxIdx - widthSamples]
            #print( "thisdelta", thisMax - thisMin, PKTHRESH, len( trigIdx), len( pulseTrig ))
            if (thisMax - thisMin) > PKTHRESH:
                print( "thisPK t = {:.3f}, pk = {:.3f}".format( idx, thisMax ) )
                pks.append(thisMax - thisMin)
                pkIdx.append( idx )

    print( "NUMPKS = ", len( pkIdx ), len( trigIdx ) )
    return pkIdx, pks




def main():
    global freq
    global moogList
    global firstPlotEntry
    global stimList
    global numPulses
    global numSq
    global pInter_CA1
    global pulseTrig
    global pulseThresh

    parser = argparse.ArgumentParser( description = "Deliver patterned stims to neuron with short-term-plasticity in E and I synapses" )
    parser.add_argument( "-f", "--freq", type = float, help = "Specify pattern stim freq, default 20 Hz", default = 20 )
    parser.add_argument( "-nf", "--numFrames", type = int, help = "Number of frames to generate for input, default = 16", default = 16 )
    parser.add_argument( "-nsq", "--numSq", type = int, help = "Number of squares in input pattern. Default = 5", default = 5 )
    parser.add_argument( '-seq', '--sequence', type = str, help ='Sequence of patterns to display. Can be any of 01ABCDEFGHI.', default = "A" )
    parser.add_argument( '-g', '--graphics3D', action="store_true", help ='Flag: when set, display 3-d graphs using rmoogli.' )
    parser.add_argument( '-spk', '--spiking', action="store_true", help ='Flag: when set, use high Na/K channel densities in soma to get spiking.' )
    parser.add_argument( '-v', '--voltage_clamp', action="store_true", help ='Flag: when set, do voltage clamp for glu and GABA currents respectively.')
    parser.add_argument( '-d', '--deterministic', action="store_true", help ='Flag: when set, use deterministic ODE solver. Normally uses GSSA stochastic solver.')
    parser.add_argument( "-m", "--modelName", type = str, help = "Optional: specify name of presynaptic model file, assumed to be in ./Models dir.", default = "BothPresyn75.g" )
    parser.add_argument( "-s", "--seed", type = int, help = "Optional: Seed to use for random numbers both for Python and for MOOSE.", default = 1234 )
    parser.add_argument( "-vglu", "--volGlu", type = float, help = "Optional: Volume scaling factor for Glu synapses. Default=1", default = 1.0 )
    parser.add_argument( "-vGABA", "--volGABA", type = float, help = "Optional: Volume scaling factor for GABA synapses. Default=0.5", default = 0.5 )
    parser.add_argument( "-ic", "--inhibitoryConvergence", type = float, help = "Optional: Convergence of projections from Inh interneurons to inhibitory synapses on CA1. Default=1.0 ", default = 1.0 )
    args = parser.parse_args()

    freq = args.freq
    if not args.graphics3D:
        moogList = []

    numPulses = args.numFrames
    if numPulses != len( args.sequence ):
        print( "Warning: numPulses ({}) != sequence length {}. Will cycle through.".format( numPulses, len( args.sequence ) ) )
    numSq = args.numSq
    pInter_CA1 = args.inhibitoryConvergence / 256.0

    if args.voltage_clamp:
        moogList = []   # Shut off 3d display for this dual run.
        stimList = [['soma', '1', '.', 'vclamp', '-0.070' ]]
        firstPlotEntry = ['soma', '1', 'vclamp', 'current','Vclamp current']
        # Possibly here we need to increase the vclamp gain.
    else:
        firstPlotEntry = ['soma', '1', '.', 'Vm', 'Membrane potential']

    generatePatterns( args.sequence, args.seed )

    rdes = buildModel( args.modelName, args.seed, not args.deterministic, args.volGlu, args.volGABA, args.spiking )
    pr = moose.PyRun( "/model/stims/stimRun" )
    #pr.initString = 'print(freq)'
    pr.runString = 'stimFunc()'
    pr.tick = 14

    makeNetwork( rdes )

    # Set up the stimulus timings
    df = pandas.read_hdf( patternData )
    # NOTE: Iterate over cells here before going to individual pulses.
    pulseTrig = np.array(df.iloc[SWEEP, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    pulseThresh = ( min( pulseTrig ) + max( pulseTrig ) ) / 2.0
    print( "LEN PULSE TR=", len( pulseTrig ) )
    epsc = np.array(df.iloc[SWEEP, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ])
    tepsc = np.linspace( 0, SAMPLE_TIME, len(epsc) )
    ipsc = np.array(df.iloc[SWEEP + 24, SAMPLE_START:SAMPLE_START+NUM_SAMPLES ] )

    pkIdx, pks = findPeaks( pulseTrig, pulseThresh, epsc )
    fitEPSP = np.zeros( len( epsc ) )
    for idx, pp in zip( pkIdx, pks ):
        fitEPSP[idx:idx+ALPHAWINDOW] += ALPHA*pp

    #moose.reinit()
    '''
    print( "NumGluR = ", len( moose.vec('/model/chem/glu/glu') ),
        "  NumGABA sites = ", 
        len(moose.vec( '/model/chem/GABA/GABA' )) )
    print( "NumSpines = ", len(moose.wildcardFind( '/model/elec/head#' )),
            "  NumGABA sites = ", 
            len(moose.wildcardFind( '/model/chem/GABA/GABA[]' )) )
    '''
    #runtime = settleTime + numPulses / freq + postStim
    runtime = settleTime + SAMPLE_TIME + postStim
    if len( moogList ) > 0:
        rdes.displayMoogli( 0.001, runtime, rotation = 0.00, mergeDisplays=True, colormap = "plasma" )


    plt.rcParams.update( {"font.size": 24} )
    for freq in [ args.freq ]:
        #runtime = settleTime + numPulses / freq + postStim
        if args.voltage_clamp:
            moose.element( "/model/stims/stim0" ).expr = gluR_clamp_potl
            #vc = moose.element( '/model/elec/soma/vclamp' )
            #vc.gain *= 1.0
        #print( freq )
        #moose.reinit()
        #moose.start( runtime )
        #gluAuc = getAUC( '/model/graphs/plot1', freq )
        #GABAAuc = getAUC( '/model/graphs/plot2', freq )
        fig = plt.figure( "StimFreq = " + str( freq ), figsize = (11.8, 3.5) )
        #addSubplot( 1, 1, "Current", freq, args.sequence  )
        #addSubplot( 2, 3, "Ca", freq, args.sequence  )
        #addSubplot( 3, 5, "Docked", freq, args.sequence  )
        #addSubplot( 4, 7, "RR_pool", freq, args.sequence  )
        #addSubplot( 5, 9, "released", freq, args.sequence  )

        #ax1 = plt.subplot( 1, 1, 1 )
        #ax1.title.set_fontsize(12)
        #plt.title( "Freq = {} Hz, NumSq = {}, Seq = {}".format( freq, args.numSq, args.sequence ), size = 12 )
        offset = -90.0 if args.spiking else -70.0
        #plt.text(1, offset, args.sequence, fontsize = 24, fontfamily = 'monospace')
        #plot0 = moose.element( '/model/graphs/plot0' ).vector
        dt = moose.element( '/model/graphs/plot0' ).dt
        print( "DT = ", dt )
        tpt = np.arange( 0.0, len( pulseTrig ), 1 )/SAMPLE_FREQ
        #t = np.arange( 0, len( plot0 ) * dt - 1e-6, dt )
        #stimY = np.array( [(0.0 if (tt < settleTime or tt > stimEnd) else float((int((tt-1.0)*freq) % 8))) for tt in t] )
        if args.voltage_clamp:
            plot0[:10] = 0  # clear out transient
            plt.plot( tpt, pulseTrig * 100, "g", label = "pulseTrig" )
            plt.plot( t, plot0 * 1e12, "r", label = "I_gluR" )
            #y = [ min( plot0 ) * 1000, max( plot0 ) * 1000 ]
            #plt.plot( t, stimY * 1e9, "g", label = "stimStage" )
            moose.element( "/model/stims/stim0" ).expr = GABAR_clamp_potl
            moose.reinit()
            #moose.start( runtime )
            #plot0 = moose.element( '/model/graphs/plot0' ).vector
            #plot0[:10] = GABAR_clamp_offset / 1e12  # clear out transient
            #plt.plot( t, plot0 * 1e12 - GABAR_clamp_offset, "b", label = "I_GABAR" )
            plt.plot( tepsc, epsc + 50, label = "expt_EPSC" )
            plt.plot( tepsc, ipsc + 150, label = "expt_IPSC" )
            plt.ylabel( "clamp curr (pA)" )
            plt.legend()

            #plt.ylabel( "ipsc (pA)" )
        else:
            #plt.plot( t, plot0 * 1000 + 65, "r", label = "Sim EPSP" )
            plt.plot( tpt, pulseTrig * 100, "g", label = "Trigger" )
            #y = [ offset, max( plot0 ) * 1000 ]
            plt.ylabel( "Vm (mV)" )
            plt.xlim( settleTime - 0.1, runtime + 0.1 )
            plt.plot( tepsc, epsc, "b", label = "Data EPSP" )
            plt.plot( tepsc, fitEPSP, "r", label = "Fit EPSP" )
            plt.ylabel( "EPSP (mV)" )
            plt.legend()
        # Plot the vertical lines separating 8 stimuli
        for ii in range( 0, args.numFrames+1, 8):
            x = settleTime + ii/freq
            #plt.plot( [x, x], y, "g:" )
        plt.xlim( settleTime - 0.1, runtime + 0.1 )
        plt.xlabel("Time (s)")
        fig.tight_layout()
    plt.show()

    
if __name__ == "__main__":
    main()
