import moose
import pandas
import matplotlib.pyplot as plt
import numpy as np
import math
import rdesigneur as rd
import multiprocessing
from pathlib import Path
import argparse

stimAmpl = 5e-2     # mM
basalCa = 0.08e-3   # mM
GABAdelay = 5.0e-3  # seconds
width = 0.002
plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Vm'],
        ['head37', '1', 'glu/Ca_ext', 'nInit', 'Stimulus'],
        ['head37', '1', 'glu', 'Ik', 'GluR current'],
        ['head37', '1', 'glu/Ca', 'n', 'Glu bouton Ca'],
        ['head37', '1', 'glu/RR_pool', 'n', 'RR pool'],
        ['head37', '1', 'glu/Docked', 'n', 'Docked'],
        ['head#', '1', 'glu/glu', 'n', 'Glu Released'],
        ['dend0', '1', 'GABA', 'Ik', 'GABAR current'],
        ['dend0', '1', 'GABA/Ca', 'n', 'Bouton Ca'],
        ['dend0', '1', 'GABA/RR_pool', 'n', 'RR pool'],
        ['dend0', '1', 'GABA/Docked', 'n', 'Docked'],
        ['dend#', '1', 'GABA/GABA', 'n', 'GABA released'],
]

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
pInter_CA1 = 0.004

#used as globals. 
CA3_CA1 = 0.0002
CA3_Inter = 0.0008
Inter_CA1 = 0.004

interState = 0
thresh_CA3_Inter = 0.9999   # Avoid doing exact float comparisons to cross thresh.
thresh_CA3_CA1 = 0.9999
thresh_Inter_CA1 = 0.9999
repeatPatterns = False
inputs = []
stimList = []
SAMPLE_FREQ = 20000
elecDt = 0.00005
chemDt = 0.0005
SAMPLE_TIME = 11
RUNTIME = 4
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
patternDict2 = {}

PulseTrain = np.array([4001,10684,11276,11603,13433,15914,16193,17131,19457,19827,20561,21153,21578,
    22460,24407,24665,25093,25667,26213,26726,27343,28046,28625,29322,29608,31223,31729,32400,32756,
    33317,33897,35890,36496,36986,37267,37484,38755,39890,40495,41873,42970,43399,45768,46100,46695,
    46931,47430,47639,47877,48568,49189,51579,52910,53373,53643,56169,56686,57112,57467,57834,58721,
    59254,60261,60473,61816,63607,64798,66090,66291,69446,70416,70666,70898,71145,71821,72805,73201,
    74279,74777,75520,76181,77447,77966,78309,79050,79331,80383,81575,82380,82991,85548,87622,88515,
    88839,89510,89866,90977,91257,91841,92837,93249,94872,95549,96164,96975,98498,99152,99545,99795,
    100493,101582,102149,103757,107075,107600,107969,108705,109143,109875,110347,110856,113988,114470,
    115634,116946,117489,118060,119694,121243,122078,122580,124326,125053,127211,128234,128814,129380,
    129945,130884,131133,131550,132432,133262,133560,134345,134707,135065,135938,136529,137450,137806,
    139055,140234,141304,143221,143573,144296,145640,145984,146846,147856,148671,150909,152493,152852,
    153268,153931,155048,155690,156475,157345,158850,159443,159768,160600,160919,161424,161660,161956,
    163448,163758,164107,165661,166052,166540,167119,168032,169773,170130,171780,172502,173106,174142,
    174728,175182,175694,176340,177236,178437,179524,180446,183258,183781,185319,187213,189396,190365,
    190837,191267,191619,192282,192848,193144,193689,194521,195822,196751,197884,199981,200689,201095,
    202108,203280,204018,205585,206552,207234,207796,209126,209832])

TRIG = np.zeros( NUM_SAMPLES )
for ii in PulseTrain:
    TRIG[ii] = 1.0

ReducedPulseIdx = np.zeros( round(SAMPLE_TIME / chemDt ), dtype=int )
for pp in PulseTrain:
    ReducedPulseIdx[round(pp / 10)] = 1


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


def generatePatterns( args ):
    global CA3_CA1
    global CA3_Inter
    global Inter_CA1
    global patternDict2
    pd = patternDict()

    # Assumes all are random projections.
    np.random.seed( args.seed )
    CA3_Inter = (np.random.rand( 256, 256 ) < args.pCA3_Inter) * 1.0
    CA3_CA1 = (np.random.rand( numCA1Exc, 256 ) < args.pCA3_CA1) * 1.0
    Inter_CA1 = (np.random.rand( numCA1Inh, 256 ) < args.pInter_CA1) * 1.0
    px = []
    for char in ["A", "B", "C", "D", "E"]:
        temp = pd[char].reshape(8,8).repeat(2, axis = 0)
        temp2 = np.array(temp)
        temp2[:,5:] = 0
        px.append( np.append(temp2.reshape(-1), np.zeros(128) ) )

    patternDict2 = {
        46:px[0],
        47:px[1],
        48:px[2],
        49:px[3],
        50:px[4],
        52:pd["F"].reshape(8,8).repeat( 4, axis=0 ).reshape(-1),
        53:pd["G"].reshape(8,8).repeat( 4, axis=0 ).reshape(-1),
        55:pd["H"].reshape(8,8).repeat( 4, axis=0 ).reshape(-1)
    }



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
        elecDt = elecDt,
        chemDt = chemDt,
        funcDt = chemDt,
        diffDt = chemDt,
        chemPlotDt = chemDt,
        elecPlotDt = elecDt,
        diffusionLength = 1e-3, # ensure single chem compt
        useGssa = useGssa,
        # cellProto syntax: ['ballAndStick', 'name', somaDia, somaLength, dendDia, dendLength, numDendSeg]
        cellProto = [['ballAndStick', 'soma',  10e-6, 10e-6, 2e-6, 200e-6, 1]],
        spineProto = [['makeActiveSpine()', 'spine']],
        spineDistrib = [
            # Put 100 spines, 2 um apart, with size scale 1, size sdev 0.5
            ['spine', 'dend#', '2e-6', '-1e-6', '1', '0.5']  
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
            ['glu', 'head#', 'presyn_spine', '1', rGlu, 0.5 ],
            # Args: chem_model, elec_compts, mesh_type, spatial_distrib, r_absolute, r_sdev, spacing
            ['GABA', 'dend#', 'presyn_dend', '1', 0.5e-6 * rGABA, 0, 1e-6 ],
        ],
        chanDistrib = [
            ['Na', 'soma', 'Gbar', str(NaGbar) ],
            ['K', 'soma', 'Gbar', str(KGbar) ],
            ['GABA', 'dend#', 'Gbar', '3.2' ]
        ],
        adaptorList = [
            [ 'glu/glu', 'n', 'glu', 'activation', 0.0, 4000 ],
            [ 'GABA/GABA', 'n', 'GABA', 'activation', 0.0, 4000 ],
        ],
        stimList = stimList,
        plotList = plotList,
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
    moose.reinit()
    return rdes


def stimFunc( patternIdx ):
    t = moose.element( '/clock' ).currentTime
    # Need to look up if this is time to generate pulse. 
    idx = int(round( t/chemDt ) )
    if idx % int( 1.0/chemDt )  == 0:   # dot every second.
        print( ".", flush = True, end = "" )
    if idx >= len( ReducedPulseIdx ):
        return
    CA3isActive = (ReducedPulseIdx[idx] > 0.5) #Stimulus is to be delivered
    idx2 = int( round( (t - GABAdelay) / chemDt ) )
    if idx2 >= len( ReducedPulseIdx ):
        return
    InterIsActive = ( ReducedPulseIdx[idx2] > 0.5 )
    gluInput = moose.vec( "/model/chem/glu/Ca_ext" )
    gabaInput = moose.vec( "/model/chem/GABA/Ca_ext" )
    thresh = 2.0
    if CA3isActive:
        #print( "Trig CA3 at {:.3f} {} with {}".format( t, idx, patternIdx ))
        ca3cells = moose.vec( "/model/elec/CA3/soma" )
        ca3cells.Vm = patternDict2[patternIdx]
        gluInput.concInit = (np.matmul( CA3_CA1, ca3cells.Vm ) >= thresh_CA3_CA1 ) * stimAmpl
        '''
        print( "MEAN CA3_Inter = ", np.mean( np.matmul( CA3_Inter, ca3cells.Vm ) ),
            " CA3_CA1 = ", np.mean( np.matmul( CA3_CA1, ca3cells.Vm ) ),
            " Inter_CA1 = ", np.mean( np.matmul( Inter_CA1, Inter.Vm ) ),
            )
        '''
    else:
        moose.vec( "/model/elec/CA3/soma" ).Vm = 0.0
        gluInput.concInit = basalCa

    if InterIsActive:
        Inter = moose.vec( "/model/elec/Inter/soma" )
        Inter.Vm = (np.matmul( CA3_Inter, patternDict2[patternIdx]) >= thresh_CA3_Inter ) * 1.0
        gabaInput.concInit = (np.matmul( Inter_CA1, Inter.Vm ) >= thresh_Inter_CA1 ) * stimAmpl
    else:
        gabaInput.concInit = basalCa
        moose.vec( "/model/elec/Inter/soma" ).Vm = 0
    #print( "{:.4f}  CA3={}, Inter={}, pattern={}".format( t, CA3isActive, InterIsActive, patternIdx))

def makeNetwork( rdes ):
    origNeuronId = rdes.elecid
    CA3cells, CA3args = makeInputs( "CA3", 20e-6 )
    interneurons, interneuronArgs = makeInputs( "Inter", 70e-6 )

def innerMain( args ):
    patternIdx = args.pattern
    if args.voltage_clamp:
        stimList = [['soma', '1', '.', 'vclamp', '-0.070' ]]
        firstPlotEntry = ['soma', '1', 'vclamp', 'current','Vclamp current']
        # Possibly here we need to increase the vclamp gain.
    else:
        firstPlotEntry = ['soma', '1', '.', 'Vm', 'Membrane potential']

    generatePatterns( args )

    rdes = buildModel( args.modelName, args.seed, not args.deterministic, args.volGlu, args.volGABA, args.spiking )
    pr = moose.PyRun( "/model/stims/stimRun" )
    pr.runString = 'stimFunc({})'.format( patternIdx )
    pr.tick = 14 # This would be chemDt. Which is currently 0.5 ms.

    makeNetwork( rdes )

    moose.reinit()
    '''
    print( "NumGluR = ", len( moose.vec('/model/chem/glu/glu') ),
        "  NumGABA sites = ", 
        len(moose.vec( '/model/chem/GABA/GABA' )) )
    print( "NumSpines = ", len(moose.wildcardFind( '/model/elec/head#' )),
            "  NumGABA sites = ", 
            len(moose.wildcardFind( '/model/chem/GABA/GABA[]' )) )
    '''
    runtime = RUNTIME

    if args.voltage_clamp:
        moose.element( "/model/stims/stim0" ).expr = gluR_clamp_potl
    moose.reinit()
    moose.start( runtime )
    offset = -90.0 if args.spiking else -70.0
    plot0 = moose.element( '/model/graphs/plot0' ).vector
    dt = moose.element( '/model/graphs/plot0' ).dt
    #t = np.arange( 0, len( plot0 ) * dt - 1e-6, dt )
    if args.voltage_clamp:
        moose.element( "/model/stims/stim0" ).expr = GABAR_clamp_potl
        moose.reinit()
        moose.start( runtime )
        plot0 = moose.element( '/model/graphs/plot0' ).vector
        plot0[:10] = GABAR_clamp_offset / 1e12  # clear out transient

    displayGraphs()

    moose.delete( "/model" )
    moose.delete( "/library" )
    return (plot0, patternIdx, args.repeatIdx, args.seed)

def displayGraphs():
    numPlots = len( plotList )
    '''
    iplots = moose.element( '/model/graphs/plot14' )
    moose.showfield( iplots )
    nd = iplots.numData
    iplots = moose.vec( '/model/graphs/plot14' )
    for idx, ii in enumerate(iplots):
        print( idx, " sum = ", sum( ii.vector ) )
    iplots = moose.vec( '/model/graphs/plot2' )
    for idx, ii in enumerate(iplots):
        print( idx, " sum2 = ", sum( ii.vector ) )
    '''
    t = np.arange( 0, RUNTIME+chemDt, chemDt )
    fig = plt.figure( figsize = (10, 15) )
    plt.rcParams.update( {"font.size": 16} )
    for idx, [compt, ii, path, field, label] in enumerate( plotList ): 
        if idx < 2:
            color = "black"
        elif idx < 7:
            color = "blue"
        else:
            color = "red"
        ax = plt.subplot( numPlots, 1, idx+1 )
        if idx == 0:
            ax.text( -0.05, 1.05, "C", fontsize = 20, weight = 'bold', transform=ax.transAxes )
        if idx == 1:
            dat= -moose.element( '/model/graphs/plot{}'.format(idx) ).vector
        elif idx == 6:
            dat = moose.element( '/model/graphs/plot6[37]'.format(idx) ).vector
        elif idx == 11:
            dat = moose.element( '/model/graphs/plot11[53]'.format(idx) ).vector
        else:
            dat = moose.element( '/model/graphs/plot{}'.format(idx) ).vector
        ax.plot( np.linspace(0, RUNTIME, len(dat) ), dat, color = color )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_yaxis().set_ticks([])
        if idx == 4:
            ax.text( -0.05, 0.15, label, fontsize = 16, weight = None, transform=ax.transAxes )
        else:
            ax.text( -0.05, 0.65, label, fontsize = 16, weight = None, transform=ax.transAxes )
        if idx < 11:
            ax.spines['bottom'].set_visible(False)
            ax.get_xaxis().set_ticks([])
        else:
            ax.set_xlabel( "Time (s)" )

    fig.tight_layout()
    plt.show()


def main():

    parser = argparse.ArgumentParser( description = "Deliver patterned stims to neuron with short-term-plasticity in E and I synapses" )
    parser.add_argument( "-n", "--numProcesses", type = int, help = "Number of processes to launch, default = 1", default = 1 )
    parser.add_argument( "-nr", "--numRepeats", type = int, help = "Number of repeats for each pattern, default = 1", default = 1 )
    parser.add_argument( "-p", "--pattern", type = int, help = "Index of pattern. 5-sq patterns are 46 to 50. 15 sq patterns are 52, 53, 55. Default = 46", default = 46 )
    parser.add_argument( '-spk', '--spiking', action="store_true", help ='Flag: when set, use high Na/K channel densities in soma to get spiking.' )
    parser.add_argument( '-v', '--voltage_clamp', action="store_true", help ='Flag: when set, do voltage clamp for glu and GABA currents respectively.')
    parser.add_argument( '-d', '--deterministic', action="store_true", help ='Flag: when set, use deterministic ODE solver. Normally uses GSSA stochastic solver.')
    parser.add_argument( "-m", "--modelName", type = str, help = "Optional: specify name of presynaptic model file, assumed to be in ./Models dir.", default = "BothPresyn75.g" )
    parser.add_argument( "-s", "--seed", type = int, help = "Optional: Seed to use for random numbers both for Python and for MOOSE.", default = 1234 )
    parser.add_argument( "-vglu", "--volGlu", type = float, help = "Optional: Volume scaling factor for Glu synapses. Default=2", default = 2.0 )
    parser.add_argument( "-vGABA", "--volGABA", type = float, help = "Optional: Volume scaling factor for GABA synapses. Default=0.5", default = 0.5 )
    parser.add_argument( "--pInter_CA1", type = float, help = "Optional: Probability of a given Interneuron connecting to the CA1 cell. Default=0.004 ", default = 0.004 )
    parser.add_argument( "--pCA3_CA1", type = float, help = "Optional: Probability of a given CA3 cell connecting to the CA1 cell. Default=0.002 ", default = 0.0002 )
    parser.add_argument( "--pCA3_Inter", type = float, help = "Optional: Probability of a given CA3 cell connecting to an interneuron. Default=0.002 ", default = 0.0008 )

    parser.add_argument( "-o", "--outputFile", type = str, help = "Optional: specify name of output file, in hdf5 format.", default = "simData.h5" )
    args = parser.parse_args()
    argdict = vars( args )
    argdict["pattern"] = 50
    argdict["repeatIdx"] = 0
    argdict["seed"] = args.seed + 5014 * args.numRepeats + ii
    innerArgs = argparse.Namespace( **argdict )
    innerMain( innerArgs )

    
if __name__ == "__main__":
    main()