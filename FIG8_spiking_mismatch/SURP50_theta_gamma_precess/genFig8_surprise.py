import moose
import pandas
import matplotlib.pyplot as plt
import numpy as np
import math
import rdesigneur as rd
import multiprocessing
from pathlib import Path
import argparse

freq = 80.0 # Hz
stimDuration = 0.002   # seconds
numPulses = 16
stimAmpl = 5e-2     # mM
basalCa = 0.08e-3   # mM
GABAdelay = 5.0e-3  # seconds
stimWidth = 0.002
firstPlotEntry = ['soma', '1', '.', 'Vm', 'Membrane potential']
numSq = 15

gluStimStr = "8e-5"
GABAStimStr = "8e-5"
gluR_clamp_potl = "-0.07"
GABAR_clamp_potl = "0.0"
GABAR_clamp_offset = 0.1    # nA
gluConductanceScale = 2   # Relative to default value in the spine proto
gluTau2Scale = 2   # Relative to default value in the spine proto

numCA1Exc = 100
numCA1Inh = 200
pCA3_CA1 = 0.0002
pCA3_Inter = 0.0008
pInter_CA1 = 0.004
#fracDesensitize = 1.2 # 1/4 of CA3 cells are active early but desensitize

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
pulseTrig = []
pulseThresh = 0.001
patternData = "../../../2022/VC_DATA/all_cells_SpikeTrain_CC_long.h5"
SAMPLE_FREQ = 20000
elecDt = 0.00005
chemDt = 0.0005
SAMPLE_TIME = 5
#SAMPLE_TIME = 2
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
patternDict2 = {}

## Here are params for the ChR2 desensitization
tauCell = 0.010         # Charging tau through gLeak
ChR2chanOpenTime = 0.001     # Balances out tauChR2chan
tauChR2chan = 0.005     # Charging tau through ChR2chan.
tauChR2recovery = 1.5   # Tau for recovery
ChR2decrement = 0.0004  # Scaling for decrement
ChR2_basal_desensitization = 0.01
Erest = 0               # Using baseline as zero.
EChR2 = 60              # Reversal potl in mV relative to baseline.
FracChR2active = np.zeros( NUM_SAMPLES // 10 )

def updatePulseTrain( freq ):
    ReducedPulseIdx = np.zeros( round(SAMPLE_TIME / chemDt ), dtype=int )
    ReducedPulseIdx[int(round(0.2 / chemDt) ) ] = 1

    for pp in range( 32 ):
        ReducedPulseIdx[int(round( (0.5 + pp/freq) / chemDt)) ] = 1

    return ReducedPulseIdx

def desensitization( events, dt ):
    # Turn off desensitization
    return np.ones_like( events ) * EChR2
    '''
    ret = []
    lastStim = -10.0    # Long time since last event.
    dy = ChR2_basal_desensitization
    Vm = Erest
    G = 1.0
    for idx, rr in enumerate( events ):
        t = idx * dt
        if rr:
            dy =  ChR2_basal_desensitization + ChR2decrement/(t-lastStim)
            #G -= dt * G*( 1-np.exp(-dy) )
            G -= G*( 1-np.exp(-dy) )
            Vm += ChR2chanOpenTime * ( EChR2-Vm ) * G / tauChR2chan
            lastStim = t
        G += dt * ( 1-G )/tauChR2recovery
        Vm += dt * (Erest-Vm)/tauCell
        ret.append( Vm )
    return np.array( ret )
    '''

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
             0,1,1,1,1,1,1,1,
             1,1,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,
             0,1,1,1,1,1,1,1,]

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
             1,1,1,1,1,1,1,1,
             1,1,1,1,1,1,1,1,
             1,1,0,0,0,0,0,0,
             1,1,1,1,1,1,1,0,
             1,1,1,1,1,1,1,0,
             1,1,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,]

    patternG =  [
             0,0,1,1,1,1,1,0,
             0,1,1,1,1,1,1,1,
             0,1,1,0,0,0,0,0,
             0,1,1,0,0,1,1,0,
             0,1,1,0,0,0,1,1,
             0,1,1,0,0,0,1,1,
             0,1,1,1,1,1,1,1,
             0,0,1,1,1,1,1,0,]

    patternH =  [
             1,1,0,0,0,1,1,0,
             1,1,0,0,0,1,1,0,
             1,1,0,0,0,1,1,0,
             1,1,1,1,1,1,1,0,
             1,1,1,1,1,1,1,0,
             1,1,0,0,0,1,1,0,
             1,1,0,0,0,1,1,0,
             1,1,0,0,0,1,1,0,]

    patternI =  [
             0,1,1,1,1,1,1,0,
             0,1,1,1,1,1,1,0,
             0,0,0,1,1,1,0,0,
             0,0,0,1,1,1,0,0,
             0,0,0,1,1,1,0,0,
             0,0,0,1,1,1,0,0,
             1,1,1,1,1,1,1,1,
             1,1,1,1,1,1,1,1,]

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

pandasColumnNames = [
                 'cellID',              'sex',         'ageAtInj',
              'ageAtExpt',       'incubation',             'unit',
               'location',         'protocol',          'exptSeq',
                 'exptID',            'sweep',         'stimFreq',
                  'numSq',        'intensity',       'pulseWidth',
              'clampMode',   'clampPotential',        'condition',
                     'AP',               'IR',              'tau',
          'sweepBaseline',      'numPatterns',      'patternList',
              'numPulses',  'pulseTrainStart',  'probePulseStart',
       'frameChangeTimes',       'pulseTimes',      'sweepLength',
           'baselineFlag',           'IRFlag',           'RaFlag',
            'spikingFlag',         'ChR2Flag',        'fieldData',
             'peaks_cell',  'peaks_cell_norm',         'auc_cell',
             'slope_cell',       'delay_cell',      'peaks_field',
       'peaks_field_norm',         'cell_fpr',        'field_fpr',
               'cell_ppr',        'cell_stpr',        'field_ppr',
             'field_stpr'
]

colIdx = { nn:idx for idx, nn in enumerate( pandasColumnNames )}

def makeRow( data, args ):
    TRIG = np.zeros( NUM_SAMPLES )
    TRIG[int( round( 0.2*SAMPLE_FREQ ) )] = 1.0
    for ii in range(50):
        idx = int( round (0.5 * SAMPLE_FREQ + ii*SAMPLE_FREQ/args.freq) )
        if not(ii < 4 or ii in range(9,17) or ii in range(22,30) or ii in range( 35,43 ) or ii > 47):
            TRIG[idx] = 1.0
            print( "TRIG idx = ", ii, idx )
        #print( ii, "TRIG = ", idx )
    print( "NUM TRIG = ", sum(TRIG) )
    row = [0]*SAMPLE_START + list( data[:NUM_SAMPLES] ) + [0.0]*NUM_SAMPLES + list(TRIG) + [0.0]*NUM_SAMPLES
    row[colIdx['exptSeq']] = 0
    row[colIdx["patternList"]] = [46,47,48,49]
    #row[colIdx["patternList"]] = [52,53,54,55]
    row[colIdx["numSq"]] = 5
    row[colIdx["sweep"]] = args.repeatIdx
    row[colIdx["stimFreq"]] = args.freq
    row[colIdx["clampMode"]] = "VC" if args.voltage_clamp else "CC"
    # Do clampPotential
    row[colIdx["intensity"]] = 100
    row[colIdx["protocol"]] = "surprise"

    return row

def evenOutConnectivity( connMtx, px, patternDict2 ):
    # patTemplate whose nonzero entries are the pattern # of unique pts
    patList = [46,47,48,49,50]
    #patList = [52,53,54,55]
    patTemplate = np.zeros( len( px[0].flat), dtype = int ) 
    for pp, qq in zip( px, patList):
        patTemplate += np.array(pp.flat, dtype = int ) * qq
    for idx, pp in enumerate( patTemplate ):
        if pp > 55:
            patTemplate[idx] = 0
    #print( "UNIQUE = ", np.unique( patTemplate ) )
    # for each pattern, form a list of the indices of unique entries.
    uniqueEntries = { pp:[idx for idx, tt in enumerate(patTemplate) if tt == pp] for pp in patList}
    #print( uniqueEntries ) 
    meanNumSyn = 0
    maxNumSyn = 0
    numSyn = {}
    for pp in patList:
        pat = patternDict2[pp]
        numSyn[pp] = sum( np.matmul( connMtx, pat ) > 0 )
        meanNumSyn += numSyn[pp]
        maxNumSyn = max( numSyn[pp], maxNumSyn )
        print( "pat {}: len = {}".format( pp, numSyn[pp] )  )
    meanNumSyn = int( meanNumSyn / len( patList ) ) # Go for the mean.
    print( "Mean = ", meanNumSyn, "Max = ", maxNumSyn )
    for pp in patList:
        uu = uniqueEntries[pp]
        for qq in range( maxNumSyn - numSyn[pp] ):
        #for qq in range( meanNumSyn - numSyn[pp] ):
            if qq < len( uu ):
                # Now we need to find untouched synapses for this pixel
                syns = connMtx[:,uu[qq]].flat
                #print( "doing", qq, uu[qq], len(syns) )
                zz = np.arange( len(syns), dtype=int)[syns == 0 ]
                if qq < len( zz ):
                    syns[zz[qq]] = 1
    for pp in patList:
        pat = patternDict2[pp]
        ret = sum( np.matmul( connMtx, pat ) > 0 )
        print( "After fix: pat {}: len = {}".format( pp, ret )  )


def generatePatterns( args ):
    global CA3_CA1
    global CA3_Inter
    global Inter_CA1
    global patternDict2
    pd = patternDict()

    # Assumes all are random projections.
    np.random.seed( args.seedConnections )
    CA3_Inter = (np.random.rand( 256, 256 ) < args.pCA3_Inter) * 1.0
    CA3_CA1 = (np.random.rand( numCA1Exc, 256 ) < args.pCA3_CA1) * 1.0
    Inter_CA1 = (np.random.rand( numCA1Inh, 256 ) < args.pInter_CA1) * 1.0
    px = []
    for char in ["A", "B", "C", "D", "E"]:
        #temp = pd[char].reshape(8,8).repeat(2, axis = 0)
        temp = np.array(pd[char]).reshape(8,8).repeat(4, axis=0).reshape(-1)
        # Now put in a mask that zeros 3/4 of the values
        zero_indices = np.random.choice(256, args.zeroIndices,replace=False)
        temp[zero_indices] = 0
        px.append( temp )
        #px.append( np.append(temp2.reshape(-1), np.zeros(128) ) )

    patternDict2 = {
        0:np.zeros(256),
        46:px[0],
        47:px[1],
        48:px[2],
        49:px[3],
        50:px[4],
        51:px[0].transpose(),
        52:px[1].transpose(),
        53:px[2].transpose(),
        54:px[3].transpose(),
        55:px[4].transpose(),
    }
    print( "Even out CA3_CA1::" )
    evenOutConnectivity( CA3_CA1, px, patternDict2 )
    print( "\nEven out CA3_Inter::" )
    evenOutConnectivity( CA3_Inter, px, patternDict2 )

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


#def buildModel( presynModelName, seed, useGssa, vGlu, vGABA, spiking ):
def buildModel( args ):
    useGssa = not args.deterministic
    rGlu = pow( args.volGlu, 1.0/3.0)
    rGABA = pow( args.volGABA, 1.0/3.0)
    NaGbar = 400.0 if args.spiking else 6.0
    KGbar = 450.0 if args.spiking else 3.5

    rdes = rd.rdesigneur(
        elecDt = elecDt,
        chemDt = chemDt,
        funcDt = 0.0005,
        diffDt = 0.0005,
        chemPlotDt = 0.0005,
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
            ['Models/{}'.format( args.modelName ), 'chem'],
        ],
        chanProto = [
            ['make_Na()', 'Na'],
            ['make_K_DR()', 'K'],
            ['make_GABA()', 'GABA']
        ],
        # Assign passive properties of the cell.
        passiveDistrib = [
            ['#', 'CM', '0.01', 'Em', '-0.065', 'RM', '1.0']
        ],
        chemDistrib = [
            # Args: chem_model, elec_compts, mesh_type, spatial_distrib, r_scalefactor, radius_sdev
            ['glu', 'head#', 'presyn_spine', '1', rGlu, 0.5 ],
            # Args: chem_model, elec_compts, mesh_type, spatial_distrib, r_absolute, r_sdev, spacing
            ['GABA', 'dend#', 'presyn_dend', '1', 0.5e-6 * rGABA, 0, 1e-6 ],
        ],
        chanDistrib = [
            ['Na', 'soma', 'Gbar', str(NaGbar) ],
            ['K', 'soma', 'Gbar', str(KGbar) ],
            ['GABA', 'dend#', 'Gbar', str( args.wtGABA ) ]
        ],
        adaptorList = [
            [ 'glu/glu', 'n', 'glu', 'activation', 0.0, 1500 ],
            [ 'GABA/GABA', 'n', 'GABA', 'activation', 0.0, 1500 ],
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
        ]
    )
    moose.seed( args.seedStochastic ) 
    gluReceptor = moose.element( '/library/spine/head/glu' )
    gluReceptor.Gbar *= args.wtGlu # Tweak conductance
    #print ( "GLUGGGGG = ", vGlu, gluConductanceScale )
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

def isPulse( freq, t, precession ):
    if t < 0.1999:
        return False, 0
    elif t < (0.1999 + stimWidth):
        return True, 0
    pulseNum = int( round( (t - 0.5) * freq ) )
    if pulseNum > 49:
        return False, 0
    pulseT = 0.49999 + pulseNum / freq
    if ( t > 0.49999 and t > pulseT and t < (pulseT + stimWidth )) :
        #print( "pulseNum ={}, t = {:.4f}, freq={:.4f}, patIdx = {}".format( pulseNum, t, freq, pulseNum//8 ) )
        if pulseNum < 4 or pulseNum in range(9,17) or pulseNum in range(22,30) or pulseNum in range( 35,43 ) or pulseNum > 47:
            return True, 0
        else: 
            #return True, 46 + pulseNum // 12 # This increments each theta
            # Precession of 13 keeps it fixed, 12 advances, 14 retreats.
            # There are only 10 patterns.
            return True, 46 + ((pulseNum-4) % precession)%10
    return False, 46


def stimFunc( freq, ChR2AmplScale, precession ):
    t = moose.element( '/clock' ).currentTime
    stimWidthIdx = int( round( stimWidth / chemDt ) )
    CA3isActive, patternIdx = isPulse( freq, t, precession )
    InterIsActive, patternIdx2 = isPulse( freq, t - GABAdelay, precession )
    idx = int(round( t/chemDt ) )
    if idx >= len( FracChR2active ):
        return
    # Stim is to be delivered for the entire duration of stimWidth.
    #if idx >= stimWidthIdx:
    if CA3isActive:
        chr2Ampl = max(FracChR2active[1+idx-stimWidthIdx:idx+1]) * ChR2AmplScale
        #print( "chr2amp if idx > stimwidthIdx: ", chr2Ampl, idx )
    else:
        chr2Ampl = FracChR2active[idx] * ChR2AmplScale
        #print( "chr2amp if idx <= stimwidthIdx: ", chr2Ampl, idx )
    idx2 = int( round( (t - GABAdelay) / chemDt ) )
    if idx2 >= len( FracChR2active ):
        return
    #if idx2 > stimWidthIdx:
    if InterIsActive:
        chr2Ampl2 = max(FracChR2active[1+idx2-stimWidthIdx:1+idx2]) * ChR2AmplScale
    else:
        chr2Ampl2 = FracChR2active[idx2] * ChR2AmplScale
    gluInput = moose.vec( "/model/chem/glu/Ca_ext" )
    gabaInput = moose.vec( "/model/chem/GABA/Ca_ext" )
    if CA3isActive:
        #inputPattern = patternDict2[ patternIdx]
        print( "Trig CA3 at {:.3f} {} with {} @ {:.3f}".format( t, idx, patternIdx, chr2Ampl ))
        ca3cells = moose.vec( "/model/elec/CA3/soma" )
        pd = patternDict2[patternIdx]
        amplIdx = min( len( pd ), int( chr2Ampl * len( pd ) ) )
        pd =  np.append( pd[:amplIdx], np.zeros( len(pd)-amplIdx ) )
        ca3cells.Vm = pd
        '''
        #ca3cells.Vm = np.where( np.random.rand(len(pd)) < chr2Ampl, pd, 0 )
        amplIdx = min( len( pd ), int( chr2Ampl * len( pd ) ) )
        ca3cells.Vm = np.append( pd[:amplIdx], np.zeros( len(pd)-amplIdx ) )
        '''

        gluInput.concInit = (np.matmul( CA3_CA1, ca3cells.Vm ) >= thresh_CA3_CA1 ) * stimAmpl
        print( "{}  t={:.5f} idx={}  NUMGlu={:.1f}    chr2Ampl={:.3f}".format( patternIdx, t, idx, sum( gluInput.concInit ) / stimAmpl, chr2Ampl ), flush=True )
    else:
        moose.vec( "/model/elec/CA3/soma" ).Vm = 0.0
        gluInput.concInit = basalCa

    if InterIsActive:
        pd = patternDict2[patternIdx]
        amplIdx = min( len( pd ), int( chr2Ampl2 * len( pd ) ) )
        #pd = np.where( np.random.rand( len( pd ) ) < chr2Ampl, pd, 0 )
        pd =  np.append( pd[:amplIdx], np.zeros( len(pd)-amplIdx ) )

        Inter = moose.vec( "/model/elec/Inter/soma" )
        Inter.Vm = (np.matmul( CA3_Inter, pd) >= thresh_CA3_Inter ) * 1.0
        gabaInput.concInit = (np.matmul( Inter_CA1, Inter.Vm ) >= thresh_Inter_CA1 ) * stimAmpl
        if patternIdx == 0:
            print( "{:.4f} {}  NUMGABA={:.1f}    chr2Ampl={:.3f}".format( t, patternIdx, sum( gabaInput.concInit) / stimAmpl, chr2Ampl ), flush=True )
    else:
        gabaInput.concInit = basalCa
        moose.vec( "/model/elec/Inter/soma" ).Vm = 0
    #print( "{:.4f}  CA3={}, Inter={}, pattern={}".format( t, CA3isActive, InterIsActive, patternIdx))

def makeNetwork( rdes ):
    origNeuronId = rdes.elecid
    CA3cells, CA3args = makeInputs( "CA3", 20e-6 )
    interneurons, interneuronArgs = makeInputs( "Inter", 70e-6 )

def innerMain( args, ReducedPulseIdx ):
    global FracChR2active
    FracChR2active = desensitization( ReducedPulseIdx, chemDt )
    patternIdx = args.pattern
    if args.voltage_clamp:
        stimList = [['soma', '1', '.', 'vclamp', '-0.070' ]]
        firstPlotEntry = ['soma', '1', 'vclamp', 'current','Vclamp current']
        # Possibly here we need to increase the vclamp gain.
    else:
        firstPlotEntry = ['soma', '1', '.', 'Vm', 'Membrane potential']

    generatePatterns( args )

    rdes = buildModel( args )
    pr = moose.PyRun( "/model/stims/stimRun" )
    pr.runString = 'stimFunc({}, {}, {})'.format( args.freq, args.ChR2_ampl, args.precession )
    pr.tick = 14 # This would be chemDt. Which is currently 0.5 ms.

    makeNetwork( rdes )

    moose.reinit()
    numGluRR = moose.vec('/model/chem/glu/RR_pool').nInit
    numGABARR = moose.vec('/model/chem/GABA/RR_pool').nInit
    print( "NumGluRR = {}, mean = {:.3f}, sdev = {:.3f}".format( 
        len( numGluRR ), np.mean(numGluRR), np.std( numGluRR ) ) )
    print( "NumGABARR = {}, mean = {:.3f}, sdev = {:.3f}".format ( 
        len( numGABARR ), np.mean(numGABARR), np.std( numGABARR ) ) )
    '''
    print( "NumGluR = ", len( moose.vec('/model/chem/glu/glu') ),
        "  NumGABA sites = ", 
        len(moose.vec( '/model/chem/GABA/GABA' )) )
    print( "NumSpines = ", len(moose.wildcardFind( '/model/elec/head#' )),
            "  NumGABA sites = ", 
            len(moose.wildcardFind( '/model/chem/GABA/GABA[]' )) )
    '''
    runtime = SAMPLE_TIME

    if args.voltage_clamp:
        moose.element( "/model/stims/stim0" ).expr = gluR_clamp_potl
    # This time we reseed it so that the stoch runs differ.
    moose.seed( args.seedStochastic + args.freq * args.numRepeats + args.repeatIdx )
    moose.reinit()
    moose.start( runtime )
    offset = -90.0 if args.spiking else -70.0
    plotE = np.zeros( NUM_SAMPLES )
    plot0 = moose.element( '/model/graphs/plot0' ).vector
    for ee in moose.vec( '/model/graphs/plot1' ):
        plotE += ee.vector[:NUM_SAMPLES]
    #print( "numSpinePlots = ", len( moose.vec( '/model/graphs/plot1' ) ) )
    #plotE = moose.vec( '/model/graphs/plot1' )[gluSynIdx].vector
    plotI = moose.vec( '/model/graphs/plot2' )[0].vector # Only one dend rec
    dt = moose.element( '/model/graphs/plot0' ).dt
    #t = np.arange( 0, len( plot0 ) * dt - 1e-6, dt )
    if args.voltage_clamp:
        moose.element( "/model/stims/stim0" ).expr = GABAR_clamp_potl
        moose.reinit()
        moose.start( runtime )
        plot0 = moose.element( '/model/graphs/plot0' ).vector
        plot0[:10] = GABAR_clamp_offset / 1e12  # clear out transient

    moose.delete( "/model" )
    moose.delete( "/library" )
    return (plot0, args.freq, args.repeatIdx )

def runSession( args, whichArg ):
    changedValue = 0
    if whichArg != "orig":
        changedValue = getattr(args, whichArg)
    fname = "{}_{}_{}.h5".format( Path( args.outputFile ).stem, whichArg, changedValue )
    print( "Working on: ", fname )
    pool = multiprocessing.Pool( processes = args.numProcesses )
    ret = []
    data = []
    argdict = vars( args )
    for freq in [100]:
        ReducedPulseIdx = updatePulseTrain( freq )
        for ii in range( args.numRepeats ):
            argdict["repeatIdx"] = ii
            argdict["freq"] = freq
            argdict["seedConnections"] = args.seedConnections
            argdict["seedStochastic"] = args.seedStochastic
            print( "Launching {}.{}".format( freq, ii ) )
            innerArgs = argparse.Namespace( **argdict )
            ret.append( pool.apply_async( innerMain, args = (innerArgs, ReducedPulseIdx )))
    for rr in ret:
        ( plot0, args.freq, args.repeatIdx ) = rr.get()
        data.append( makeRow( plot0*1e3, args ) )
    df = pandas.DataFrame(data, columns=pandasColumnNames + [str(i) for i in range( NUM_SAMPLES *4 ) ] )
    #df.to_hdf( args.outputFile, "SimData", mode = "w", format='table', complevel=9)
    df.to_hdf( fname, "SimData", mode = "w", complevel=9)

def main():
    global freq
    global firstPlotEntry
    global stimList
    global numPulses
    global numSq
    global pulseTrig
    global pulseThresh

    parser = argparse.ArgumentParser( description = "Deliver patterned stims to neuron with short-term-plasticity in E and I synapses" )
    parser.add_argument( "-n", "--numProcesses", type = int, help = "Number of processes to launch, default = 1", default = 1 )
    parser.add_argument( "-nr", "--numRepeats", type = int, help = "Number of repeats for each pattern, default = 1", default = 1 )
    parser.add_argument( "-p", "--pattern", type = int, help = "Index of pattern. 5-sq patterns are 46 to 50. 15 sq patterns are 52, 53, 54, 55. Default = 46", default = 46 )
    parser.add_argument( '-spk', '--spiking', action="store_true", help ='Flag: when set, use high Na/K channel densities in soma to get spiking.' )
    parser.add_argument( '-v', '--voltage_clamp', action="store_true", help ='Flag: when set, do voltage clamp for glu and GABA currents respectively.')
    parser.add_argument( '-d', '--deterministic', action="store_true", help ='Flag: when set, use deterministic ODE solver. Normally uses GSSA stochastic solver.')
    parser.add_argument( "-m", "--modelName", type = str, help = "Optional: specify name of presynaptic model file, assumed to be in ./Models dir.", default = "BothPresyn86.g" )
    parser.add_argument( "-s", "--seed", type = int, help = "Optional: Seed to use for random numbers both for Python and for MOOSE.", default = 1234 )
    parser.add_argument( "-sc", "--seedConnections", type = int, help = "Optional: Seed to use for random numbers for setting up connections in Python.", default = 18 )
    parser.add_argument( "-ss", "--seedStochastic", type = int, help = "Optional: Seed to use for random numbers for MOOSE stochastic calculations.", default = 1234 )
    parser.add_argument( "-vglu", "--volGlu", type = float, help = "Optional: Volume scaling factor for Glu synapses. Default=0.5", default = 0.5 )
    parser.add_argument( "-wglu", "--wtGlu", type = float, help = "Optional: weight scaling factor for Glu synapses. Default=3.0", default = 3.0 )
    parser.add_argument( "-vGABA", "--volGABA", type = float, help = "Optional: Volume scaling factor for GABA synapses. Default=0.5", default = 0.5 )
    parser.add_argument( "-wGABA", "--wtGABA", type = float, help = "Optional: Weight of GABA synapses. Default=10", default = 10 )
    parser.add_argument( "--pInter_CA1", type = float, help = "Optional: Probability of a given Interneuron connecting to the CA1 cell. Default=0.01 ", default = 0.01 )
    parser.add_argument( "--pCA3_CA1", type = float, help = "Optional: Probability of a given CA3 cell connecting to the CA1 cell. Default=0.02 ", default = 0.02 )
    parser.add_argument( "--pCA3_Inter", type = float, help = "Optional: Probability of a given CA3 cell connecting to an interneuron. Default=0.01 ", default = 0.01 )
    parser.add_argument( "--ChR2_ampl", type = float, help = "Optional: Scale factor for ChR2 stimulus amplitude. Default=0.1", default = 0.1 )
    parser.add_argument( "--precession", type = int, help = "Optional: Determines precession of gamma patterns wrt theta. Default of 13 does not precess.", default = 13 )
    parser.add_argument( "-z", "--zeroIndices", type = int, help = "Optional: Number of optical inputs to zero out, range 0 to 256. Default=192.", default = 192 )

    parser.add_argument( "-o", "--outputFile", type = str, help = "Optional: specify name of output file, in hdf5 format.", default = "P.h5" )
    args = parser.parse_args()
    orig = args.precession
    for args.precession in [13, 12, 14]:
        runSession( args, "precession" )
    args.precession = orig
    
if __name__ == "__main__":
    main()
