import matplotlib.pyplot as plt
import numpy as np
import math
from pathlib import Path
from scipy.ndimage import zoom
from scipy.ndimage import shift
import argparse



proj = []
freq = 80.0 # Hz
stimDuration = 0.002   # seconds
numPulses = 16
CaStimAmpl = 2e-2     # mM. This is 20 micromolar.
basalCa = 0.08e-3   # mM
GABAdelay = 5.0e-3  # seconds
width = 0.002
numSq = 15
ChR2_sigma = 2.0  # Width of uniform deviate multiple for ChR2 stimulus
ChR2_background = 0.05 # Width of uniform deviate for ChR2 stimulus


seed = 1234
numCA3 = 256
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
CA3_noise = 0.0
Inter_noise = 0.0
CA3_spk_thresh = 10
#CA3_thresholded = np.zeros( numCA3 )

repeatPatterns = False
inputs = []
stimList = []
pulseTrig = []
pulseThresh = 0.001
SAMPLE_FREQ = 20000
elecDt = 0.00005
chemDt = 0.0005
SAMPLE_TIME = 11
#SAMPLE_TIME = 2
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 16
patternDrivenCellActivity = {}


def patternDict():
    patternZeros = [0]*64
    patternOnes = [1]*64

    patternA =  [
             0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,
             0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,
             0,0,1,0,0,0,1,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,0,0,]

    patternB =  [
             0,0,0,0,0,0,0,0,
             0,1,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,0,0,
             0,0,0,0,0,0,0,0,
             0,1,0,0,0,0,1,0,
             0,0,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,]

    patternC =  [
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,
             0,1,0,0,0,0,0,0,]

    patternD =  [
             0,0,0,0,0,0,0,0,
             0,0,1,0,0,1,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,0,
             0,1,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,1,0,0,]

    patternE =  [
             0,1,0,0,0,0,0,0,
             0,0,0,0,0,0,1,0,
             0,0,1,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,0,0,
             0,0,0,0,0,0,0,0,]

    patternF =  [ 
             0,0,0,0,0,0,1,0,
             0,1,0,1,0,0,0,0,
             0,0,0,1,0,1,0,0,
             1,0,0,0,1,0,0,0,
             0,0,1,0,0,0,1,0,
             0,1,0,0,0,0,1,0,
             0,0,0,1,0,0,1,1,
             0,1,0,0,0,0,0,0,]

    patternG =  [
             0,0,0,0,0,0,1,0,
             0,1,1,0,0,1,0,0,
             0,0,0,1,0,0,0,0,
             1,0,0,0,1,0,1,0,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,0,1,0,
             0,0,0,1,0,0,0,1,
             0,1,0,0,0,1,0,0,]

    patternH =  [
             0,1,0,0,0,0,1,0,
             0,0,1,0,0,1,1,0,
             0,0,1,1,0,0,0,0,
             1,0,0,0,0,0,1,0,
             1,1,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,0,1,
             0,1,0,0,0,1,0,0,]

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


def opticProjDict( patternDict2 ):
    global proj
    # First, we define a few template optical projections. These are
    # the somatic responses to stimuli at the specified points in space
    proj = [
        [
            0,0,0,0,1,0,0,0,
            0,0,0,1,2,0,1,0,
            0,0,0,2,0,0,0,0,
            0,0,1,3,1,1,0,0,
            0,0,0,9,1,1,0,0,
            0,0,6,4,1,0,0,0,
            0,0,2,2,1,0,0,0,
            0,0,3,2,2,1,1,0,
        ],
        [
            0,0,0,0,1,0,0,0,
            0,0,0,1,1,1,0,0,
            0,0,0,2,2,0,0,0,
            0,0,1,3,2,0,0,0,
            0,0,2,9,1,1,0,0,
            0,0,2,5,3,2,0,0,
            0,0,1,4,2,1,0,0,
            0,0,1,2,2,2,0,0,
        ],
        [
            0,0,0,1,0,1,1,0,
            0,0,0,1,0,0,1,0,
            0,0,0,2,1,2,0,0,
            0,0,1,9,7,0,0,0,
            0,0,2,3,4,1,0,0,
            0,0,1,2,3,1,0,0,
            0,0,1,1,2,1,0,0,
            0,0,1,1,0,0,0,0,
        ],
        [
            0,0,0,0,0,1,0,0,
            0,0,2,0,1,0,1,0,
            0,0,2,0,1,1,1,0,
            0,0,1,3,1,1,0,0,
            0,0,0,9,0,0,0,0,
            0,0,0,2,6,0,0,0,
            0,0,0,2,2,0,0,0,
            0,0,0,1,4,1,1,0,
        ]
    ]
    '''
    # Double the number of templates by getting mirror images.
    for pp in proj:
        proj2.append( np.array( pp ).reshape(8,8) )
        proj2.append( np.fliplr( proj2[-1] ) )
    # Expand these 8x8 arrays into 16x16 ones.
    '''

    proj2 = []
    # Zoom in each template neuron to 16x16, and add mirror image of each.
    for pp in proj:
        proj2.append( zoom( np.array( pp ).reshape(8,8), 2, order=1 ) )
        proj2.append( np.fliplr( proj2[-1] ) )
        #proj2.append( zoom( np.array( pp ).reshape(8,8), 2, order=1 ) )

    # Provide offsets to move these templates around in a bigger 16x16 grid.
    # We have above 8 templates. We need to shuffle these 32 times
    # to get a total of 256 CA3 neurons.
    proj3 = []
    for pp in proj2:
        for dx in [-7, -5, -3, -1, 1, 3, 5, 7]:
            for dy in [ -3, -1, 1, 3]:
                proj3.append( shift(pp, shift=[dx, dy], mode='constant', cval=0) )

    # Add noise to each of the templates above to get the full list.
    for idx, pp in enumerate( proj3):
        pp = pp*np.random.uniform( size=(16,16) ) * ChR2_sigma + np.random.uniform( size=(16,16) ) * ChR2_background
        #print( idx, np.mean( pp ) )

    # multiply and sum each input pattern with the template for each neuron,
    # aka a dot product, to get an input-output matrix
    patternStim = {}
    for pat, val in patternDict2.items():
        patternStim[pat] = np.array([ np.sum( val * pp.flatten() ) for pp in proj3 ])

    return patternStim


def generatePatterns():
    global CA3_CA1
    global CA3_Inter
    global Inter_CA1
    global patternDrivenCellActivity
    pd = patternDict()

    # Assumes all are random projections.
    np.random.seed( seed )
    CA3_Inter = (np.random.rand( numCA1Inh, numCA3 ) < pCA3_Inter) * 1.0
    CA3_CA1 = (np.random.rand( numCA1Exc, numCA3 ) < pCA3_CA1) * 1.0
    Inter_CA1 = (np.random.rand( numCA1Inh, numCA1Inh ) < pInter_CA1) * 1.0
    px = {}
    for char in ["A", "B", "C", "D", "E", "F", "G", "H", "I"]:
        orig = pd[char].reshape(8,8)
        temp = np.zeros( (16,16), dtype = float )
        temp[::2, ::2] = orig       # Python slicing is cool.
        #print( char, sum( orig ), sum( temp ) )
        px[char] = temp.flatten()

    patternDict2 = {
        46:px["A"],
        47:px["B"],
        48:px["C"],
        49:px["D"],
        50:px["E"],
        52:px["F"],
        53:px["G"],
        55:px["H"],
    }

    patternDrivenCellActivity = opticProjDict( patternDict2 )

def stimFunc( patternIdx, ChR2AmplScale ):
    global history
    #global CA3_thresholded
    t = moose.element( '/clock' ).currentTime
    # Need to look up if this is time to generate pulse. 
    idx = int(round( t/chemDt ) )
    '''
    if idx % int( 1.0/chemDt )  == 0:   # dot every second.
        print( ".", flush = True, end = "" )
    '''
    if idx >= len( ReducedPulseIdx ):
        return
    CA3isActive = (ReducedPulseIdx[idx] > 0.5) #Stimulus is to be delivered
    assert( len( FracChR2active ) == len( ReducedPulseIdx ) )
    chr2Ampl = FracChR2active[idx] * ChR2AmplScale
    idx2 = int( round( (t - GABAdelay) / chemDt ) )
    if idx2 >= len( ReducedPulseIdx ):
        return
    InterIsActive = ( ReducedPulseIdx[idx2] > 0.5 )
    gluInput = moose.vec( "/model/chem/glu/Ca_ext" )
    gabaInput = moose.vec( "/model/chem/GABA/Ca_ext" )
    if t == chemDt:
        print( "{}: AmplScale = {:5.2f}  ".format( patternIdx, ChR2AmplScale ) )
    if CA3isActive:
        history.CA3ActiveTime = t
        #print( "Trig CA3 at {:.3f} {} with {}".format( t, idx, patternIdx ))
        # Here I look up the activity of all cells, scale by chr2Ampl for
        # desensitization, add some noise to each cell, and decide if the
        # result of all this is over threshold for each cell.
        pdca = patternDrivenCellActivity[patternIdx]
        history.CA3_thresholded = (
            (pdca*chr2Ampl +
            np.random.normal( size=len(pdca), loc=0, scale=CA3_noise)
            ) > CA3_spk_thresh ) * 1
        #print( "crampl={:.3f}, sumglu={:.3f}, #ca3cells={:.2f}, pdca={:.2f}".format( chr2Ampl, sum( gluInput.concInit ), np.sum( CA3_thresholded ), max( pdca ) ) )
        if idx in [4, 1068, 2056, 3122, 4049, 4050]: # stimtrig values.
            print( "{:7.2f}, {:<3d}".format( chr2Ampl, int(np.sum(history.CA3_thresholded)) ), end = "", flush = True )

        if idx == 4050:    # end it.
            print( "  std={:.2f}  #zero={:d}".format(np.std(history.CA3_thresholded), np.sum(history.CA3_thresholded < 0.5) ) )

    if (t - history.CA3ActiveTime) < (stimDuration - 1e-6):
        #print( "CA1Stim", t )
        ca3cells = moose.vec( "/model/elec/CA3/soma" )
        ca3cells.Vm = history.CA3_thresholded 
        gluInput.concInit = (np.matmul( CA3_CA1, history.CA3_thresholded ) >= thresh_CA3_CA1 ) * CaStimAmpl
    else:
        moose.vec( "/model/elec/CA3/soma" ).Vm = 0.0
        gluInput.concInit = basalCa

    if InterIsActive:
        history.InterActiveTime = t
        history.Inter_thresholded = (
            (
                np.matmul( CA3_Inter, history.CA3_thresholded) + 
                np.random.normal(
                    size=len(history.Inter_thresholded), loc=0, scale=Inter_noise)
            ) >= thresh_CA3_Inter ) * 1.0

def panelB():
    patternA = patternDict()["A"]
    patternA_matrix = np.array(patternA).reshape(8, 8)
    # Reshape proj and patternA into 8x8 matrices
    pm = [np.array(proj[i]).reshape(8, 8) for i in range(len(proj))]
    # Shift proj entries by -3, -1, 1, 3, padding with zeros.
    m0 = np.pad(pm[0][:, 3:], ((0, 0), (0, 3)), mode='constant')
    m1 = np.pad(pm[1][:, 1:], ((0, 0), (0, 1)), mode='constant')
    m2 = np.pad(pm[2][:, :-1:], ((0, 0), (1, 0)), mode='constant')
    m3 = np.pad(pm[3][:, :-3:], ((0, 0), (3, 0)), mode='constant')
    pm2 = [m0, m1, m2, m3]
    # Plot the heatmaps
    plt.figure(figsize=(16, 4))
    for i in range(len(pm2)):
        plt.subplot(1, len(pm2), i+1)
        plt.imshow(pm2[i], cmap='hot', interpolation='nearest')
        # Find non-zero indices in patternA_matrix and draw a square open box
        for x in range(8):
            for y in range(8):
                if patternA_matrix[x, y] != 0:
                    rect = plt.Rectangle((y - 0.5, x - 0.5), 1, 1, fill=False, edgecolor='cyan', linewidth=4)
                    plt.gca().add_patch(rect)
        plt.title(f'Projection {i+1}')
        plt.axis('off')
    
    plt.tight_layout()
    #plt.show()


def panelC():
    patternA = patternDict()["A"]
    patternH = patternDict()["H"]
    patternA_matrix = np.array(patternA).reshape(8, 8)
    patternH_matrix = np.array(patternH).reshape(8, 8)

    aOut = patternDrivenCellActivity[46].reshape( 16, 16 )
    print( "aout shape = ", aOut.shape )
    hOut = patternDrivenCellActivity[52].reshape( 16, 16 )

    aOutThresh = (aOut > CA3_spk_thresh )*1
    hOutThresh = (hOut > CA3_spk_thresh )*1

    # Plot the heatmaps
    plt.figure(figsize=(8, 12))
    plt.subplot( 3, 2, 1 )
    plt.imshow(patternA_matrix, cmap='hot', interpolation='nearest')
    plt.subplot( 3, 2, 2 )
    plt.imshow(patternH_matrix, cmap='hot', interpolation='nearest')
    plt.subplot( 3, 2, 3 )
    plt.imshow(aOut, cmap='hot', interpolation='nearest')
    plt.subplot( 3, 2, 4 )
    plt.imshow(hOut, cmap='hot', interpolation='nearest')
    plt.subplot( 3, 2, 5 )
    plt.imshow(aOutThresh, cmap='hot', interpolation='nearest')
    plt.subplot( 3, 2, 6 )
    plt.imshow(hOutThresh, cmap='hot', interpolation='nearest')
    plt.axis('off')
    #plt.title(f'Projection {i+1}')
    plt.tight_layout()


def main():
    generatePatterns()
    panelB()
    panelC()
    plt.show()

    
if __name__ == "__main__":
    main()

