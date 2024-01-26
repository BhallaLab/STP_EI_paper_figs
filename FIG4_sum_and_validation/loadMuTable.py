import pandas
import numpy as np
import matplotlib.pyplot as plt
import moose

#patternData = "/home1/bhalla/adityaa/Lab/Projects/EI_Dynamics/Analysis/parsed_data/all_cells_SpikeTrain_VC_long.h5"
patternData = "../../../2022/VC_DATA/all_cells_SpikeTrain_VC_long.h5"

SAMPLE_FREQ = 20000
SAMPLE_TIME = 11
NUM_SAMPLES = SAMPLE_FREQ * SAMPLE_TIME
SAMPLE_START = 49
SWEEP = 39

def main():
    df = pandas.read_hdf( patternData )
    # 48 rows x 880049 columns
    stimTTL = df.iloc[SWEEP, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ]
    stimTrig = np.array(stimTTL[::10])
    print( "Trig = ", stimTrig, sum( stimTrig ) )
    quit()
    # 2 ms pulses, sample rate 20KHz. pk=0.5
    #print( df )
    print( df.columns[:55] )
    print( "LEN COL = ", len( df.columns ) )
    #dataStartColumn = len( df.columns ) - 80000
    print("PL = \n", df['patternList'])
    print("numpat = \n", df['numPatterns'])
    print("numpulse = \n", df['numPulses'])
    print("pulseTimes = \n", df['pulseTimes'])
    #print( df.iloc[:,7:33])
    print( df.iloc[:,:22])
    t = np.linspace( 0, SAMPLE_TIME, NUM_SAMPLES, endpoint = False )
    plt.figure( figsize = ( 15, 12 ) )
    ax = plt.subplot(4, 1, 1)
    ax.plot( t, df.iloc[SWEEP, SAMPLE_START:SAMPLE_START + NUM_SAMPLES ] )
    ax.set_title( "Block 1: EPSC?" )

    ax = plt.subplot(4, 1, 2)
    ax.plot( t, df.iloc[SWEEP, SAMPLE_START + NUM_SAMPLES:SAMPLE_START+2*NUM_SAMPLES ] )
    ax.set_title( "Block 2: Trigger?" )


    ax = plt.subplot(4, 1, 3)
    ax.plot( t, df.iloc[SWEEP, SAMPLE_START + 2*NUM_SAMPLES:SAMPLE_START+3*NUM_SAMPLES ] )
    ax.set_title( "Block 3: TTL? Doesn't look like photodiode" )

    ax = plt.subplot(4, 1, 4)
    ax.plot( t, df.iloc[SWEEP, SAMPLE_START + 3*NUM_SAMPLES: ] )
    ax.set_title( "Block 4: Field??" )


    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
