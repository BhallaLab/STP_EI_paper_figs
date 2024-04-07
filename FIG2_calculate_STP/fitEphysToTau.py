import moose
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sci
import rdesigneur as rd
gluCa = 5e-3   # 50 uM.
chemDt = 1e-4
chemPlotDt = 2e-4
elecDt = 50e-6
diffusionLength = 1e-3 # ensure single chem compt
settleTime = 0.050
stimTime = 0.002
postStimTime = 0.050
longestISI = 0.5
eps = 1e-6
runTime = settleTime + longestISI + postStimTime

def buildModel( tau ):
    #RM = str( tau / 10.0 )
    RM = str(1.0)
    #RA = str( 1.0 )
    RA = RM
    rdes = rd.rdesigneur(
        turnOffElec = False,
        verbose = False,
        chemDt = chemDt,
        chemPlotDt = chemPlotDt,
        elecDt = elecDt,
        elecPlotDt = elecDt,
        # Tau is in ms.
        #diffusionLength = 1e-3, # Default diffusion length is 2 microns
        cellProto = [['ballAndStick', 'soma',  10e-6, 10e-6, 2e-6, 200e-6, 1]],
        passiveDistrib = [
            ['#', 'CM', '0.01', 'Em', '-0.065', 'RM', RM, 'RA', RA],
        ],
        chanProto = [
            ['make_glu()', 'glu'],
            ['make_GABA()', 'GABA'],
        ],
        chanDistrib = [
            ['glu', 'dend#', 'Gbar', '2'],
            ['GABA', 'dend#', 'Gbar', '2'],
        ],
        chemProto = [['Models/BothPresyn77.g', 'chem']],
        chemDistrib = [
            ['glu', 'dend#', 'dend', '1', diffusionLength ],
            ['GABA', 'dend#', 'dend', '1', diffusionLength ]
        ],
        adaptorList = [
            [ 'glu/glu', 'n', 'glu', 'activation', 0.0, 0.1 ],
            [ 'GABA/GABA', 'n', 'GABA', 'activation', 0.0, 0.1 ],
        ],
        plotList = [
            #['dend#', '1', 'glu/glu', 'conc', 'glu Conc'],
            #['dend#', '1', 'GABA/GABA', 'conc', 'GABA Conc'],
            #['soma', '1', '.', 'Vm', 'Vm'],
            #['soma', '1', 'vclamp', 'current', 'ISoma', 'time', -0.01e-9, 0]
            ['soma', '1', 'vclamp', 'current', 'ISoma'],
            ['dend#', '1', 'glu/glu', 'conc', 'glu Conc'],
        ],
        stimList=[
            ['soma', '1', '.', 'vclamp', '-0.07' ],
            ['dend#', '1', 'glu/Ca_ext', 'concInit', '80e-6 + (t>=0.2 && t < 0.202)*'+str(gluCa)],
            ['dend#', '1', 'GABA/Ca_ext', 'concInit', '80e-6'],
        ]
    )
    gluReceptor = moose.element( '/library/glu' )
    #gluReceptor.tau2 *= 0.4
    gluReceptor.tau2 = tau / 1000
    rdes.buildModel()
    return rdes


def main():
    stimIdx = int( round( (0.2) /elecDt ) )
    pkIdx = int( round( (0.2 + 0.005) /elecDt ) )
    for tau in [ 10, 20, 30, 40, 50, 60, 70, 80 ]:
        rdes = buildModel( tau )    # tau in ms.
        moose.reinit()
        moose.start( 0.5 )
        plotEPSC = moose.element( '/model/graphs/plot0' )
        baseline = np.mean( plotEPSC.vector[stimIdx - 10:stimIdx] )
        y = plotEPSC.vector[pkIdx:] - baseline
        x = np.linspace( 0, len(y)*elecDt, len(y), endpoint = False )
        ret, cov = sci.curve_fit(lambda t,a,tau: a*np.exp(-t/tau), x, y, p0=(min(y),0.02) )
        print ( "tau = {:.3f}, {:.3f}, pk={:.2f}pA, stim={:.1f}uM".format( 
            tau, ret[1] * 1000, min(y)*1e12, gluCa*1000 ) )
        rdes.display()
        moose.delete( "/model" )
        moose.delete( "/library" )


    #rdes.display()


if __name__ == "__main__":
    main()


