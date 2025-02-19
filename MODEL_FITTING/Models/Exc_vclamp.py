import moose
import matplotlib.pyplot as plt
import numpy as np
import rdesigneur as rd
chemDt = 1e-4
chemPlotDt = 2e-4
elecDt = 50e-6
diffusionLength = 1e-3 # ensure single chem compt

def load( scaleParam = [], chemFile = None ):
    if chemFile == None:
        chemFile = 'Models/BothPresyn79.g' #Use the default

    tau2 = TAU  # Replace tau with a cell-specific value obtained from Look up fitEphysToTau.py
    RM = '1.0'
    RA = '1.0'
    for idx in range( 0, len( scaleParam ), 3 ):
        if scaleParam[idx+1] == 'tau2':
            tau2 = scaleParam[idx+2]
        elif scaleParam[idx+1] == 'RM':
            RM = str(scaleParam[idx+2])
        elif scaleParam[idx+1] == 'RA':
            RA = str(scaleParam[idx+2])

    rdes = rd.rdesigneur(
        turnOffElec = True,
        verbose = False,
        chemDt = chemDt,
        chemPlotDt = chemPlotDt,
        elecDt = elecDt,
        elecPlotDt = elecDt,
        diffusionLength = diffusionLength,
        cellProto = [['ballAndStick', 'soma',  10e-6, 10e-6, 2e-6, 200e-6, 1]],
        passiveDistrib = [
            ['#', 'CM', '0.01', 'Em', '-0.065', 'RM', RM, 'RA', RA],
        ],
        chanProto = [
            ['make_glu()', 'gluR'],
            ['make_GABA()', 'GABAR'],
        ],
        chanDistrib = [
            ['gluR', 'dend#', 'Gbar', '50'],
            ['GABAR', 'dend#', 'Gbar', '0.00000000001'],
        ],
        chemProto = [[chemFile, 'chem']],
        chemDistrib = [
            ['glu', 'dend#', 'presyn_dend', '1', 0.5e-6, 0, 200e-6],
            ['GABA', 'dend#', 'presyn_dend', '1', 0.5e-6, 0, 200e-6]
        ],
        adaptorList = [
            [ 'glu/glu', 'n', 'gluR', 'activation', 0.0, 50 ],
            [ 'GABA/GABA', 'n', 'GABAR', 'activation', 0.0, 0.1e-9 ],
        ],
    )
    gluReceptor = moose.element( '/library/gluR' )
    gluReceptor.tau2 = tau2
    return rdes

def build( rdes ):
    rdes.buildModel()

def main():
    rdes = load()
    build(rdes)

if __name__ == "__main__":
    main()
