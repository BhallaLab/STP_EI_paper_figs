#########################################################################
## This program is part of 'MOOSE', the
## Messaging Object Oriented Simulation Environment.
##           Copyright (C) 2014 Upinder S. Bhalla. and NCBS
## It is made available under the terms of the
## GNU Lesser General Public License version 2.1
## See the file COPYING.LIB for the full notice.
#########################################################################

import moose
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab
import numpy
import sys
import os

scriptDir = os.path.dirname( os.path.realpath( __file__ ) )

preTime = 0.1
stimTime = 0.002
postTime = 0.040

def main():
    """
This example illustrates loading, and running a kinetic model
for synaptic release.
    """

    solver = "gsl"  # Pick any of gsl, gssa, ee..
    #solver = "gssa"  # Pick any of gsl, gssa, ee..
    runtime = 2000.0
    fname = "Models/ExcPresyn76.g"
    if ( len( sys.argv ) == 2 ):
        fname = sys.argv[1]
    modelId = moose.loadModel( fname, 'model', solver )
    plot1 = moose.Table( '/model/graphs/Ca' )
    moose.connect( plot1, 'requestOut', '/model/kinetics/glu/Ca', 'getN')
    plot2 = moose.Table( '/model/graphs/glu' )
    moose.connect( plot2, 'requestOut', '/model/kinetics/glu/glu','getN')
    # Increase volume so that the stochastic solver gssa
    # gives an interesting output
    compt = moose.element( '/model/kinetics' )
    compt.volume = 5e-19
    for ii in range( 8, 20):
        moose.setClock( ii, 1e-4 ) # for the plots

    moose.reinit()
    moose.element( '/model/kinetics/glu/Ca_ext' ).concInit = 80e-6
    moose.start( preTime )
    moose.element( '/model/kinetics/glu/Ca_ext' ).concInit = 50e-3
    moose.start( stimTime )
    moose.element( '/model/kinetics/glu/Ca_ext' ).concInit = 80e-6
    moose.start( postTime )

    # Display all plots.
    startIdx = int(numpy.round( preTime * 0.9 / plot1.dt ))
    t = numpy.array(numpy.arange( 0, plot1.vector.size, 1 ) * plot1.dt)
    t = t[startIdx:]
    fig = plt.figure( figsize=(12, 10 ) )
    ax = fig.add_subplot( 311 )
    ax.plot( t, plot1.vector[startIdx:], 'b-', label=plot1.name )
    ax = fig.add_subplot( 312 )
    y = plot2.vector[startIdx:]
    ax.plot( t, y, 'c-', label=plot2.name )
    ax = fig.add_subplot( 313 )
    tot = 0.0
    y2 = []
    for vv in y:
        tot += vv
        y2.append( tot )
    ax.plot( t, y2 , 'r-', label=plot2.name )
    plt.ylabel( '# molecules' )
    plt.xlabel( 'Time (s)' )
    pylab.legend()
    pylab.show()

# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
        main()
