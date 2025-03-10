# Single Neurons Detect Spatiotemporal Activity Transitions through STP and EI imbalance

## Code for generating supplementary movie for simulation.

Copyright (C) 2025 Aditya Asopa and Upinder S. Bhalla, 
National Centre for Biological Sciences,
Tata Institute of Fundamental Research, Bangalore, India.

All code in this project is licensed under GPL 3.0

There is a preliminary version of this paper out on bioRxiv:

"Short-term plasticity of EI balance at single neurons can detect pattern transitions"

Aditya Asopa and Upinder S. Bhalla

bioRxiv 2024.10.30.621034; doi: https://doi.org/10.1101/2024.10.30.621034

## About
This is the code and data repository for the supplementary movie to show
mismatch detection using STP and EI balance.

## Dependencies

All the code here uses **Python** version 3.10.

All the simulation code uses **MOOSE** version 4.1.0 Jhangri.

The MOOSE graphics in this figure depends on the 3-D rendering package **vpython**.

The conversion from a series of pngs to the movie is performed by **ffmpeg**:

	ffmpeg version 4.4.2-0ubuntu0.22.04.1

[MOOSE](https://github.com/BhallaLab/moose-core) is the Multiscale Object 
Oriented
Simulation environment. It is good for running ODE-based signaling models
defined in the standard SBML, as well as multiscale models combining electrical
and chemical signaling. The current project depends on its multiscale 
capabilities.

## CONTENTS

This directory contains the code for the movie. Note that this is a real-time
display of a running simulation. Each frame is stored as a .png and then
the frames are put together into a movie.

## GENERATING THE MOVIE

`python movie_mismatch_simulation.py`

To convert the pngs into the mp4: 

`ffmpeg -framerate 30 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p mismatch_Asopa_and_Bhalla2025.mp4`



