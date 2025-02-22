

# Single Neurons Detect Spatiotemporal Activity Transitions through STP and EI imbalance

## Code for simulations, data analysis, and figure generation.

Copyright (C) 2025 Aditya Asopa and Upinder S. Bhalla, 
National Centre for Biological Sciences,
Tata Institute of Fundamental Research, Bangalore, India.

All code in this project is licensed under GPL 3.0

There is a preliminary version of this paper out on bioRxiv:

"Short-term plasticity of EI balance at single neurons can detect pattern transitions"

Aditya Asopa and Upinder S. Bhalla

bioRxiv 2024.10.30.621034; doi: https://doi.org/10.1101/2024.10.30.621034

## About
This is the code and data repository for this study, notably including the
entire model specification, figure generation, and data analysis code.


## Dependencies

All the code here uses **Python** version 3.10.
All the simulation code uses **MOOSE** version 4.1.0 Jhangri
The model optimization code also depends on  **HOSS** and **FindSim** 
Some of the MOOSE graphics depends on the 3-D rendering package **vpython**

[MOOSE](https://github.com/BhallaLab/moose-core) is the Multiscale Object 
Oriented
Simulation environment. It is good for running ODE-based signaling models
defined in the standard SBML, as well as multiscale models combining electrical
and chemical signaling. The current project depends on its multiscale 
capabilities.

[FindSim](https://github.com/BhallaLab/FindSim) is the Framework for Integrating
Neuronal Data and Signaling Models. 

[HOSS](https://github.com/BhallaLab/HOSS) HOSS provides a set of methods for performing hierarchical optimization of signaling and other models.


## CONTENTS

This repository contains the following subdirectories:

-	`FIG3_Presyn_model/`	- Analyzes hdf5-formatted voltage clamp data, converts to a Pandas-hdf5 file with extracted peaks and timings.
-	`FIG4_Model_layout/`	- Generates dynamic 3-D images of the multiscale model in your browser
-	`FIG5_burst/`		- Simulates burst data and plots Fig 5.
-	`FIG6_Poisson_train/`	- Simulates Poisson train data and plots Fig 6.
-	`FIG7_mismatch/`	- Simulates mismatch responses and plots Fig 7
-	`FIG8_spiking_mismatch/`- Simulates spiking mismatch and plots Fig 8
-	`ANIMATION/`		- Uses MOOSE 3D graphics to generate a movie of the simulation of mismatch detection.
-	`MODEL_FITTING/`	- Code for fitting model to extracted peak and timing data. Also code for Supplementary figures 1, 2 and 3 which use this.
-	`SUPP/`			- Code for plotting Extended figures 3 to 9. Also simulated data output for Poisson train stimulus, for reference model, model without STP, and model without NMDA.


