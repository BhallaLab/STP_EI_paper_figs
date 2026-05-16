

# Single Neurons Detect Spatiotemporal Activity Transitions through STP and EI imbalance

## Code for simulations, data analysis, and figure generation.

Copyright (C) 2026 Aditya Asopa and Upinder S. Bhalla, 
National Centre for Biological Sciences,
Tata Institute of Fundamental Research, Bangalore, India.

All code in this project is licensed under GPL 3.0

This paper is submitted to eLife 2026. A preliminary version is on bioRxiv:

"Short-term plasticity of EI balance at single neurons can detect pattern transitions"

Aditya Asopa and Upinder S. Bhalla

bioRxiv 2024.10.30.621034; doi: https://doi.org/10.1101/2024.10.30.621034

Experimental data are available on Zenodo: https://doi.org/10.5281/zenodo.20179470

## About
This is the code and data repository for this study, notably including the
entire model specification, figure generation, and data analysis code.

Large experimental data files (>5 MB) are available on Zenodo rather than in
this repository. Scripts that require them expect the data files in the same
directory as the plotting script.

## Dependencies

All the code here uses **Python** version 3.10.

All the simulation code uses **MOOSE** version 4.1.0 Jhangri.

The model optimization code also depends on **HOSS** and **FindSim**.

Some of the MOOSE graphics depends on the 3-D rendering package **vpython**.

[MOOSE](https://github.com/BhallaLab/moose-core) is the Multiscale Object
Oriented Simulation environment. It is good for running ODE-based signaling
models defined in the standard SBML, as well as multiscale models combining
electrical and chemical signaling. The current project depends on its
multiscale capabilities.

[FindSim](https://github.com/BhallaLab/FindSim) is the Framework for
Integrating Neuronal Data and Signaling Models.

[HOSS](https://github.com/BhallaLab/HOSS) provides a set of methods for
performing hierarchical optimization of signaling and other models.

[Jardesigner](https://github.com/upibhalla/jardesigner) is three things: A JSON format for storing multiscale neuronal model definitions; a Python library for building models using this format, and a Javascript app for building and editing jardesigner models in a GUI. The Python library is included in the Zenodo repository with this paper. The working version of jardesigner on GITHUB is a8c8373.
All the simulation scripts here use a model defined in the Jardesigner JSON format.


## CONTENTS

This repository contains the following subdirectories:

-	`FIG1_basic_data/`		- Loads and plots basic experimental EI dynamics data for Fig 1. Also contains `datapaths.py` (used by all figure scripts) and the `eidynamics` Python package.
-	`FIG2_E_I_component_heatmaps/`	- Plots E and I synaptic component heatmaps for Fig 2.
-	`FIG3_EI_balance_expts/`	- Plots EI balance experimental data for Fig 3 and the supplementary panel for Fig 3.
-	`FIG4_Model_layout/`		- Generates dynamic 3-D images of the multiscale model in your browser.
-	`FIG5_burst_EPSP/`		- Simulates burst EPSP responses (`jgenFig5.py`) and plots Fig 5 (`plotFig5.py`).
-	`FIG6_Poisson_stim/`		- Simulates Poisson stimulus responses and plots Fig 6, including no-STP and no-NMDA control variants.
-	`FIG7_Mismatch_Heatmaps/`	- Simulates mismatch responses for panels AB, C, and E (jitter), and plots Fig 7 heatmaps.
-	`FIG8_spiking/`			- Simulates spiking mismatch responses including frequency sweep and theta rhythm variants, and plots Fig 8.
-	`FIG9_param_sensitivity/`	- Simulates and plots parameter sensitivity analysis for Fig 9.
-	`Fig10_summary_schematic/`	- Source file (SVG/ODG) for the Fig 10 summary schematic.
-	`ANIMATION/`			- Uses MOOSE 3D graphics to generate a movie of the simulation of mismatch detection.
-	`MODEL_FITTING/`		- Code for fitting the model to extracted voltage-clamp peak and timing data. Also contains code for Supplementary figures 1, 2, and 3.
-	`SUPP/`				- Code for plotting Extended figures 3 to 9. Also contains simulated data output for the Poisson train stimulus for the reference model, model without STP, and model without NMDA.
-	`SUPP_TABLES_PARAMS/`		- Supplementary methods PDF, tabulated chemical-electrical model parameters, and model equations.
-	`polygonProtocols/`		- Photostimulation protocol definition files used in the experiments.


In the directories for FIG5 to FIG9, any file of the form jgen\*.py is for performing simulations to generate the output .h5 files. Most of these, with the exception of FIG5, require to be run on a moderate server with 24 or more threads.


## DATA LOCATION AND SIZE

The various plotFig\*.py scripts are for plotting the figures. They need simulation data and in many cases experimental data. These data files are in .h5 format. The paths to the data will need to be updated in the scripts, as per your system. The required experimental data files are in the zenodo repository. You'll need to run the jgen scripts to get the simulation data files. A few of the smaller .h5 simulation files are provided in this repository. For reference, here are some of the file sizes:

Zenodo repository: ~7 GB

FIG5: ~7 MB
FIG6: ~100 MB
FIG7: ~143 MB
FIG8: ~781 MB
FIG9: ~1181 MB



