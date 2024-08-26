# 1D wave model

This repository is the simple 1D wave model written in MATLAB to model the progression of wave heights from offshore towards the land. The model presented in the manuscript "Insert title". 

The model captures four different mechanisms that affect wave propagation in the cross-shore direction, including wave breaking, vegetation drag, shoaling, and bed friction. In particular, the vegetation drag is computed based on plant morphology without any emprirically calibrated parameters. At each spatial step, the model takes wave height and the effects of the four mechanisms to compute the wave height at the next spatial step. 

## File access
- The 'func' folder contains the MATLAB functions to run the 1-D wave model example script 'wave_model_example.m'.
- The 'script' folder contains an example of how to apply the 1-D wave model in the script 'wave_model_example.m'. The inputs required are 1-D bathymetry profile at a resolution that matches with the step interval between each wave height modeling, vegetation species and their morphology, and offshore wave period and wave height. 
