# 1D wave model

This repository is the simple 1D wave model written in MATLAB to model the progression of wave heights from offshore to shore. The model presented in the manuscript "Insert title". 

The model captures four different mechanisms that affect wave propagation in the cross-shore direction, including wave breaking, vegetation drag, shoaling, and bed friction. In particular, the vegetation drag is computed based on plant morphology without any emprirically calibrated parameters. At each spatial step, the model takes wave height and the effects of the four mechanisms to compute the wave height at the next spatial step. 

# File access
- The 'srcs' folder includes the MATLAB functions to construct the 1D wave model.
- The 'script' folder includes an example of how to apply the 1D wave model. The inputs required are 1D bathymetry profile at a resolution that matches with the step interval between each wave height modeling, vegetation species and their morphology, offshore wave period and wave height. 
