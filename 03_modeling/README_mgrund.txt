sws_tools/03_modeling

Instructions how to perform a modeling run for a single seismic station

M. Grund 2020/06/02

offical MATLAB Toolboxes required for running the modeling routines:

- Deep Learning Toolbox
- Mapping Toolbox

======================================================================================
- MSAT package required! First download from https://www1.gly.bris.ac.uk/MSAT/

1) Install MSAT.

2) Add the whole MSAT package to your MATLAB path.

3) optional: if you want to get familiar with the modeling in MSAT and behaviour 
   of different settings change to the directory MSAT/examples/splitting_model/ 
   in which the script split_model.m is located (this script is described 
   in the paper of Walker and Wookey (2012)). There you can play around
   with different parameters and see how these affect the splitting parameters.

======================================================================================

First we precompute the different models (required disc space for exemplary settings 
are given below):

- two layer models : function XXXX (disc space for phis between -90 to 90 and 5° steps and dts between 0 and 4 s, steps 0.2 s for each layer: XXXX GB)
- one and dipping layer models: function XXXX (disc space for XXX, XXX and XXXX: XXXX GB)

This step is time consuming since all possible variations are computed. 

======================================================================================

Data input is expected in SplitLab and/or StackSplit output format. 


Then we merge all these models in a single file (function XXXX) with columns:
1) phi_eff: effective phi values over backazimuth
2) dt_eff: effective dt values over backazimuth
3) two-layer parameters (if current model is two-layer)
4) dipping layer parameters (if current model is dipping layer)

======================================================================================

Then the measured data is used to compare it against all pre computed models and
a RMSE is calculated

======================================================================================

Finally, we take the minimum RMSE or something else metric to get the model which best 
describes the data.
