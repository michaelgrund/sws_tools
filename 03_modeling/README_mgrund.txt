sws_tools

M. Grund 2020/06/02
======================================================================================
======================================================================================
Instructions how to perform a modeling run using the modified MSAT functions
for a single seismic station.

MATLAB Toolboxes required for running the modeling and other plotting stuff in sws_tools:
- Deep Learning Toolbox
- Mapping Toolbox


Data input is expected in SplitLab format. Other formats can be adjusted....

1) First download MSAT from XXXXXX

2) Please place the whole content of the directory XXXX in MSATS example directory.

3) Add the whole MSAT package to your MATLAB path.

4) change to directory XXX and start modeling by running function XXXXX

5) optional: if you want to get familiar with the modeling in MSAT and behaviour of different settings
   change to the directory XXX in which the script split_model of is placed (this script is described 
   in the paper of Walker and Wookey (2012)). There you can play around
   with different parameters and see what's going on.....




======================================================================================

First we precompute the different models (reuqired disc space for exemplary settings are given below):

- two layer models : function XXXX (disc space for phis between -90 to 90 and 5° steps and dts between 0 and 4 s, steps 0.2 s for each layer: XXXX GB)
- one and dipping layer models: function XXXX (disc space for XXX, XXX and XXXX: XXXX GB)

This step is time consuming since all possible variations are computed. 

======================================================================================

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

======================================================================================
======================================================================================
visualize content/generate stereoplotls of GR2020...


======================================================================================
======================================================================================
make stereoplots with your own data....