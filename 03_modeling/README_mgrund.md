# Instructions how to perform a modeling run for a single seismic station

- Created: M. Grund 2020/06/02
- Completed: Y. Fröhlich 2023/07/31
- Details: [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of [Grund & Ritter (2020)](https://doi.org/10.1093/gji/ggaa388)

Official MATLAB Toolboxes required for running the modeling routines:

- Deep Learning Toolbox
- Mapping Toolbox

--------------------------------------------------------------------------------------

The _MATLAB Seismic Anisotropy Toolkit_ (MSAT) is required:

1) Download MSAT from https://www1.gly.bris.ac.uk/MSAT/ or https://github.com/andreww/MSAT.

2) Add the whole MSAT package to your MATLAB path.

3) Optional: If you want to get familiar with the modeling in MSAT and behavior
   of different settings change to the directory `MSAT/examples/splitting_model/`
   in which the script `split_model.m` is located (this script is described in the
   paper of [Walker and Wookey (2012))](https://doi.org/10.1016/j.cageo.2012.05.031).
   There you can play around with different parameters and see how these affect the
   splitting parameters.

--------------------------------------------------------------------------------------

First we pre-compute the different models. Required disc space for exemplary settings:

- one-layer models: function `SWS_modeling_precomp_single_layer.m`
  - phi = [-90:5:90] deg and dt = [0.2:0.2:4] s
  - disc space: 7 MB
- two-layer models: function `SWS_modeling_precomp_twolayers.m`
  - phi = [-90:5:90] deg and dt = [0.2:0.2:4] s in each layer
  - disc space: 5.1 GB
- dipping layer models: functions `SWS_modeling_precomp_dippinglayer.m` and `SWS_modeling_calc_dipping.m`
  - downdip direction = [0:5:360] deg, dip angle = [5:5:75] deg, layer thickness = [25:25:500] km
  - disc space: 500 MB

&rarr; Total disc space: 5.7 GB

This step is time consuming since all possible variations are computed.

Then we merge all these models in a single MATLAB structure (function `SWS_modeling_precomp_models_main.m`) with fields:
1) `phi_eff`: effective phi values over backazimuth
2) `dt_eff`: effective dt values over backazimuth
3) `mod_paras`: model parameters depending on model type
4) `type`: string corresponding to the model type ('single_layer', 'two_layers', 'dipping')

--------------------------------------------------------------------------------------

Then the measured data is used to compare it against all pre-computed models and a RMSE
is calculated. Data input is expected in _SplitLab_ and _StackSplit_ output formats.

Finally, we take the minimum RMSE or something else metric to get the model which best
describes the data.
