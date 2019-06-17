# SWStools
MATLAB tools for analyzing and plotting of shear wave splitting data

If you make use of the content in this repository please acknowledge my PhD thesis in whose framework several of the scripts and functions were developed:

- **_Grund, M. (2019)_**, Exploring geodynamics at different depths with shear wave splitting, Dissertation, Karlsruhe Institute of Technology (KIT), https://doi.org/10.5445/IR/1000091425 

Content of the directories:

[01_stereoplots](https://github.com/michaelgrund/SWStools/tree/master/01_stereoplots)
Generate stereoplots of splitting results based on output of **SplitLab** and **StackSplit**.

- color-code bars of splits with respect to backazimuth (see below for supported colormaps)
- shade full background or specific backazimuthal wedges
- produce publication-ready plots in different formats (pdf, png, jpg...)

Supported colormaps 
  1) standard MATLAB colormaps: parula, winter, summer, copper....

  2) MatPlotLib 2.0 Colormaps: Perceptually Uniform and Beautiful 
    (see and download here: https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-2-0-colormaps-perceptually-uniform-and-beautiful)

  3) Scientific colour-maps by F. Crameri (Zenodo. http://doi.org/10.5281/zenodo.1243862)
    (see and download here: http://www.fabiocrameri.ch/colourmaps.php)
    
Run function `SWS_Analysis_BASICS_stereoplots` in the MATLAB command window to generate the stereoplots shown in the figure below using the provided test data set of station VAF.    


![PLOT_res_ALL](https://user-images.githubusercontent.com/23025878/56903070-dfe03a00-6a9b-11e9-9cc0-606d9c2a4173.png)

##[02_data_handling](https://github.com/michaelgrund/SWStools/tree/master/02_data_handling)

Read and visualize the shear wave splitting data available from: https://doi.org/10.5445/IR/1000091427

- generate individual stereoplots `SWS_read_evstruct('stereo')` 
- generate stereoplots of all 266 analyzed stations `SWS_read_evstruct('stereoall')`
- generate histogram `SWS_read_evstruct('histo')`
- generate output structs separated in split and null measurements `[SPLITS,NULLS]=SWS_read_evstruct`

![PLOTS_merged](https://user-images.githubusercontent.com/23025878/57140316-c864c200-6db7-11e9-933b-b96c452f8349.png)
