function [SPLITS,NULLS]=SWS_modeling_read_evstruct(varargin)
%
% This function reads structs with data of the full data set published by
% Grund & Ritter (2020). Functionality with other data sets is not
% guranteed. For the required mat-file structure see the data available
% from KITopenData via: 
%
%  https://publikationen.bibliothek.kit.edu/1000091427
%
% How to use:
%
% 1) switch to directory that contains the content of the downloaded
%    and uzipped splitting data set: <<< mgrund_diss_2019_ELAPP >>>
%
% 2) run this m-file in the command window to generate outputs:
%
%    [SPLITS,NULLS]=SWS_modeling_read_evstruct 
%
% 2020-06-04 -MG- (michael.grund@kit.edu)
%
%===============================================================================

close all 
clc

% check for directory
dirfold=dir('*_2019_ELAPP');

if length(dirfold) > 1
   warning('More than one folder in current directory!')   
end

if ~isempty(dirfold) && isfolder(dirfold.name)
    cd(dirfold.name)
end

% check for structs downloaded from KITopenData 
cd('02_SWS_splits')
dirS=dir('03_SWS_SA_splits_full.mat');
cd ..
cd('03_SWS_nulls')
dirN=dir('03_SWS_SA_nulls_full.mat');
cd ..

if ~isempty(dirN) && ~isempty(dirS)

    disp(' ')
    disp('Load structs...')
    
    cd('03_SWS_nulls')
    load 03_SWS_SA_nulls_full.mat 
    NULLS=dataNULLS;
    cd ..
  
    cd('02_SWS_splits')
    load 03_SWS_SA_splits_full.mat 
    SPLITS=dataSPLITS;

    clear dataNULLS dataSPLITS;

else
    disp(' ')
    disp('Files not found, please check directory!')
    disp(' ')
end