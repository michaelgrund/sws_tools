function [SPLITS,NULLS]=SWS_read_evstruct(varargin)
%
% This function reads structs with data of the full data set published by
% Grund & Ritter (2020). Functionality with other data sets is not
% guaranteed. For the required mat-file structure see the data available
% from KITopenData via: 
%
%  https://publikationen.bibliothek.kit.edu/1000091427
%
% How to use:
%
% 1) add the directory < 02_data_handling > to the matlab path using
%    addpath('.../SWStools/02_data_handling') or the Set Path dialog box
%
% 2) switch to directory that contains the content of the downloaded
%    and unzipped splitting data set: <<< mgrund_diss_2019_ELAPP >>>
%
% 3) run this m-file in the command window to generate different outputs
%    (depending on your choice): 
%
%    (a) [SPLITS,NULLS]=SWS_read_evstruct => gives only the full structs 
%        (e.g. for further analysis)
%    (b) SWS_read_evstruct('stereo') => allows to directly generate (single) 
%        stereoplots for the analyzed stations via further user input
%    (c) SWS_read_evstruct('stereoall') => automatically generate
%        stereoplots for all 266 seismic stations
%    (d) SWS_read_evstruct('histo') => generates a histogram showing the 
%        overall distribution for fast axis (phi) and delay time (dt),
%        similar to Fig. 7.6 in the dissertation
%
% 2019-04-10 -MG- (https://orcid.org/0000-0001-8759-2018)
%
% see also functions: SWS_stereoplot_col, SWS_histogram
%===============================================================================

close all 
clc

% check for directory
currdir=pwd;
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


%==================================================================================
% STEREOLOTS individual
if ~isempty(varargin) && strcmp(varargin,'stereo')
    
    % make new directory to store generated stereoplots of the current
    % session
    cd(currdir)
    savedir=['PLOTTING_stereoplots_' datestr(now,'dd-mmm-yyyy_HH_MM_SS')];
    mkdir(savedir)

    %............................
    %subfunction
    func_makeplot(NULLS,SPLITS,savedir)
    %............................
end
%==================================================================================

%==================================================================================
% STEREOLOTS all
if ~isempty(varargin) && strcmp(varargin,'stereoall')
    
    % make new directory to store generated stereoplots of the current
    % session
    cd(currdir)
    savedir=['PLOTTING_stereoplots_' datestr(now,'dd-mmm-yyyy_HH_MM_SS')];
    mkdir(savedir)   

    disp(' ') 
    disp('Generate stereoplots of all 266 seismic stations...');
    pause(2)
    
    %............................
    %subfunction
    func_makeplot_all(NULLS,SPLITS,savedir)
    %............................
    
    disp(' ') 
    disp('Done!');

end
%==================================================================================

%==================================================================================
% HISTOGRAM

if ~isempty(varargin) && strcmp(varargin,'histo')

    cd(currdir)
    
    disp(' ') 
    disp('Generate histogram for given data...');
    pause(2)

    phisall=[SPLITS.phiSC];
    phisall=phisall(2:3:end);

    dtsall=[SPLITS.dtSC];
    dtsall=dtsall(2:3:end);

    %............................
    SWS_histogram(phisall,dtsall)
    %............................

end
%==================================================================================

end %EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

function func_makeplot(NULLS,SPLITS,savedir)

% Type in station name
disp(' ')
findname=input('Make stereoplot. Please insert station name (e.g. PVF):','s');

indxNULL= strcmp({NULLS.staname},findname);
indxSP= strcmp({SPLITS.staname},findname);

NULLSsel=NULLS(indxNULL);
SPLITSsel=SPLITS(indxSP);

%================================================================

if ~isempty(SPLITSsel) || ~isempty(NULLSsel)
    
    % Want to plot bars in color?
    disp(' ')
    fast_col=input('Color-code bars with respect to fast axis  ([0]=no, [1]=yes):');
    disp(' ')
        
    %............................
    SWS_stereoplot_col(fast_col,SPLITSsel,NULLSsel,savedir)
    %............................

else
    disp(' ')
    disp('Station name not available! Try another one!')
    disp(' ')
end

%================================================================
disp(' ')
plotother=input('Plot other station?  ([0]=no, [1]=yes):');
disp(' ')

if plotother == 1
    func_makeplot(NULLS,SPLITS,savedir)
else
    return
end
%================================================================

end %EOsubF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_makeplot_all(NULLS,SPLITS,savedir)

% Want to plot bars in color? 
disp(' ') 
fast_col=input('Color-code bars with respect to fast axis  ([0]=no, [1]=yes):');
disp(' ')

staunq=unique(horzcat({NULLS.staname},{SPLITS.staname}));

for ii=1:length(staunq)
    
    findname=staunq(ii);
    indxNULL= strcmp({NULLS.staname},findname);
    indxSP= strcmp({SPLITS.staname},findname);

    NULLSsel=NULLS(indxNULL);
    SPLITSsel=SPLITS(indxSP);
    
    %............................
    SWS_stereoplot_col(fast_col,SPLITSsel,NULLSsel,savedir)
    %............................

end

%================================================================

end %EOsubF
