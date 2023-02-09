function modsall_sort = SWS_modeling_calc_misfit_custom(modelsin, modrange_low, modrange_upp)
%
% fit customly prepared SWS parameters to synthetic models
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% below all data is loaded from structs which have to be 
% prepared using function SWS_modeling_prep_custom_data
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% INPUT:
% modelsin: string with the filename which consists all pre-computed models
%           (e.g. 'sws_modout_domper8s.mat'), to create models first see 
%           function SWS_modeling_precomp_models_main 
% modrange_low: lower BAZ limit to model (e.g. 10)
% modrange_upp: upper BAZ limit to model (e.g. 270)
%
% .........................................................................
% .........................................................................
% EXAMPLE 
%
% 1) change to a working directory of your choice
%
% 2) prepare your data using SWS_modeling_prep_custom_data
%
% 3) define input variables
% 
%    modelsin = 'sws_modout_domper8s.mat' % (be sure to have that file in the
%                                         % current directory!)
%    modrange_low = 3
%    modrange_upp = 90
%
% 4) run misfit routine
%    
%    modsall_sort = SWS_modeling_calc_misfit_custom(modelsin, modrange_low, modrange_upp)
%
% .........................................................................
% .........................................................................
%
% be sure to e.g. exclude discrepant pairs (SKS-SKKS) from your dataset 
% before running this function 
%
% feel free to modify/adjust the code for your needs or submit improvements
%
% bugs etc. can be reported by opening a "New issue" in the GitHub
% repository
%
% LICENSE
%
% Copyright (C) 2020  Michael Grund, Karlsruhe Institute of Technology (KIT), 
% ORCID: https://orcid.org/0000-0001-8759-2018
% GitHub: https://github.com/michaelgrund/sws_tools
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% TERMS OF USE
%
% The modeling routines are provided "as is" and without any warranty. 
% The author cannot be held responsible for anything that happens to you 
% or your equipment. Use it at your own risk.
%
%============================================================== 
%============================================================== 

clc
close all

% model_out: preprocessed models based on parameters
disp(['Loading model file <' modelsin '>...'])
models = load(modelsin);
model_out = models.splitmods;

%================================================================
% INITIAL SETTINGS
% mainly for later plotting etc.
% adjust for your needs

% plot maximum XX best models in final figure
plot_mod_max=20;
keep_mods = 500; % finally keep only the XX bests models 

% only model a specific range of the data (baz1:baz2)
modrange_col = [219,219,219]./256;
modrange_edcol = modrange_col;

% BAZ range to plot
BAZ=0:1:360;

% MODELS plotting
colmod_bf_1=[0.6350 0.0780 0.1840];  % color best model 
colmod_bf_2max=[175 175 175]./256; % color others if plot_mod_max > 1      

% SPLITS plotting
fs=8; 
ms=7;
fs_RMSE=7;
lw_symb=1; % width of edge of symbols
lw_mod=1.2;

% NULLS plotting
ms_null=6.5;
lw_symb_null=lw_symb; % width of edge of symbols
col_edge_null='k';
col_face_null='w';

%================================================================
% error checks

if ~isfield(model_out,'phi_eff') && ~isfield(model_out,'dt_eff') 
    error('Required fields <phi_eff> & <dt_eff> do not exist in variable <model_out> (synthetic models)! Check input struct!')
end

%================================================================
% load custom data 

dir_res_split=dir('RES_split_cust.mat');
dir_res_nulls=dir('RES_nulls_cust.mat');
dir_res_stack=dir('RES_stack_cust.mat');

if ~isempty(dir_res_split)
    Splits = load(dir_res_split.name);
    RES_split = Splits.RES_split;
    staname_split=RES_split(1).staname;
else
    RES_split = [];
end
    
if ~isempty(dir_res_nulls)
    Nulls = load(dir_res_nulls.name);
    RES_nulls = Nulls.RES_nulls;
    staname_split=RES_nulls(1).staname;
else
    RES_nulls = [];
end


if ~isempty(dir_res_stack)
    Stacks = load(dir_res_stack.name);
    RES_stack =Stacks.RES_stack;
    staname_split=RES_stack(1).staname;
else
    RES_stack = [];   
end

%================================================================
% fitting method:

disp(' ')
whichfit=input('Fitting method: [1] only phi (RMSE), [2] joint phi/dt (RMSE)');

if isempty(whichfit)
    whichfit = 2;
end

disp(' ')

%================================================================
% Ask for phase results that should be plotted

phaselist_split=RES_split.phase;
phaselist_null=RES_nulls.phase;
phaselist_all=unique(horzcat(phaselist_split,phaselist_null));

disp(' ')
disp('Available phase results: ')
disp(' ')

for ii=1:length(phaselist_all)+1
    if ii==1
      disp('   [0] ALL') 
    else  
        disp(['   [' num2str(ii-1) '] only ' phaselist_all{ii-1}])
    end
end

disp(' ')
wphases=input('Which phases you wanna use (Enter number [XX] from the given list)?:');
disp(' ')

if isempty(wphases)
    wphases=0;
elseif ~isempty(wphases) && ~ismember(wphases,[0:1:length(phaselist_all)])
    error('You selected a non-available entry ;) Try again!')
end

if wphases~=0
    phases2use=phaselist_all(wphases);
    find_index_split= strcmp(phaselist_split,phases2use);
    find_index_null= strcmp(phaselist_null,phases2use);
    
    RES_split=RES_split(find_index_split);
    RES_nulls=RES_nulls(find_index_null);
end

%================================================================
% prepare variables

% splits
meas_phiSC=[RES_split.phiSC_err_min; RES_split.phiSC; ...
    RES_split.phiSC_err_max]'; % create vector with [lower err phi; phi; upper err phi]
meas_phiSC(:,4)=1; % 1 in 4th column means single splits for later coloring
meas_dtSC=[RES_split.dtSC_err_min; RES_split.dtSC; ...
    RES_split.dtSC_err_max]'; % create vector with [lower err dt; dt; upper err dt]
meas_BAZ=[RES_split.baz]';

% round corresponding BAZs of measured SWS parameters (no big difference), 
% since theoretical values are only available as integers, otherwise no
% comparison possible!
meas_BAZ_floor=floor(meas_BAZ); 
               
% same for the nulls                                
meas_phiSC_null=[RES_nulls.phiSC_err_min; RES_nulls.phiSC; ...
    RES_nulls.phiSC_err_max]'; 
meas_dtSC_null=[RES_nulls.dtSC_err_min; RES_nulls.dtSC; ...
    RES_nulls.dtSC_err_max]';     
meas_BAZ_null=[RES_nulls.baz]';
meas_BAZ_floor_null=floor(meas_BAZ_null);  

% same for the stacks 
if ~isempty(RES_stack)
    meas_phiSC_stack=[RES_stack.phiSTACK_err_min;...
        RES_stack.phiSTACK; RES_stack.phiSTACK_err_max]';
    meas_phiSC_stack(:,4)=2; % 2 in 4th column means stacked splits for later coloring
    meas_dtSC_stack=[RES_stack.dtSTACK_err_min; ...
        RES_stack.dtSTACK; RES_stack.dtSTACK_err_max]';     
    meas_BAZ_stack=[RES_stack.meanbaz]';
    meas_BAZ_floor_stack=floor(meas_BAZ_stack);  
    
    % merge single splits with stacked splits   
    merged_phiSC=vertcat(meas_phiSC,meas_phiSC_stack);
    merged_dtSC=vertcat(meas_dtSC,meas_dtSC_stack);
    merged_BAZ_floor=vertcat(meas_BAZ_floor,meas_BAZ_floor_stack);
  
    % rename original variables
    meas_phiSC=merged_phiSC;
    meas_dtSC=merged_dtSC;
    meas_BAZ_floor=merged_BAZ_floor;

    % for plotting only
    meas_phiSC4plot=merged_phiSC;
    meas_dtSC4plot=merged_dtSC;
    meas_BAZ_floor4plot=merged_BAZ_floor;
 
else
    meas_phiSC4plot=meas_phiSC;
    meas_dtSC4plot=meas_dtSC;
    meas_BAZ_floor4plot=meas_BAZ_floor;
end
%================================================================
% sort to model range

findvals=find(meas_BAZ_floor > modrange_low & meas_BAZ_floor < modrange_upp);

meas_BAZ_floor=meas_BAZ_floor(findvals);
meas_phiSC=meas_phiSC(findvals,:);
meas_dtSC=meas_dtSC(findvals,:);

%================================================================
% calc deviations of measured SWS parameters from theoretical distributions
% for each synthetic model, then compute RMSE

count_mods=1;

res_dt=zeros(length(meas_phiSC),1);
res_phi=res_dt;

for ii=1:length(model_out)

    curr_mod_phi=model_out(ii).phi_eff;
    curr_mod_dt=model_out(ii).dt_eff;
    curr_mod_type=model_out(ii).type;

    sizeSC = size(meas_phiSC);
    
    for kk=1:sizeSC(1)

        % find theoretical value for BAZ of corresponding measured value
        find_theo_phi=curr_mod_phi(meas_BAZ_floor(kk));
        find_theo_dt=curr_mod_dt(meas_BAZ_floor(kk));
       
        %calc residuum
        res_phi(kk,1)=abs(find_theo_phi-meas_phiSC(kk,2));
        
        % if there is a residuum > 90, calc 180-res_phi to consider 
        % the flip from -90° to +90°
        if res_phi(kk,1) > 90  
            res_phi(kk,1)=180-res_phi(kk,1);
        end
        
        res_dt(kk,1)=abs(find_theo_dt-meas_dtSC(kk,2));
    end

         if strcmp(curr_mod_type,'two_layers')
             
            modsall(count_mods).phi_eff=curr_mod_phi;
            modsall(count_mods).dt_eff=curr_mod_dt;
            modsall(count_mods).mod_type=curr_mod_type;
            modsall(count_mods).phi=model_out(ii).mod_paras.phi_in;
            modsall(count_mods).dt=model_out(ii).mod_paras.dt_in;      
            modsall(count_mods).modrange_low=modrange_low; 
            modsall(count_mods).modrange_upp=modrange_upp; 
            modsall(count_mods).RMSE_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
            modsall(count_mods).RMSE_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));
            modsall(count_mods).staname=staname_split;
            modsall(count_mods).azi4plot=[];
            modsall(count_mods).fast4plot=[];
            modsall(count_mods).dt4plot=[];
          
         elseif  strcmp(curr_mod_type,'single_layer')

            modsall(count_mods).phi_eff=curr_mod_phi;
            modsall(count_mods).dt_eff=curr_mod_dt;
            modsall(count_mods).mod_type=curr_mod_type;
            modsall(count_mods).phi=model_out(ii).mod_paras.phi_in;
            modsall(count_mods).dt=model_out(ii).mod_paras.dt_in;     
            modsall(count_mods).modrange_low=modrange_low; 
            modsall(count_mods).modrange_upp=modrange_upp; 
            modsall(count_mods).RMSE_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
            modsall(count_mods).RMSE_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));
            modsall(count_mods).staname=staname_split;
            modsall(count_mods).azi4plot=[];
            modsall(count_mods).fast4plot=[];
            modsall(count_mods).dt4plot=[];
  
         elseif strcmp(curr_mod_type,'dipping')
    
            modsall(count_mods).phi_eff=curr_mod_phi;
            modsall(count_mods).dt_eff=curr_mod_dt;
            modsall(count_mods).mod_type=curr_mod_type;
            modsall(count_mods).downdipdir=model_out(ii).mod_paras.downdipdir;
            modsall(count_mods).dip=model_out(ii).mod_paras.dip;
            modsall(count_mods).thick=model_out(ii).mod_paras.thick; 
            modsall(count_mods).modrange_low=modrange_low; 
            modsall(count_mods).modrange_upp=modrange_upp; 
            modsall(count_mods).RMSE_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
            modsall(count_mods).RMSE_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));
            modsall(count_mods).staname=staname_split;
            modsall(count_mods).azi4plot=model_out(ii).mod_paras.azi4plot;
            modsall(count_mods).fast4plot=model_out(ii).mod_paras.fast4plot;
            modsall(count_mods).dt4plot=model_out(ii).mod_paras.dt4plot;
          
         end

          if whichfit==2
          % use phi and dt for joint fitting
            modsall(ii).RMSE= modsall(ii).RMSE_phi/90+modsall(ii).RMSE_dt/4;  
            nameend='joint';
          
          elseif whichfit==1

          % only use phi for fitting
            modsall(ii).RMSE=modsall(ii).RMSE_phi/90;
            nameend='only_phi';

          end
          
    count_mods=count_mods+1;
end

%================================================================
% sort models based on total RMSE
[b,index]=sort([modsall.RMSE]); 
modsall_sort=modsall(index); % entry 1 corresponds to minimum overall RMSE

% keep only XXX best models
modsall_sort=modsall_sort(1:keep_mods);

%================================================================
% plotting results 

% single splits
colorsedge(1,:)=[0 0 0];
colorsfill(1,:)=[66 91 169]./256;

% multievent
colorsedge(2,:)=[0 0 0];
colorsfill(2,:)=[79 197 104]./256;

%###########################################

% model fit and model distribution
SWS_modeling_plot_results(BAZ,modsall_sort,plot_mod_max,...
    meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,colmod_bf_1,colmod_bf_2max,lw_mod,...
    colorsfill,colorsedge,ms,ms_null,lw_symb,...
    lw_symb_null,col_face_null,col_edge_null,fs,...
    fs_RMSE,modrange_col,modrange_edcol,meas_BAZ_floor4plot,...
    meas_phiSC4plot,meas_dtSC4plot,staname_split,nameend)

% stereoplot displaying the splitting parameter distribution of the
% best models (based on lowest RMSE, see above)
plotnum = 1; % e.g. 1:5 to plot stereoplot for 5 best models 

for ii=plotnum
    SWS_modeling_plot_stereo_synthetic(modsall_sort,ii)
end

%###########################################

% EOF  
