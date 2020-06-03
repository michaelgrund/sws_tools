function GR_SWS_modeling_calc_misfit(modelsin)
% modelsin: string with the filename which consists all pre-computed models
%
% function to fit measured SWS parameters to synthetic models
% 
% be sure to e.g. exclude discrepant pairs (SKS-SKKS) from your dataset 
% before running this function 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR TESTING
modelsin = 'sws_modout_domper8s.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc

% model_out: preprocessed models based on parameters
models = load(modelsin);
model_out = models.splitmods;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% struct >>> structsnew <<< have first to
% be loaded from <<< SS_models_out_layers2_dfreq0.125_stepphi5_stepdt0.25_WORKS.mat >>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO check, for what really used?
stepsizephi=2.5;
stepsizedt=0.25;

%================================================================
%================================================================
% INITIAL SETTINGS, mainly for later plotting etc.
% adjust for your needs

% plot maximum XX best models in final figure
plot_mod_max=20;
keep_mods = 500; % finally keep only the XX bests models 

% only model a specific range of the data (baz1:baz2)
modrange_low = 15;
modrange_upp = 340;
modrange_col = [219,219,219]./256;
modrange_edcol = modrange_col;

% BAZ range to plot
BAZ=0:1:360;

% MODELS plotting
color_SS_bf_1=[0.6350 0.0780 0.1840];  % color best model 
color_SS_bf_2max=[175 175 175]./256; % color others if plot_mod_max > 1      

% SPLITS plotting
myfontsize=10; 
mymarkersize=7;
fonsize_subletters=8;
linewidth_symbols=1; % width of edge of symbols
linewidth_SS=1.2;

% NULLS plotting
mymarkersize_null=6.5;
linewidth_symbols_null=linewidth_symbols; % width of edge of symbols
color_edge_null='k';
color_face_null='w';

%================================================================
%================================================================

%================================================================
% error checks
if ~isfield(model_out,'phi_eff') && ~isfield(model_out,'dt_eff') 
    error('Required fields <phi_eff> & <dt_eff> do not exist in variable <model_out> (synthetic models)! Check input struct!')
end

%================================================================
% fitting method:

disp(' ')
whichfit=input('Fitting method: [1] only phi (RMSE), [2] joint phi/dt (RMSE)');
disp(' ')

%================================================================
% measured data input

dir_res_split=dir('splitresults_PERM_FIN_KEF.txt');
dir_res_nulls=dir('splitresultsNULL_PERM_FIN_KEF.txt');
dir_res_stack=dir('KEF_stackresults.mat');

%%%%%%%%%% READ in SL data results

use_QUAL=2; % only good & fair, no query from function >>> SWS_modelling_read_data <<< appears

[RES_split, RES_nulls, RES_stack]=SWS_modelling_read_data(dir_res_split,dir_res_nulls,dir_res_stack,use_QUAL);
staname_split=RES_split(1).staname;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================
% Ask for phase results that should be plotted
phaselist_split={RES_split.phase};
phaselist_null={RES_nulls.phase};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare variables

% splits
meas_phiSC=[RES_split.phiSC_err_min; RES_split.phiSC; RES_split.phiSC_err_max]'; % create vector with [lower err phi; phi; upper err phi]
meas_phiSC(:,4)=1; % 1 in 4th column means single splits for later coloring
meas_dtSC=[RES_split.dtSC_err_min; RES_split.dtSC; RES_split.dtSC_err_max]';     % create vector with [lower err dt; dt; upper err dt]
meas_BAZ=[RES_split.baz]';
meas_BAZ_floor=floor(meas_BAZ); % round corresponding BAZs of measured SWS parameters (no big difference), 
                                % since theoretical values are only available as integers, otherwise no
                                % comparison possible!
                                
% same for the nulls                                
meas_phiSC_null=[RES_nulls.phiSC_err_min; RES_nulls.phiSC; RES_nulls.phiSC_err_max]'; 
meas_dtSC_null=[RES_nulls.dtSC_err_min; RES_nulls.dtSC; RES_nulls.dtSC_err_max]';     
meas_BAZ_null=[RES_nulls.baz]';
meas_BAZ_floor_null=floor(meas_BAZ_null);  


% same for the stacks , TODO check that only non-NULL are included   

if ~isempty(RES_stack)

    meas_phiSC_stack=[RES_stack.phiSTACK_err_min; RES_stack.phiSTACK; RES_stack.phiSTACK_err_max]'; 
    meas_phiSC_stack(:,4)=2; % 2 in 4th column means stacked splits for later coloring
    meas_dtSC_stack=[RES_stack.dtSTACK_err_min; RES_stack.dtSTACK; RES_stack.dtSTACK_err_max]';     
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

%
%================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================
% calc deviations of measured SWS parameters from theoretical distributions
% for each synthetic model

count_mods=1;

res_dt=zeros(length(meas_phiSC),1);
res_phi=res_dt;

for ii=1:length(model_out)

    curr_mod_phi=model_out(ii).phi_eff;
    curr_mod_dt=model_out(ii).dt_eff;
    curr_mod_type=model_out(ii).type;

    for kk=1:length(meas_phiSC)
        % find theoretical value for BAZ of corresponding measured value
        find_theo_phi=curr_mod_phi(meas_BAZ_floor(kk));
        find_theo_dt=curr_mod_dt(meas_BAZ_floor(kk));
       
        %calc residuum
        res_phi(kk,1)=abs(find_theo_phi-meas_phiSC(kk,2));
        if res_phi(kk,1) > 90  % if there is a residuum > 90, calc 180-res_phi to consider the flip from -90° to +90°
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

          modsall(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          modsall(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));

          modsall(count_mods).staname=staname_split;
          
         elseif  strcmp(curr_mod_type,'single_layer')

          modsall(count_mods).phi_eff=curr_mod_phi;
          modsall(count_mods).dt_eff=curr_mod_dt;
          
          modsall(count_mods).mod_type=curr_mod_type;
          
          modsall(count_mods).phi=model_out(ii).mod_paras.phi_in;
          modsall(count_mods).dt=model_out(ii).mod_paras.dt_in;
                    
          modsall(count_mods).modrange_low=modrange_low; 
          modsall(count_mods).modrange_upp=modrange_upp; 

          modsall(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          modsall(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));

          modsall(count_mods).staname=staname_split;
  
         elseif strcmp(curr_mod_type,'dipping')
    
          modsall(count_mods).phi_eff=curr_mod_phi;
          modsall(count_mods).dt_eff=curr_mod_dt;
          
          modsall(count_mods).mod_type=curr_mod_type;

          modsall(count_mods).downdipdir=model_out(ii).mod_paras.downdipdir;
          modsall(count_mods).dip=model_out(ii).mod_paras.dip;
          modsall(count_mods).thick=model_out(ii).mod_paras.thick; 
          
          %BEST_models(count_mods).fast_eff4plot=model_out(ii).mod_diplayer_fast_eff4plot;
          %BEST_models(count_mods).tlag_eff4plot=model_out(ii).mod_diplayer_tlag_eff4plot; 
          %BEST_models(count_mods).azi4plot=model_out(ii).mod_diplayer_azi4plot; 

          modsall(count_mods).modrange_low=modrange_low; 
          modsall(count_mods).modrange_upp=modrange_upp; 

          modsall(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          modsall(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));
          
          modsall(count_mods).staname=staname_split;
          
         end
          
          
  
          if whichfit==2
          % use phi and dt for joint fitting
            modsall(ii).RMS= modsall(ii).RMS_phi/90+modsall(ii).RMS_dt/4;  
            
            nameend='joint';
          
          elseif whichfit==1

          % only use phi for fitting
            modsall(ii).RMS=modsall(ii).RMS_phi/90;
            
            nameend='only_phi';

          end
          %%% only use dt for fitting
          %BEST_models(ii).RMS=BEST_models(ii).RMS_dt/4;
          
          
    count_mods=count_mods+1;
    
    


    clear count_fits
end

%%


[b,index]=sort([modsall.RMS]); % sort models based on total RMSE
modsall_sort=modsall(index); % entry 1 corresponds to minimum overall RMSE

% if ~exist('BEST_models','var')
%     warning('No models fit the measured values with the current SETTINGS (count_fits_min)!')
% end

% keep only XXX best models
modsall_sort=modsall_sort(1:keep_mods);

%================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===========================================================================================================
% plotting

% colors Liddell paper
% single splits
colorsedge(1,:)=[0 0 0];
colorsedge(2,:)=[0 0 0];

% multievent
colorsfill(1,:)=[66 91 169]./256;
colorsfill(2,:)=[0.9922    0.5529    0.2353];%./256;

%###########################################

GR_SWS_modeling_plot_results(BAZ,modsall_sort,plot_mod_max,...
    meas_BAZ_floor,meas_phiSC,meas_dtSC,meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,color_SS_bf_1,color_SS_bf_2max,linewidth_SS,colorsfill,colorsedge,...
    mymarkersize,mymarkersize_null,linewidth_symbols,linewidth_symbols_null,color_face_null,...
    color_edge_null,myfontsize,fonsize_subletters,modrange_col,modrange_edcol,meas_BAZ_floor4plot,meas_phiSC4plot,meas_dtSC4plot)

%###########################################

filename=['PLOT_modelling_DIPPING_' staname_split '_phi_dt_' num2str(plot_mod_max) '_best_models_' nameend];
    print ('-depsc', '-painters','-r600', [filename '.eps'])
    dir_eps_file=dir([filename '.eps']);
    [status,cmdout]=system(['epstopdf ' dir_eps_file.name]);

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot parameter distribtuion of models

SWS_modelling_PLOT_parameters(modsall_sort,stepsizedt,stepsizephi)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
