function SWS_modelling_calc_misfit_TOTAL_models(model_out,plot_mod_max)
%
% function to fit measured SWS parameters to synthetic models
% discrepant pairs (SKS-SKKS) are sorted out within the processing 
% in this function
%
%
%
%%%%%%%%%%%%
% DONE
% check for discrepant pairs from SKS-SKKS database and remove them from 
% the analysis at this point, however plot them in other color!
%%%%%%%%%%%
%
%
%

%%

clc
cd /home/mgrund/matlab/MSAT/examples/splitting_model
%cd /data_local/00_LITHOS-CAPP/02_SKS_Modelling/00000000_FINAL_modelling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% struct >>> structsnew <<< have first to
% be loaded from <<< SS_models_out_layers2_dfreq0.125_stepphi5_stepdt0.25_WORKS.mat >>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clearvars -except model_out

clearvars -except structsnew model_out MODELLstruct model_gesamt model_out_NEW
%model_out=structsnew;
%model_out=model_gesamt;
%model_out=MODELLstruct;
model_out=model_out_NEW;

stepsizephi=2.5;
stepsizedt=0.25;

plot_mod_max=50;
show_densplot=0;


% only model a specific range of the data (baz1:baz2)
modrange_low=15;
modrange_upp=340;
modrange_col=[219,219,219]./256;
modrange_edcol=modrange_col;

% model_out: preprocessed models based on parameters
% mod_max: plot maximum XX models in final figure

%================================================================
%#################################
% INITIAL SETTINGS

% number of minimum fitted points to plot model in the end
count_fits_min=1; 

% BAZ range 2 plot
BAZ=0:1:360;

%#################################
% MODEL plotting 

color_SS_bf_1=[0.6350    0.0780    0.1840];  % best model 
color_SS_bf_2max=[175 175 175]./256;         % others if plot_mod_max > 1      

%#################################
% SPLITS plotting

myfontsize=10; 
mymarkersize=7;
fonsize_subletters=8;
linewidth_symbols=1; % width of edge of symbols
linewidth_SS=1.2;

color_edge=[0 51 102]./256;
color_face=[0    0.4470    0.7410];
color_error=[0 51 102]./256;

%#################################
% NULLS plotting

mymarkersize_null=6.5;
linewidth_symbols_null=linewidth_symbols; % width of edge of symbols

color_edge_null='k';
color_face_null='w';

%#################################

%================================================================
% check input
% if ~isfield(model_out,'phi_eff') && ~isfield(model_out,'dt_eff') 
%     error('Required fields <phi_eff> & <dt_eff> do not exist in variable <model_out> (synthetic models)! Check input struct!')
% end

if ~isfield(model_out,'phi_eff') && ~isfield(model_out,'dt_eff') 
    error('Required fields <phi_eff> & <dt_eff> do not exist in variable <model_out> (synthetic models)! Check input struct!')
end


%================================================================
% fitting method:

disp(' ')
whichfit=input('Fitting method: [1] only phi (RMS), [2] joint phi/dt (RMS)');
disp(' ')
%================================================================
%
% disp all result files that were found in the current folder

% dir_splits=dir('splitresults_*.txt');
% 
% if isempty(dir_splits)
%    disp(' ')
%    disp('  #####################################################')
%    disp('  No < splitresults_*.txt > files found in that folder!')
%    disp('  #####################################################')
%    return
% end
% 
% disp(' ')
% disp('START modelling your SWS data...')
% disp(' ')
% disp('Found results files: ')
% disp(' ')
% 
% % modelling starts only if at least one single splitting measurement was
% % made
% 
% for ii=1:length(dir_splits)
%     disp(['   [' num2str(ii) '] ' dir_splits(ii).name])
% end
% 
% disp(' ')
% wfile=input('Which file you wanna use (Enter number [XX] from the given list)?:');
% disp(' ')
% 
% if isempty(wfile)
%     error('You selected a non-available entry (potentially ENTER) ;) Try again!')
% elseif ~ismember(wfile,[1:length(dir_splits)])
%     error('You selected a non-available entry ;) Try again!')
% end

% dir_res_split=dir('splitresults_PERM_SWE_GRAU.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_SWE_GRAU.txt');
% dir_res_stack=dir('GRAU_stackresults.mat');



dir_res_split=dir('splitresults_PERM_FIN_RAF.txt');
dir_res_nulls=dir('splitresultsNULL_PERM_FIN_RAF.txt');
dir_res_stack=dir('RAF_stackresults.mat');


dir_res_split=dir('splitresults_PERM_FIN_KEF.txt');
dir_res_nulls=dir('splitresultsNULL_PERM_FIN_KEF.txt');
dir_res_stack=dir('KEF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_KAF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_KAF.txt');
% dir_res_stack=dir('KAF_stackresults.mat');
% 
% dir_res_split=dir('splitresults_PERM_FIN_OUL.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_OUL.txt');
% dir_res_stack=dir('OUL_stackresults.mat');
% 
% dir_res_split=dir('splitresults_PERM_RUS_LVZ.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_RUS_LVZ.txt');
% dir_res_stack=dir('LVZ_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_RNF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_RNF.txt');
% dir_res_stack=dir('RNF_stackresults.mat');
% 
% dir_res_split=dir('splitresults_PERM_FIN_SGF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_SGF.txt');
% dir_res_stack=dir('SGF_stackresults.mat');
% 
% dir_res_split=dir('splitresults_SA64.txt');
% dir_res_nulls=dir('splitresultsNULL_SA64.txt');
% dir_res_stack=dir('SA64_stackresults.mat');
% 
% dir_res_split=dir('splitresults_SA19.txt');
% dir_res_nulls=dir('splitresultsNULL_SA19.txt');
% dir_res_stack=dir('SA19_stackresults.mat');




% dir_res_split=dir('splitresults_PERM_FIN_FIA1.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_FIA1.txt');
% dir_res_stack=dir('FIA1_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_SWE_FINU.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_SWE_FINU.txt');
% dir_res_stack=dir('FINU_stackresults.mat');
% 
% dir_res_split=dir('splitresults_PERM_LAT_SLIT.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_LAT_SLIT.txt');
% dir_res_stack=dir('SLIT_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_JOF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_JOF.txt');
% dir_res_stack=dir('JOF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_VAF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_VAF.txt');
% dir_res_stack=dir('VAF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_PVF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_PVF.txt');
% dir_res_stack=dir('PVF_stackresults.mat');

%%%%%%%%%% READ in SL data results

use_QUAL=2; % only good & fair, no query from function >>> SWS_modelling_read_data <<< appears

[RES_split, RES_nulls, RES_stack]=SWS_modelling_read_data(dir_res_split,dir_res_nulls,dir_res_stack,use_QUAL);
staname_split=RES_split(1).staname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out discrepant SKS-SKKS pairs

clc

[RES_split,RES_nulls]=SWS_modelling_sort_out_disc_SKS_SKKS(RES_split,RES_nulls);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================%========================================================================
%========================================================================
%========================================================================
%========================================================================
% if two stations should be plotted together and modelled together
% uncomment the following code
%========================================================================
%========================================================================
%========================================================================
%========================================================================%========================================================================

% dir_res_split=dir('splitresults_PERM_FIN_FIA1.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_FIA1.txt');
% dir_res_stack=dir('FIA1_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_KEF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_KEF.txt');
% dir_res_stack=dir('KEF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_KAF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_KAF.txt');
% dir_res_stack=dir('KAF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_JOF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_JOF.txt');
% dir_res_stack=dir('JOF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_VAF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_VAF.txt');
% dir_res_stack=dir('VAF_stackresults.mat');

% dir_res_split=dir('splitresults_PERM_FIN_PVF.txt');
% dir_res_nulls=dir('splitresultsNULL_PERM_FIN_PVF.txt');
% dir_res_stack=dir('PVF_stackresults.mat');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% READ in SL data results
% 
% use_QUAL=2; % only good & fair
% 
% [RES_split2, RES_nulls2, RES_stack2]=SWS_modelling_read_data(dir_res_split,dir_res_nulls,dir_res_stack,use_QUAL);
% staname_split=RES_split2(1).staname;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sort out discrepant SKS-SKKS pairs
% 
% [RES_split2,RES_nulls2]=SWS_modelling_sort_out_disc_SKS_SKKS(RES_split2,RES_nulls2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % merge the results of the two stations
% RES_splitTOT=horzcat(RES_split,RES_split2);
% RES_nullsTOT=horzcat(RES_nulls,RES_nulls2);
% 
% RES_split=RES_splitTOT;
% RES_nulls=RES_nullsTOT;

%================================================================

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

res_split_all=[RES_split.phiSC];

%================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================

% [1] place panels for phi and dt beneath each other
% [2] place panels beside each other (EGU 2017 poster style)

% DEFAULT is [1]

clc
disp(' ')
wformat=input('Which format?: [1] place vertical [2] place horizontal');
disp(' ')

if isempty(wformat)
   wformat=1;
end

%
%================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================

% RES_split=1
% RES_nulls=1

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

% %
% 
% 
% %average results with a gliding window in 5° steps
% 
% 
% windowcenters=[2:2:358];
% takevalmeanphi=nan(length(windowcenters),4);
% takevalmeandt=nan(length(windowcenters),3);
% 
% for ii=1:length(windowcenters)
%     
%     currwin=[windowcenters(ii)-2 windowcenters(ii)+2];
%     
%     GG=1;
%     for jj=1:length(meas_BAZ_floor)
%         
%         if currwin(1) < meas_BAZ_floor(jj) && currwin(2) > meas_BAZ_floor(jj)
%             takevalphi(GG)=meas_phiSC(jj,2);
%             takevaldt(GG)=meas_dtSC(jj,2);
%             GG=GG+1;
%         end
% 
%     end
%     
%     %
%     takevalmeanphi(ii,4)=1;
%     if exist('takevalphi','var') && length(takevalphi) > 1
%         takevalmeanphi(ii,:)=[1 mean(takevalphi) 1 1];
%         takevalmeandt(ii,:)=[1 mean(takevaldt) 1];
%     end
% 
%     clear takevalphi
% 
% end
% 
% % close all
% % plot(takevalmeanphi,'o')
% % ylim([-90 90])
% 
% 
% 
% 
% 
% 
%  % rename original variables
%  meas_phiSC=takevalmeanphi;
%  meas_dtSC=takevalmeandt;
%  meas_BAZ_floor=windowcenters;
%  
%  % for plotting only
%  meas_phiSC4plot=takevalmeanphi;
%  meas_dtSC4plot=takevalmeandt;
%  meas_BAZ_floor4plot=windowcenters;
% 
% %

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
    
    count_fits=0;
    
    curr_mod_phi=model_out(ii).phi_eff;
    curr_mod_dt=model_out(ii).dt_eff;
    curr_mod_type=model_out(ii).mod_type;
    

%     for jj=1:length(meas_phiSC)
%             % find theoretical value for BAZ of corresponding measured value
%             find_theo_phi=curr_mod_phi(meas_BAZ_floor(jj)); 
%             find_theo_dt=curr_mod_dt(meas_BAZ_floor(jj));
% 
%             % only if measured phi and dt fit theoretical values correct within errorbounds, result is counted as ok
%             if find_theo_phi >= meas_phiSC(jj,1) && find_theo_phi <= meas_phiSC(jj,3) && find_theo_dt >= meas_dtSC(jj,1) && find_theo_dt <= meas_dtSC(jj,3)
%                 count_fits=count_fits+1;
%             end
%     end
%     
%     
%         if count_fits >= count_fits_min
%           BEST_models(count_mods).phi_eff=curr_mod_phi;
%           BEST_models(count_mods).dt_eff=curr_mod_dt;
%           BEST_models(count_mods).meas_vals=length(res_split_all);
%           BEST_models(count_mods).corr_vals=count_fits;
%           BEST_models(count_mods).corr_vals_perc=(100/length(res_split_all))*count_fits;
%           
%           BEST_models(count_mods).phi=model_out(ii).phi_in;
%           BEST_models(count_mods).dt=model_out(ii).dt_in;
%           
% %           BEST_models(AA).fast=current_model_phi;
% %           BEST_models(AA).dt=current_model_dt;
%           count_mods=count_mods+1;
%     end
    
    
    
    
    
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
    
           
         if strcmp(curr_mod_type,'2_layer')
             
             
          BEST_models(count_mods).phi_eff=curr_mod_phi;
          BEST_models(count_mods).dt_eff=curr_mod_dt;
          
          BEST_models(count_mods).mod_type=curr_mod_type;
          
          BEST_models(count_mods).phi=model_out(ii).mod_2layer_phi_in;
          BEST_models(count_mods).dt=model_out(ii).mod_2layer_dt_in;
                    
          BEST_models(count_mods).modrange_low=modrange_low; 
          BEST_models(count_mods).modrange_upp=modrange_upp; 

          BEST_models(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          BEST_models(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));

          BEST_models(count_mods).staname=staname_split;
          
         elseif  strcmp(curr_mod_type,'1_layer')

          BEST_models(count_mods).phi_eff=curr_mod_phi;
          BEST_models(count_mods).dt_eff=curr_mod_dt;
          
          BEST_models(count_mods).mod_type=curr_mod_type;
          
          BEST_models(count_mods).phi=model_out(ii).mod_1layer_phi_in;
          BEST_models(count_mods).dt=model_out(ii).mod_1layer_dt_in;
                    
          BEST_models(count_mods).modrange_low=modrange_low; 
          BEST_models(count_mods).modrange_upp=modrange_upp; 

          BEST_models(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          BEST_models(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));

          BEST_models(count_mods).staname=staname_split;
  
         elseif strcmp(curr_mod_type,'dip_layer')
    
          BEST_models(count_mods).phi_eff=curr_mod_phi;
          BEST_models(count_mods).dt_eff=curr_mod_dt;
          
          BEST_models(count_mods).mod_type=curr_mod_type;

          BEST_models(count_mods).downdipdir=model_out(ii).mod_diplayer_downdipdir;
          BEST_models(count_mods).dip=model_out(ii).mod_diplayer_dip;
          BEST_models(count_mods).thick=model_out(ii).mod_diplayer_thick; 
          
          BEST_models(count_mods).fast_eff4plot=model_out(ii).mod_diplayer_fast_eff4plot;
          BEST_models(count_mods).tlag_eff4plot=model_out(ii).mod_diplayer_tlag_eff4plot; 
          BEST_models(count_mods).azi4plot=model_out(ii).mod_diplayer_azi4plot; 

          BEST_models(count_mods).modrange_low=modrange_low; 
          BEST_models(count_mods).modrange_upp=modrange_upp; 

          BEST_models(count_mods).RMS_phi=sqrt(sum(res_phi.^2)/length(meas_phiSC));
          BEST_models(count_mods).RMS_dt=sqrt(sum(res_dt.^2)/length(meas_phiSC));
          
          BEST_models(count_mods).staname=staname_split;
          
         end
          
          
  
          if whichfit==2
          % use phi and dt for joint fitting
            BEST_models(ii).RMS= BEST_models(ii).RMS_phi/90+BEST_models(ii).RMS_dt/4;  
            
            nameend='joint';
          
          elseif whichfit==1

          % only use phi for fitting
            BEST_models(ii).RMS=BEST_models(ii).RMS_phi/90;
            
            nameend='only_phi';

          end
          %%% only use dt for fitting
          %BEST_models(ii).RMS=BEST_models(ii).RMS_dt/4;
          
          
    count_mods=count_mods+1;
    
    


    clear count_fits
end
%
% if ~exist('BEST_models','var')
%     warning('No models fit the measured values with the current SETTINGS (count_fits_min)!')
%     BEST_models=[];
%     return
% end
% 
% %[b,index]=sort([BEST_models.corr_vals]);
% [b,index]=sort([BEST_models.RMS]);
% BEST_models_sort=BEST_models(index);
% BEST_models_sort=fliplr(BEST_models_sort);
%%

[b,index]=sort([BEST_models.RMS]);
BEST_models_sort=BEST_models(index);
clear BEST_models b index

% if ~exist('BEST_models','var')
%     warning('No models fit the measured values with the current SETTINGS (count_fits_min)!')
% end

BEST100=BEST_models_sort(1:500);
clear BEST_models_sort
BEST_models_sort=BEST100;

% BEST900=BEST_models_sort(1:900);
% clear BEST_models_sort
% BEST_models_sort=BEST900;
% clear BEST900


%



%================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================
% TODO revise

%===========================================================================================================
% plotting

% colorsedge(1,:)=[0 51 102]./256;
% colorsedge(2,:)=[139,26,26]./256;
	
% colorsfill(1,:)=[0    0.4470    0.7410];
% colorsfill(2,:)=[238,44,44]./256;

% new vega colors 
% % single splits
% colorsedge(1,:)=[0.1922    0.6392    0.3294];
% colorsedge(2,:)=[0.9020    0.3333    0.0510];
% 
% % multievent
% colorsfill(1,:)=[0.4549    0.7686    0.4627];
% colorsfill(2,:)=[0.9922    0.5529    0.2353];

% colors Liddell paper
% single splits
colorsedge(1,:)=[0 0 0];
colorsedge(2,:)=[0 0 0];

% multievent
colorsfill(1,:)=[66 91 169]./256;
colorsfill(2,:)=[0.9922    0.5529    0.2353];%./256;






%###########################################
if wformat==1

 SWS_modelling_PLOT_results_vert(BAZ,BEST_models_sort,plot_mod_max,...
    meas_BAZ_floor,meas_phiSC,meas_dtSC,meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,color_SS_bf_1,color_SS_bf_2max,linewidth_SS,colorsfill,colorsedge,...
    mymarkersize,mymarkersize_null,linewidth_symbols,linewidth_symbols_null,color_face_null,...
    color_edge_null,myfontsize,fonsize_subletters,modrange_col,modrange_edcol,meas_BAZ_floor4plot,meas_phiSC4plot,meas_dtSC4plot)

elseif wformat==2 % TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  % TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  % TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
 SWS_modelling_PLOT_results_horz()
    
else
    error('Wrong input!')
end
%###########################################



filename=['PLOT_modelling_DIPPING_' staname_split '_phi_dt_' num2str(plot_mod_max) '_best_models_' nameend];
    print ('-depsc', '-painters','-r600', [filename '.eps'])
    dir_eps_file=dir([filename '.eps']);
    [status,cmdout]=system(['epstopdf ' dir_eps_file.name]);

% filename='test_forward_model_2';
%     print ('-dpng', '-painters','-r600', [filename '.png'])

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot parameter distribtuion of models

SWS_modelling_PLOT_parameters(BEST_models_sort,stepsizedt,stepsizephi)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

