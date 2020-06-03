function SWS_mgrund_systematic_modelling_ALL_stations_auto()


% go to each station, check for results and make modelling, based on 1,2
% and dipping layers
%
%

clear all
close all
clc

cd /home/mgrund/matlab/MSAT/examples/splitting_model


load models_out_DIPPING_layer_domper8s737317.6001.mat 

%%

% first load struct with all precalculated modells

    %======================= 
    % necessary paths
    MAIN_dir='/data_local/00_LITHOS-CAPP/00_DATA';
    SAVEDIR='/home/mgrund/matlab/MSAT/examples/splitting_model/FINAL_models/';
    
    %RES_dir_EV_DISTR='/data_local/00_LITHOS-CAPP/00_DATA/114_PLOT_TOTAL_EV_DISTR/';
    %=======================
    
    cd(MAIN_dir)
    
        cd 00_PERMANENT_STATIONS
    
        %============================
        cd 01_FINLAND
         
            dir_fol=dir('*_PERM_FIN_*');
        
                for ii=1:length(dir_fol)
                    
                    ii
                    
                    if isdir(dir_fol(ii).name)
                        
                        cd(dir_fol(ii).name)
                        
                            cd 02_results
                        
                            dir_res_split=dir('splitresults_*.txt');
                            dir_res_null=dir('splitresultsNULL_*.txt');
                            dir_res_stack=dir('*_stackresults.mat');

                            if ~isempty(dir_res_split) && length(dir_res_split) ==1

                                 if ~isempty(dir_res_stack) && length(dir_res_stack) ==1
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,dir_res_stack);
                                 else
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,[]);

                                 end
                                 
                                 % make stereoplot from best model
                                 if ~isempty(BEST_models_sort)
                                    SWS_modelling_stereoplot_THEO_dipping_4_systematic_modelling(BEST_models_sort)
                                 end
                                
                                cd(CURR_DIR)

                                cd ../..
                                                            
                            elseif ~isempty(dir_res_null) && length(dir_res_null) ==1 &&...
                                     isempty(dir_res_split) && length(dir_res_split) ==1
                                
                                 disp('NULL station (no modelling)!!!')
                                 pause(2)
                                 cd ../..
                            else
                                 
                                 warning('Problems in this folder!!!!!!!!!!!')
                                pause(2)
                                cd ../..
                            end  
                    end
                end
                
                cd ..
 %%
        %============================
        cd 02_NORWAY
         
            dir_fol=dir('*_PERM_NOR_*');

            for ii=21:length(dir_fol)
                    
                    ii
                    
                    if isdir(dir_fol(ii).name)
                        
                        cd(dir_fol(ii).name)
                        
                            cd 02_results
                        
                            dir_res_split=dir('splitresults_*.txt');
                            dir_res_null=dir('splitresultsNULL_*.txt');
                            dir_res_stack=dir('*_stackresults.mat');

                            if ~isempty(dir_res_split) && length(dir_res_split) ==1

                                 if ~isempty(dir_res_stack) && length(dir_res_stack) ==1
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,dir_res_stack);
                                 else
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,[]);

                                 end
                                 
                                 % make stereoplot from best model
                                 if ~isempty(BEST_models_sort)
                                        SWS_modelling_stereoplot_THEO_dipping_4_systematic_modelling(BEST_models_sort)
                                 end
                                
                                cd(CURR_DIR)

                                cd ../..
                                                            
                            elseif ~isempty(dir_res_null) && length(dir_res_null) ==1 &&...
                                     isempty(dir_res_split) && length(dir_res_split) ==1
                                
                                 disp('NULL station (no modelling)!!!')
                                 pause(2)
                                 cd ../..
                            else
                                 
                                 warning('Problems in this folder!!!!!!!!!!!')
                                pause(2)
                                cd ../..
                            end  
                    end
                end
                
                cd ..
            
            
            
            
            
            
            
        
   %%
   
   DDDDD=1;
   
   
        %============================
        cd 03_SWEDEN
         
            dir_fol=dir('*_PERM_SWE_*');
            
                   for ii=1:length(dir_fol)
                    
                    ii
                    
                    if isdir(dir_fol(ii).name)
                        
                        cd(dir_fol(ii).name)
                        
                            cd 02_results
                        
                            dir_res_split=dir('splitresults_*.txt');
                            dir_res_null=dir('splitresultsNULL_*.txt');
                            dir_res_stack=dir('*_stackresults.mat');

                            if ~isempty(dir_res_split) && length(dir_res_split) ==1

                                 if ~isempty(dir_res_stack) && length(dir_res_stack) ==1
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,dir_res_stack);
                                 
                                 else
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,[]);
 
                                 end
                                 
                                 % make stereoplot from best model
                                 if ~isempty(BEST_models_sort)
                                        SWS_modelling_stereoplot_THEO_dipping_4_systematic_modelling(BEST_models_sort)
                           
                                        COLLECT_BEST_models_sort(DDDDD)=BEST_models_sort(1);
                                        DDDDD=DDDDD+1;

                                 end
                                
                                cd(CURR_DIR)

                                cd ../..
                                                            
                            elseif ~isempty(dir_res_null) && length(dir_res_null) ==1 &&...
                                     isempty(dir_res_split) && length(dir_res_split) ==1
                                
                                 disp('NULL station (no modelling)!!!')
                                 pause(2)
                                 cd ../..
                            else
                                 
                                 warning('Problems in this folder!!!!!!!!!!!')
                                pause(2)
                                cd ../..
                            end  
                    end
                end
                
                cd ..
            
  
                
            %%%%%
            % make scatter plot mit allen downdipdirs etc
%  %%           
%  close all
%             downdipird=[COLLECT_BEST_models_sort.downdipdir];
%                 
%            plot(downdipird,'o')     
%                 
            
          %%  
            
            
            
        
              
                %cd ../..
            %======================================================    
            %======================================================  
  
    cd 04_SCANarray    
    
             dir_fol=dir('SA*');

                   for ii=1:length(dir_fol)
                    
                    ii
                    
                    if isdir(dir_fol(ii).name)
                        
                        cd(dir_fol(ii).name)
                        
                            cd 02_results
                        
                            dir_res_split=dir('splitresults_*.txt');
                            dir_res_null=dir('splitresultsNULL_*.txt');
                            dir_res_stack=dir('*_stackresults.mat');

                            if ~isempty(dir_res_split) && length(dir_res_split) ==1

                                 if ~isempty(dir_res_stack) && length(dir_res_stack) ==1
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,dir_res_stack);
                                 
                                 else
                                    CURR_DIR=pwd;
                                    BEST_models_sort=SWS_modelling_calc_misfit_DIPPING_layer_4_systematic(SAVEDIR,MODELLstruct,dir_res_split,dir_res_null,[]);
 
                                 end
                                 
                                 % make stereoplot from best model
                                 if ~isempty(BEST_models_sort)
                                        SWS_modelling_stereoplot_THEO_dipping_4_systematic_modelling(BEST_models_sort)
                           
                                        COLLECT_BEST_models_sort(DDDDD)=BEST_models_sort(1);
                                        DDDDD=DDDDD+1;

                                 end
                                
                                cd(CURR_DIR)

                                cd ../..
                                                            
                            elseif ~isempty(dir_res_null) && length(dir_res_null) ==1 &&...
                                     isempty(dir_res_split) && length(dir_res_split) ==1
                                
                                 disp('NULL station (no modelling)!!!')
                                 pause(2)
                                 cd ../..
                            else
                                 
                                 warning('Problems in this folder!!!!!!!!!!!')
                                pause(2)
                                cd ../..
                            end  
                    end
                end
                
                cd ..
             
             
             
             
    %%         
             
             
             
      
    

    %=====================================================                     
    cd(RES_dir_EV_DISTR)
    
        SWS_Analysis_BASICS_PLOT_total_evdistr()
    
    cd(MAIN_dir)
    %=====================================================


