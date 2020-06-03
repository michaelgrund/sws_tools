function SWS_modelling_gen_models_2_layers(num_layers,dfreq)
%
% Systemcatic generation of theoretical SWS parameters for two and three layer models using equations 
% of Silver & Savage (1994) 
%
% Note: modified function MS_effective_splitting_N of MSAT toolbox (Walker & Wookey, 2012)
%       is needed to run the following code
%
% USE dominant frequency NOT period (so for 8 s use 0.125 as input for dfreq) !!! 
%
% single layer models are already included due to the parametrization,
% remove doublets at the end (TODO)
%
% PARFOR is used to parallelize process =>  speed up
% function parfor_progress(N) from matlab-exchange is NECESSARY if progress
% of calculation should be shown
%
% This function/program is part of the KaSP toolbox and free software under the GNU licencse!
%
% #####################################################################################################

%##########################################################################################################
%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
%#%                                                                                                      %#
%#%   Copyright (c) 2016 by the KaSP-Team.                                                               %#
%#%   This file is part of the Karlsruhe Seismology Processing (KaSP) Toolbox for MATLAB!                %#
%#%                                                                                                      %#
%#%   The KaSP toolbox is free software under the terms of the GNU General Public License!               %#
%#%                                                                                                      %#
%#%   Please see the copyright/licensce notice distributed together with the KaSP Toolbox!               %#
%#%                                                                                                      %#
%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
%##########################################################################################################
%
% #########################################################################################################
%
% Main author: Michael Grund (michael.grund@kit.edu)
%
%  created: 2016-07-25 -MG-
%  mod:     2017-11-29 -MG- (added second parfor loop for parallelization to speed up calculations)

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic stuff
% NOTE: for small stepsizes < stepphi > &  <stepdt > computation time increases significantly!!!

% BAZ=0:1:360; % BAZ range to calculate theoretical SWS parameters
% 
% stepphi=2.5;   
% stepdt=0.25;

stepphi=5;   
stepdt=0.15;

stepphi=45;   
stepdt=1;


% using 5° stepsize for paper of Sanz-Alonso et al. (2017)
BAZ=0:1:360; % BAZ range to calculate theoretical SWS parameters

% stepphi=25;   
% stepdt=0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model input check, 2 or 3 layers ? 

if num_layers == 2

    fast1=-90:stepphi:90;
    delaydt1=0:stepdt:4;

    fast2=-90:stepphi:90;
    delaydt2=0:stepdt:4;

    comb_vecs=combvec(fast1,fast2,delaydt1,delaydt2);
    
elseif num_layers == 3
    
%     fast1=0:stepphi:180;
%     delaydt1=0:stepdt:4;
% 
%     fast2=0:stepphi:180;
%     delaydt2=0:stepdt:4;
%     
%     fast3=0:stepphi:180;
%     delaydt3=0:stepdt:4;
    
    fast1=-90:stepphi:90;
    delaydt1=0:stepdt:4;

    fast2=-90:stepphi:90;
    delaydt2=0:stepdt:4;
    
    fast3=-90:stepphi:90;
    delaydt3=0:stepdt:4;
    
    
    

    comb_vecs=combvec(fast1,fast2,fast3,delaydt1,delaydt2,delaydt3);
end

disp(' ')
disp(['Total number of models to generate: ' num2str(length(comb_vecs))])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate models

% ZZ=1;
% dispstatus=1:1000:length(comb_vecs);

N = length(comb_vecs);   
model_out = repmat(struct('phi_eff',zeros(1,360),... 
    'dt_eff',zeros(1,360), ...
     'mod_paras', struct('phi_in',zeros(1,2), ...
    'dt_in',zeros(1,2),...
    'counter',zeros(1,1)), 'type', zeros(1,1)), N, 1 );

starttime=now;

%parfor_progress(N); % Initialize 

parfor ii=1:N

    current_model_in=comb_vecs(:,ii);

    model_phis=current_model_in(1:length(current_model_in)/2,:);
    model_dts=current_model_in(length(current_model_in)/2+1:end,:);

    
   [phi_eff_out,dt_eff_out]=inner_circle(dfreq,BAZ,model_phis,model_dts);
    
%     % modified MSAT functions are necessary, see matlab basic folder!!!
%     parfor (jj=1:length(BAZ),8) % parallel computing to speed up calculation, if problems occur use normal "for" loop      
%     %for jj=1:length(BAZ) % parallel computing to speed up calculation, if problems occur use normal "for" loop      
%         [phi_eff_out(jj),dt_eff_out(jj),a1,a2,a3,a4] = ...     
%          MS_effective_splitting_N(dfreq,BAZ(jj), ...  
%         [model_phis'],[model_dts'],'mode', 'S&S');
%     end

    
%     % check and remove NaNs
%    if sum(isnan(phi_eff_out)) == 0 && sum(isnan(dt_eff_out)) == 0
        model_out(ii).phi_eff=phi_eff_out;
        model_out(ii).dt_eff=dt_eff_out;
        model_out(ii).mod_paras.phi_in=model_phis;
        model_out(ii).mod_paras.dt_in=model_dts;
        model_out(ii).mod_paras.counter=1;
        model_out(ii).type='two_layers';
        
%         ZZ=ZZ+1;
%     end
%
%     clear phi_eff_out dt_eff_out
   

%parfor_progress; % Count 
   
   
end

%parfor_progress(0); % Clean up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save models

save(['SS_models_out_layers' num2str(num_layers) '_dfreq' num2str(dfreq) '_stepphi' num2str(stepphi) '_stepdt' num2str(stepdt) '.mat'],'model_out', '-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

%endtime=now;

% datestr(starttime)
% datestr(endtime)

% durtime=endtime-starttime;
% clc
% durtime*(24*60*60)



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%############################################################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BONF

%%% outsource second parfor loop

function [phi_eff_out,dt_eff_out]=inner_circle(dfreq,BAZ,model_phis,model_dts)


    % modified MSAT functions are necessary, see matlab basic folder!!!
    parfor (jj=1:length(BAZ),8) % parallel computing to speed up calculation, if problems occur use normal "for" loop      
    %for jj=1:length(BAZ) % parallel computing to speed up calculation, if problems occur use normal "for" loop      
        [phi_eff_out(jj),dt_eff_out(jj)] = ...     
         MS_effective_splitting_N(dfreq,BAZ(jj), ...  
        [model_phis'],[model_dts'],'mode', 'S&S');
    end


end

% EONF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% TESTS for 3 layers on GPI compute pool server


% for jj=1:length(BAZ)
% 
%         [phi_eff_out(jj),dt_eff_out(jj),a1,a2,a3,a4] = ...     
%          MS_effective_splitting_N(1/8,BAZ(jj), ...  
%         [a(1:3)'],[a(4:6)'],'mode', 'S&S');
% 
% 
% end











