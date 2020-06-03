function modout=GR_SWS_modeling_precomp_twolayers(dfreq, stepphi, stepdt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate two-layer splitting models

% basic stuff
% NOTE: for small stepsizes < stepphi > &  <stepdt > computation time increases significantly!!!

BAZ=0:1:360; % BAZ range to calculate theoretical SWS parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model input check, 2 or 3 layers ? 

phi1=-90:stepphi:90;
dt1=stepdt:stepdt:4; % avoid dt=0 s

phi2=-90:stepphi:90;
dt2=stepdt:stepdt:4; % avoid dt=0 s

% get all possible combinantions
comb_vecs=combvec(phi1,phi2,dt1,dt2);
    
disp(' ')
disp(['Total number of two-layer models to generate: ' num2str(length(comb_vecs))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate models

N = length(comb_vecs);  
% preallocate structs
modout = repmat(struct('phi_eff',zeros(1,360),... 
    'dt_eff',zeros(1,360), ...
    'mod_paras', struct('phi_in',zeros(1,2), ...
    'dt_in',zeros(1,2),...
    'counter',zeros(1,1)), 'type', zeros(1,1)), N, 1 );

parfor ii=1:N
    
    currmod_in = comb_vecs(:,ii);
    modphis = currmod_in(1:length(currmod_in)/2,:);
    moddts = currmod_in(length(currmod_in)/2+1:end,:);
    [phi_eff_out, dt_eff_out] = nest(dfreq, BAZ, modphis, moddts);

    modout(ii).phi_eff = phi_eff_out;
    modout(ii).dt_eff = dt_eff_out;
    modout(ii).mod_paras.phi_in = modphis;
    modout(ii).mod_paras.dt_in = moddts;
    modout(ii).mod_paras.counter = 1;
    modout(ii).type = 'two_layers';
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%############################################################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BONF

%%% outsource second parfor loop

function [phi_eff_out, dt_eff_out]=nest(dfreq, BAZ, modphis, moddts)

    parfor (jj=1:length(BAZ), 8) % parallel computing to speed up calculation
    %for jj=1:length(BAZ) % uncomment, if problems occur with parfoor      
        [phi_eff_out(jj), dt_eff_out(jj)] = ...     
         MS_effective_splitting_N(dfreq, BAZ(jj), ...  
         modphis', moddts','mode', 'S&S');
    end

end

% EONF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
