function modout=SWS_modeling_precomp_single_layer(stepphi, stepdt)

% generate one layer splitting models

phi = -90:stepphi:90;
dt = stepdt:stepdt:4; % avoid dt=0 s

% get all possible combinantions
comb_vecs = combvec(phi, dt);
N = length(comb_vecs);  

disp(' ')
disp(['Total number of single-layer models to generate: ' num2str(N)])
disp('Generate models...')

% preallocate structs
modout = repmat(struct('phi_eff',zeros(1,360),... 
    'dt_eff',zeros(1,360), ...
    'mod_paras', struct('phi_in',zeros(1,1), ...
    'dt_in',zeros(1,1),...
    'counter',zeros(1,1)), 'type', zeros(1,1)), N, 1);

for ii=1:N
  
    currmod=comb_vecs(:,ii); 
    modphis=currmod(1:length(currmod)/2,:);
    moddts=currmod(length(currmod)/2+1:end,:);

    modout(ii).phi_eff=ones(1,length(0:1:360)) * modphis;
    modout(ii).dt_eff=ones(1,length(0:1:360)) * moddts;
    modout(ii).mod_paras.phi_in=modphis;
    modout(ii).mod_paras.dt_in=moddts;
    modout(ii).mod_paras.counter=1;
    modout(ii).type = 'single_layer';

end

disp('Single-layer models done!')