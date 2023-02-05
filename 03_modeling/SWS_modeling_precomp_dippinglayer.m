function modout=SWS_modeling_precomp_dippinglayer(dfreq, stepdddir, stepdips, stepthick)

% generate dipping (one-) layer splitting models

%============================================================== 
%============================================================== 
% define incidence angle that should be used
inc=10;

% set up parameters
downdipdir=0:stepdddir:360; 
dips=stepdips:stepdips:75; % start with stepdips to avoid horizontal layer
thickness=stepthick:stepthick:500; % start with stepthick to avoid layer of zero thickness

% get all possible combinations
comb_vecs=combvec(downdipdir,dips,thickness);

disp(' ')
disp(['Total number of dipping-layer models to generate: ' num2str(length(comb_vecs))])
disp('Generate models...')

N = length(comb_vecs);   

modout = repmat(struct('phi_eff',zeros(1,360), ...
    'dt_eff',zeros(1,360),...
    'mod_paras',struct('downdipdir', 0,...
    'dip', 0,...
    'thick', 0, ...
    'azi4plot',zeros(1,360)),...
    'type', zeros(1,1)), N, 1);

parfor ii=1:N % if problems occur, replace parfor by standard for loop
    
    currmod = comb_vecs(:,ii);
    downdipdir = currmod(1,:);
    dips = currmod(2,:);
    thickness = currmod(3,:);

    [fast_eff,tlag_eff,azi4plot,fast4plot,tlag4plot]=...
        SWS_modeling_calc_dipping(inc, dips, downdipdir, thickness, dfreq);

     modout(ii).phi_eff = fast_eff;
     modout(ii).dt_eff = tlag_eff;
     modout(ii).mod_paras.downdipdir = downdipdir;
     modout(ii).mod_paras.dip = dips;
     modout(ii).mod_paras.thick = thickness;
     modout(ii).mod_paras.azi4plot = azi4plot;
     modout(ii).mod_paras.fast4plot = fast4plot;
     modout(ii).mod_paras.dt4plot = tlag4plot;
     modout(ii).type = 'dipping';
     

end

disp('Dipping-layer models done!')