function modout=SWS_modeling_precomp_dippinglayer(dfreq, stepdddir, stepdips, stepthick)

% generate dipping (one-) layer splitting models

%============================================================== 
%============================================================== 
% define incidence angle that should be used
inc=11;

% set up parameters
downdipdir=0:stepdddir:360; 
dips=stepdips:stepdips:75; % start with stepdips to avoid horizontal layer
thickness=stepthick:stepthick:500; % start with stepthick to avoid layer of zero thickness

% get all possible combinantions
comb_vecs=combvec(downdipdir,dips,thickness);

disp(' ')
disp(['Total number of dipping-layer models to generate: ' num2str(length(comb_vecs))])

N = length(comb_vecs);   

modout = repmat(struct('phi_eff',zeros(1,360), ...
    'dt_eff',zeros(1,360),...
    'mod_paras',struct('downdipdir', 0,...
    'dip', 0,...
    'thick', 0),...
    'type', zeros(1,1)), N, 1);

parfor ii=1:N % if problems occur, replace parfor by standard for loop
    
    currmod = comb_vecs(:,ii);
    downdipdir = currmod(1,:);
    dips = currmod(2,:);
    thickness = currmod(3,:);

    [fast_eff,tlag_eff]=...
        SWS_modeling_calc_dipping(inc, dips, downdipdir, thickness, dfreq);

     modout(ii).phi_eff = fast_eff;
     modout(ii).dt_eff = tlag_eff;
     modout(ii).mod_paras.downdipdir = downdipdir;
     modout(ii).mod_paras.dip = dips;
     modout(ii).mod_paras.thick = thickness;
     modout(ii).type = 'dipping';

end