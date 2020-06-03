function modout=GR_SWS_modeling_precomp_dippinglayer(dfreq, stepdddir, stepdips, stepthick)

% generate dipping (one-) layer splitting models

% define incidence angle that should be used
inc=11;

% this setting needs ~ 12 min to be created (~485 MB disk space)
%downdipdir=0:5:360;
%dips=0:5:75;
%thickness=0:25:500;

downdipdir=0:stepdddir:360; 
dips=stepdips:stepdips:75; % start with stepdips to avoid horizontal layer
thickness=stepthick:stepthick:500; % start with stepthick to avoid layer of zero thickness

% get all possible combinantions
comb_vecs=combvec(downdipdir,dips,thickness);

disp(' ')
disp(['Total number of dipping-layer models to generate: ' num2str(length(comb_vecs))])

N = length(comb_vecs);   
% MODELLstruct = repmat(struct('fast_eff',zeros(1,360), 'tlag_eff',zeros(1,360),...
%     'azi',zeros(1,360), 'downdipdir',zeros(1,1), 'dip',zeros(1,1),...
%     'thick',zeros(1,360), 'azi4plot',zeros(1,360), 'fast_eff4plot',zeros(1,360), ...
%     'tlag_eff4plot',zeros(1,360)), N, 1 );

modout = repmat(struct('phi_eff',zeros(1,360), ...
    'dt_eff',zeros(1,360),...
    'mod_paras',struct('downdipdir', 0,...
    'dip', 0,...
    'thick', 0),...
    'type', zeros(1,1)), N, 1);

parfor ii=1:N
    
    currmod = comb_vecs(:,ii);
    downdipdir = currmod(1,:);
    dips = currmod(2,:);
    thickness = currmod(3,:);

    [fast_eff,tlag_eff,fast_eff1,tlag_eff1,azi,azi4plot]=...
        GR_SWS_modeling_dipping_calc(inc, dips, downdipdir, thickness, dfreq);

     modout(ii).phi_eff = fast_eff;
     modout(ii).dt_eff = tlag_eff;
     modout(ii).mod_paras.downdipdir = downdipdir;
     modout(ii).mod_paras.dip = dips;
     modout(ii).mod_paras.thick = thickness;
     modout(ii).type = 'dipping';

end


