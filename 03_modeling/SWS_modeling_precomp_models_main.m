function SWS_modeling_precomp_models_main()

% codes to be published with SWS paper GR2020

% 1) setup 1-layer models 
% 2) setup 2-layer models
% 3) setup dipping layer models 

% MSAT package of Walker & Wookey required!

%============================================================== 
%============================================================== 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust for your needs

% dominant period to be modeled
domper = 8; % in s

% settings for single layer models
stepphis = 45; % in degrees   
stepdts = 1; % in seconds

% settings for two layer models
stepphim = 45; % in degrees   
stepdtm = 1; % in seconds

% settings for dipping layer models
stepdddir = 45; % in degrees 
stepdips = 15; % in degrees 
stepthick = 100; % in km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Model setup for shear-wave splitting modeling using two-layer and dipping layer models!')

modout1 = GR_SWS_modeling_precomp_single_layer(stepphis, stepdts);
modout2 = GR_SWS_modeling_precomp_twolayers(1/domper, stepphim, stepdtm);
modout3 = GR_SWS_modeling_precomp_dippinglayer(1/domper, stepdddir, stepdips, stepthick);

% merge models
splitmods = vertcat(modout1, modout2, modout3);

disp(' ')
disp('Merge models and save into file...' )

% save to mat-file
save(['sws_modout_domper' num2str(domper) 's.mat'],'splitmods', '-v7.3')

disp(' ')
disp('Model setup done!' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


