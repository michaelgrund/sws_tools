function SWS_modeling_precomp_models_main()
%
% modeling of shear-wave splitting measurements 
%
% run this function to pre-compute:
%   1) single-layer models 
%   2) two-layer models
%   3) dipping layer models 
%
% MSAT package of Walker & Wookey (2012) is required and can 
% be downloaded from:
%       https://www1.gly.bris.ac.uk/MSAT/
%
% LICENSE
%
% Copyright (C) 2020  Michael Grund, Karlsruhe Institute of Technology (KIT), 
% Email: michael.grund@kit.edu
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% TERMS OF USE
%
% The modeling routines are provided "as is" and without any warranty. 
% The author cannot be held responsible for anything that happens to you 
% or your equipment. Use it at your own risk.
%
%============================================================== 
%============================================================== 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust for your needs

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

modout1 = SWS_modeling_precomp_single_layer(stepphis, stepdts);
modout2 = SWS_modeling_precomp_twolayers(1/domper, stepphim, stepdtm);
modout3 = SWS_modeling_precomp_dippinglayer(1/domper, stepdddir, stepdips, stepthick);

% merge models
splitmods = vertcat(modout1, modout2, modout3);

disp(' ')
disp('Merge models and save into file...' )

% save to mat-file
save(['sws_modout_domper' num2str(domper) 's.mat'],'splitmods', '-v7.3')

disp(' ')
disp('Model setup done!' )

%============================================================== 
%============================================================== 
