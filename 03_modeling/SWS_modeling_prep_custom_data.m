function SWS_modeling_prep_custom_data()
%
% template script to prepare your own splitting data for the modeling
% routines
%
% you can use this script especially if your splitting results are not 
% available in SplitLab/StackSplit output format 
%
% after data preparation be sure to use SWS_modeling_calc_misfit_custom
% to calculate the misfit between the pre-computed splitting models 
% (see SWS_modeling_precomp_models_main) and your data
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
%%

clear all

% give the station name for which you measured splitting parameters
staname = 'dummy';

% For each individual splitting/null measurement the following parameters are
% required. Use the structs/cell arrays to store your data. Of course you
% can wrap around some routines or loops to read your data (e.g. from textfiles etc.) 
% and to store them in the structs/cell arrays.

%============================================================
% splits

% give the phase names
split_phase = {'SKS','SKS','SKKS'};

% give the corresponding backazimuths
split_baz = [30, 45, 60];

% give the corresponding lower error bounds for phis
split_phiSC_err_min = [50, 55, 20];
% give the corresponding phi values
split_phiSC = [60, 61, 30];
% give the corresponding upper error bounds for phis
split_phiSC_err_max = [70, 75, 35];

% give the corresponding lower error bounds for dts
split_dtSC_err_min = [1.1, 0.9, 0.5];
% give the corresponding dt values
split_dtSC = [1.2, 1.0, 0.8];
% give the corresponding upper error bounds for dts
split_dtSC_err_max = [1.3, 1.4, 1.3];

%============================================================
% nulls

% conventions are the same as for the splits above
% comment the following lines if you dont have nulls
null_phase = {'SKS','PKS'};
null_baz = [90,120];

null_phiSC_err_min = [0,10];
null_phiSC = [30,35];
null_phiSC_err_max = [60,80];

null_dtSC_err_min = [2.0,3.0];
null_dtSC = [3.0,3.5];
null_dtSC_err_max = [3.5,3.8];

%============================================================
% stacks

% conventions are the same as for the splits above
% comment the following lines if you dont have stacks
stack_baz = [45,250];

stack_phi_err_min = [0,10];
stack_phi = [30,35];
stack_phi_err_max = [60,80];

stack_dt_err_min = [1.4,3.0];
stack_dt = [3.0,3.5];
stack_dt_err_max = [3.5,3.8];

%============================================================
%============================================================
% store in correct format and save files

if exist('split_baz','var')

    RES_split = [];
    RES_split.staname = staname;
    RES_split.baz = split_baz;
    RES_split.phase = split_phase;  
    RES_split.phiSC_err_min = split_phiSC_err_min; 
    RES_split.phiSC = split_phiSC;
    RES_split.phiSC_err_max = split_phiSC_err_max;  
    RES_split.dtSC_err_min = split_dtSC_err_min;
    RES_split.dtSC = split_dtSC;
    RES_split.dtSC_err_max = split_dtSC_err_max; 

    save('RES_split_cust','RES_split')
end

if exist('null_baz','var')

    RES_nulls = [];
    RES_nulls.staname = staname;
    RES_nulls.baz = null_baz;
    RES_nulls.phase = null_phase;
    RES_nulls.phiSC_err_min = null_phiSC_err_min; 
    RES_nulls.phiSC = null_phiSC;
    RES_nulls.phiSC_err_max = null_phiSC_err_max; 
    RES_nulls.dtSC_err_min = null_dtSC_err_min;
    RES_nulls.dtSC = null_dtSC;
    RES_nulls.dtSC_err_max = null_dtSC_err_max; 

    save('RES_nulls_cust','RES_nulls')
end

if exist('stack_baz','var')
    
    RES_stack = [];
    RES_stack.staname = staname;
    RES_stack.meanbaz = stack_baz; 
    RES_stack.phiSTACK_err_min = stack_phi_err_min;
    RES_stack.phiSTACK = stack_phi;
    RES_stack.phiSTACK_err_max = stack_phi_err_max;
    RES_stack.dtSTACK_err_min = stack_dt_err_min;
    RES_stack.dtSTACK = stack_dt;
    RES_stack.dtSTACK_err_max = stack_dt_err_max;

    save('RES_stack_cust','RES_stack')
end

