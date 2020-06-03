function [RES_split, RES_nulls, RES_stack]=SWS_modelling_read_data(dir_res_split,dir_res_nulls,dir_res_stack,varargin)
%
% read SWS results for a single seismic station 
%
% OUTPUT: structs for splits and nulls based on the selected quality
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
%  created: 2016-08-22 -MG-
%  mod    : 2018-03-01 -MG-


% TODO: add automatic reading of multievent results

%=========================================================================================================== 

if isempty(varargin)

select_qual=input(['Select the qualities (for splits and nulls) you wanna get as output (default is 0): \n',...
    '    \n',...'
    '   [0] all \n',...'
    '   [1] good \n',...'
    '   [2] good & fair \n',...'
    '   [3] fair & poor \n',...'
    '   [4] fair \n',...'
    '   [5] poor']);

    if isempty(select_qual) %default
        select_qual=0;
    end

else
   select_qual=varargin{1};
    
end

%===========================================================================================================

scaling_factorALL=1; % if you wanna scale your delay time e.g. for GMT with a specific factor (same for all results)

%............................
if isempty(dir_res_split)
   disp('No file with splitresults in this folder!') 
   RES_split=[];
else
   scaling_factor=scaling_factorALL; % DEFAULT 1, adjust length of delay time vector if needed
   RES_split=read_results(dir_res_split,select_qual,scaling_factor);
end

%............................
if isempty(dir_res_nulls)
   disp('No file with NULL splitresults in this folder!') 
   RES_nulls=[];
else
   scaling_factor=scaling_factorALL; % DEFAULT 1,  adjust length of delay time vector for nulls, in general 1 is a good choice 
   RES_nulls=read_results(dir_res_nulls,select_qual,scaling_factor);
end

%............................
if isempty(dir_res_stack)
   disp('No file with surface stacked splitresults in this folder!') 
   RES_stack=[];
else
   scaling_factor=scaling_factorALL; % DEFAULT 1,  adjust length of delay time vector for nulls, in general 1 is a good choice 
   RES_stack=read_results_STACK(dir_res_stack,scaling_factor);
end
%............................



stanamecheck={};

if isempty(RES_split)
   disp(' ')
   disp('   >>> No match for your selected SPLIT quality choice!')
else
   disp(' ')
   disp('   >>> Found split results!')
   stanamecheck{end+1}=RES_split(1).staname;
end

if isempty(RES_nulls)
   disp(' ')
   disp('   >>> No match for your selected NULL quality choice!')
else
   disp(' ')
   disp('   >>> Found null results!') 
   stanamecheck{end+1}=RES_nulls(1).staname;
end

if isempty(RES_stack)
   disp(' ')
   disp('   >>> No match for stacked results!')
else
   disp(' ')
   disp('   >>> Found stack results!') 
   stanamecheck{end+1}=RES_stack(1).staname;
end

% final check for station
stanamecheck=unique(stanamecheck);

if length(stanamecheck) > 1
    error('Input files are from different stations!')
end


end


%===========================================================================================================
%===========================================================================================================
%===========================================================================================================
% subfunction 

function RES_out=read_results(dir_res,select_qual,scaling_factor)
           
fid=fopen(dir_res.name);
          
% allocate fields for speed
            
res_split.date_doy='NaN';     
res_split.staname='NaN';         
res_split.sta_lon=NaN;       % still empty          
res_split.sta_lat=NaN;       % still empty   
res_split.phase='NaN';         
res_split.baz=NaN;            
res_split.inc=NaN;           
res_split.filter=[NaN NaN];            
res_split.phiRC=NaN;            
res_split.dtRC=NaN;            
res_split.phiSC_err_min=NaN;            
res_split.phiSC=NaN;            
res_split.phiSC_err_max=NaN;                  
res_split.dtSC_err_min=NaN;                  
res_split.dtSC=NaN;                  
res_split.dtSC_err_max=NaN;                 
res_split.phiEV=NaN;                    
res_split.dtEV=NaN;                   
res_split.SNRSC=NaN;                   
res_split.quality_manual='NaN';                
res_split.NULL='NaN';
        
C = textscan(fid,'%s %s %s %f %f %s %f %f %s %f %f %f %s %f %s %f %f %s %f %s %f %f %f %f %s %s','headerlines',3);
fclose(fid);
    
%===========================================================================================================

for k=1:length(C{1})          
    res_split(k).date_doy=C{1,1}{k,1};            
    res_split(k).staname=C{1,2}{k,1}; 
    res_split(k).sta_lon=NaN;       % still empty          
    res_split(k).sta_lat=NaN;       % still empty  
    res_split(k).phase=C{1,3}{k,1};                    
    res_split(k).baz=C{1,4}(k);               
    res_split(k).inc=C{1,5}(k);                  
    res_split(k).filter=[C{1,7}(k) C{1,8}(k)];                    
    res_split(k).phiRC=C{1,10}(k);                    
    res_split(k).dtRC=C{1,11}(k);                    
    res_split(k).phiSC_err_min=C{1,12}(k);                    
    res_split(k).phiSC=C{1,14}(k);                    
    res_split(k).phiSC_err_max=C{1,16}(k);                    
    res_split(k).dtSC_err_min=C{1,17}(k);                    
    res_split(k).dtSC=C{1,19}(k)/scaling_factor;                   
    res_split(k).dtSC_err_max=C{1,21}(k);                   
    res_split(k).phiEV=C{1,22}(k);                   
    res_split(k).dtEV=C{1,23}(k);                 
    res_split(k).SNRSC=C{1,24}(k);                 
    res_split(k).quality_manual=C{1,25}{k,1};               
    res_split(k).NULL=C{1,26}{k,1};
    
    if strcmp(res_split(k).phase,'SKS')
        res_split(k).phase_col=[0    0.4470    0.7410];
    elseif strcmp(res_split(k).phase,'SKKS')
        res_split(k).phase_col=[0.8500    0.3250    0.0980];
    elseif strcmp(res_split(k).phase,'PKS')
        res_split(k).phase_col=[0.4660    0.6740    0.1880]; 
    else 
        res_split(k).phase_col=[0.9290    0.6940    0.1250];   
    end

end

%===========================================================================================================
% sort quality

sel_ev=[];

if select_qual==0 % all
   sel_ev=res_split;
elseif select_qual==1 % good
   find_ev=strcmp({res_split.quality_manual},'good');
   sel_ev=res_split(find_ev);
elseif select_qual==2 % good & fair
   find_ev=strcmp({res_split.quality_manual},'good') | strcmp({res_split.quality_manual},'fair');
   sel_ev=res_split(find_ev);
elseif select_qual==3 % fair & poor
   find_ev=strcmp({res_split.quality_manual},'fair') | strcmp({res_split.quality_manual},'poor');
   sel_ev=res_split(find_ev); 
elseif select_qual==4 % fair
   find_ev=strcmp({res_split.quality_manual},'fair');
   sel_ev=res_split(find_ev);  
elseif select_qual==5 % poor
   find_ev=strcmp({res_split.quality_manual},'poor');
   sel_ev=res_split(find_ev);    
end

RES_out=sel_ev;

end

%===========================================================================================================
%===========================================================================================================
%===========================================================================================================
% subfunction 

function RES_out=read_results_STACK(dir_res_stack,scaling_factor)
           
predat=load(dir_res_stack.name);
pre_RES_stack=predat.eqstack;

for k=1:length(pre_RES_stack)
    
    if strcmp(pre_RES_stack(k).results.stack_meth,'WS') || strcmp(pre_RES_stack(k).results.stack_meth,'RH')...
        || strcmp(pre_RES_stack(k).results.stack_meth,'nw')

        res_split(k).staname=pre_RES_stack(k).results.stnname; 
        res_split(k).sta_lon=pre_RES_stack(k).results.slat;                
        res_split(k).sta_lat=pre_RES_stack(k).results.slong;
        res_split(k).stack_meth=pre_RES_stack(k).results.stack_meth;   
        res_split(k).nsurf=pre_RES_stack(k).results.nsurf; 
        res_split(k).ndf=pre_RES_stack(k).results.ndf; 
        res_split(k).minbaz=pre_RES_stack(k).results.bazi_min; 
        res_split(k).maxbaz=pre_RES_stack(k).results.bazi_max; 
        res_split(k).meanbaz=pre_RES_stack(k).results.bazi_mean; 
        res_split(k).mindis=pre_RES_stack(k).results.dist_min; 
        res_split(k).maxdis=pre_RES_stack(k).results.dist_max; 
        res_split(k).meandis=pre_RES_stack(k).results.dist_mean; 
        res_split(k).phiSTACK_err_min=pre_RES_stack(k).results.phi_stack(1);   
        res_split(k).phiSTACK=pre_RES_stack(k).results.phi_stack(2);
        res_split(k).phiSTACK_err_max=pre_RES_stack(k).results.phi_stack(3); 
        res_split(k).dtSTACK_err_min=pre_RES_stack(k).results.dt_stack(1)/scaling_factor;
        res_split(k).dtSTACK=pre_RES_stack(k).results.dt_stack(2)/scaling_factor;   
        res_split(k).dtSTACK_err_max=pre_RES_stack(k).results.dt_stack(3)/scaling_factor;
        res_split(k).remark=pre_RES_stack(k).results.remark; 
        res_split(k).used_phases=pre_RES_stack(k).results.events_in; 
        
    end

end


RES_out=res_split;

end

%===========================================================================================================
%===============================================================
















