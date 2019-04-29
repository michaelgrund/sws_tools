function [RES_split, RES_nulls]=SWS_Analysis_BASICS_read_SLresults(varargin)
%===============================================================================
% read SplitLab (SL) results for a single seismic station 
%
% output: structs for splits and nulls based on the selected quality
%
% Main author: Michael Grund (michael.grund@kit.edu)
%
% created: 2016-08-22 -MG-
%===============================================================================

%clear all

if isempty(varargin)

select_qual=input(['Select the qualities you wanna get as output (default is 0): \n',...
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

if ~isempty(varargin) 

    if length(varargin) > 1
        dir_res_split=varargin{2};
        dir_res_nulls=varargin{3};
    else
        dir_res_split=dir('splitresults_*.txt');
        dir_res_nulls=dir('splitresultsNULL_*.txt'); 
    end

else
    
   dir_res_split=dir('splitresults_*.txt');
   dir_res_nulls=dir('splitresultsNULL_*.txt'); 
    
end

if isempty(dir_res_split)
   disp('No file with splitresults in this folder!') 
   RES_split=[];
else
   scaling_factor=1; % DEFAULT 1, adjust length of delay time vector if needed
   RES_split=read_results(dir_res_split,select_qual,scaling_factor);
end

if isempty(dir_res_nulls)
   disp('No file with NULL splitresults in this folder!') 
   RES_nulls=[];
else
   scaling_factor=1; % DEFAULT 1,  adjust length of delay time vector for nulls, in general 1 is a good choice 
   RES_nulls=read_results(dir_res_nulls,select_qual,scaling_factor);
end

if ~isempty(RES_split)
    disp(' ')
    disp(['This is station ' RES_split(1).staname])
elseif ~isempty(RES_nulls)
    disp(' ')
    disp(['This is station ' RES_nulls(1).staname])
else
    RES_split=[];
    RES_nulls=[];
    return
end

if isempty(RES_split)
   disp(' ')
   disp('   >>> No match for your selected SPLIT quality choice!')
   pause(2)
end

if isempty(RES_nulls)
   disp(' ')
   disp('   >>> No match for your selected NULL quality choice!')
   pause(2)
end

end

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
%===========================================================