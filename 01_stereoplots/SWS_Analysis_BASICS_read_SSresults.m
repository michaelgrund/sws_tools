function RES_out=SWS_Analysis_BASICS_read_SSresults(dir_res_stack,scaling_factor)
%===============================================================================
% read StackSplit results (stacked error surfaces)
%
% Main author: Michael Grund (https://orcid.org/0000-0001-8759-2018)
% GitHub: https://github.com/michaelgrund/sws_tools
%
% created: 2018-08-22 -MG-
%===============================================================================

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

