function vers_out=SWS_modeling_check_matlab_version()
%==========================================================================

vers=version('-release');

vers_yyyy=str2double(vers(1:4));
vers_let=vers(5);

if vers_yyyy > 2020 || (vers_yyyy == 2020 && strcmp(vers_let,'a')) % MATLAB 2020a or higher
   vers_out=1;
else
    vers_out=0;
end

%==================================================================================================================================
%==================================================================================================================================
% EOF