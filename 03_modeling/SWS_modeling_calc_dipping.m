function [fast_eff,tlag_eff]=SWS_modeling_calc_dipping(incval,dips,downdipdir,thick,dfreq)

% calculate splitting parameters for dipping layer models
%
% this function is mainly based on the original MSAT 
% function SPLIT_MODEL.M and modified to fit into the workflow

% MSAT package of Walker & Wookey required!

%============================================================== 
%============================================================== 
%  ** Handle optional arguments and defaults
   e_split_mode = 's&s'; % Default mode for effective splitting calc.
   min_azi = 0.0;
   max_azi = 360.0;
   del_azi = 1.0;
   
%==============================================================    
%  ** Setup model parameters
      aoi=incval;

%==============================================================  
%  ** Layer parameters
      L1_depth = 0; 
      L1_thick = thick ; 
      L1_dip = dips; % layer geometry
      L1_aaz =0; % a-axis azimuth (rel. to down dip direction)
      L1_faln = 0.3; % fraction aligned
      
%==============================================================     
%  ** imaging parameters
      azi = [min_azi:del_azi:max_azi]-downdipdir-180; % 0 is down dip direction (perp. to
                        % strike), note that this is the seismic 
                        % backazimuth + 180 degrees. Wave is assumed to be 
                        % polarised in this direction.  
                 
      % since here azi=0° is downdip direction (fast axis points in direction of downdip), 
      % for an estimated ~90° downdip direction one have to add 90 (subtract 270)° 
      % to the azi values to get the corresponding curves over BAZ   
      inc = ones(size(azi)).*90 - aoi ; % 90 is vertical (set aoi==0)

%  ** calculate distances
      [dist1]=distance_in_dipping_layer(L1_dip,aoi,L1_thick,azi) ;

%==============================================================       
%  ** load anisotropy, and generate an isotropic version of it. 
      [Cani,rh] = MS_elasticDB('olivine') ; 
      [Ciso] = MS_decomp(MS_axes(Cani)) ;
      Cani = MS_rot3(Cani,90,0,0) ; % orientation for dry upper mantle

%  ** generate layer elasticities:
%     This is a Voigt-Reuss-Hill average of the appropriately rotated olivine
%     tensor and its isotropic equivalent.
      [L1_C,~] = MS_VRH([L1_faln 1-L1_faln],...
         MS_rot3(Cani,0,-L1_dip,L1_aaz,'order',[3 2 1]),rh, Ciso, rh) ;

%  ** interrogate elasticities to generate splitting parameters for layer.
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L1_C, rh, inc, azi) ;   
      fast1 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag1 = dist1./vs2' - dist1./vs1' ;
                  
%  ** calculate the effective splitting for the (dipping) layer  
      fast_eff = zeros(size(azi)) ;
      tlag_eff = fast_eff ;
      for i = 1:length(azi)
         [fast_eff(i),tlag_eff(i)] = ...
            MS_effective_splitting_N(dfreq,azi(i), ...
            fast1(i),tlag1(i), 'mode', e_split_mode);
      end

%==============================================================        
      % recalculate azi to bazi for correct plotting over bazi between 0
      % and 360 degrees
      azi = min_azi:del_azi:max_azi;

      % recalculate fast_eff
      for ii=1:length(azi)
          if (fast_eff(ii)+downdipdir) > 90
              fast_eff(ii)=fast_eff(ii)+downdipdir-180;
          elseif (fast_eff(ii)+downdipdir) < -90
              fast_eff(ii)=fast_eff(ii)+downdipdir+180;
          else
               fast_eff(ii)=fast_eff(ii)+downdipdir;
          end
      end

      % now correctly sort between +-90°
      for ii=1:length(azi)
          if (fast_eff(ii)) > 90
              fast_eff(ii)=fast_eff(ii)-180;
          elseif (fast_eff(ii)) < -90
              fast_eff(ii)=fast_eff(ii)+180;
          else
               fast_eff(ii)=fast_eff(ii);
          end
      end

return

%==============================================================  
%==============================================================  

function [dist]=distance_in_dipping_layer(dip,aoi,thick,azi)

%  ** calculate apparent dip
      alp = atand(tand(dip) .* sind(azi-90)) ;
      
%  ** calculate distances
      gam = 90 - alp - aoi ;
      bet = 90 + alp ;
      dist = thick .* sind(bet) ./ sind(gam) ;

return
