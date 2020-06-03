% SPLIT_MODEL.M - Example script modelling shear-wave splitting variation 
%                 with backazimuth. 
%
% This script demonstrates simple shear-wave splitting modelling using 
% the MSAT toolset. It predicts the SKS splitting variation with azimuth 
% associated with two dipping, partially-aligned olivine layers. The 
% splitting in each layer is calculated using the Christoffel equation 
% (MS_phasvels), and, by default, combined using N-layer effective 
% splitting equations (MS_effective_splitting_N; Silver and Savage, GJI, 
% 1994). Straight raypaths through the upper mantle are assumed.
%
% The example can also calculate the effective splitting parameters 
% numerically and can limit the calculation to a subset of azimuths
% (e.g. using the 'GaussianWavelet' mode in MS_effective_splitting_N).
% These features are supported by pairs of optional arguments to the 
% example script. These are any combination of:
%
%     * split_model("mode", eff_split_mode) where eff_split_mode is a 
%       string passed into the mode argument of MS_effective_splitting_N.
%       Note that setting this to 'GaussianWavelet' results in a much 
%       longer execution time for the example and the remaining options 
%       can be used to limit this.
%     * split_model("min_azi", min_azi) the minimum azimuth considered for
%       the calculation (default 0 degrees).
%     * split_model("max_azi", max_azi) the maximum azimuth considered for
%       the calculation (default 360 degrees).
%     * split_model("del_azi", del_azi) the azimuthal spacing used for the
%       calculation (default 1 degree).
%
% See also: MS_EFFECTIVE_SPLITTING_N, MS_PHASEVELS

% (C) James Wookey and Andrew Walker, 2011
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function split_model_copy(varargin)

close all

%  ** Handle optional arguments and defaults
   e_split_mode = 's&s'; % Default mode for effective splitting calc.
   min_azi = 0.0;
   max_azi = 360.0;
   del_azi = 1.0;
   iarg = 1 ;
   while iarg <= (length(varargin))
       switch lower(varargin{iarg})
           case 'mode'
              e_split_mode = lower(varargin{iarg+1}) ; 
              iarg = iarg + 2;
           case 'min_azi'
               min_azi = varargin{iarg+1};
               iarg = iarg + 2;
           case 'max_azi'
               max_azi = varargin{iarg+1};
               iarg = iarg + 2;
           case 'del_azi'
               del_azi = varargin{iarg+1};
               iarg = iarg + 2;
           otherwise
               error(['Unknown option']);
       end
   end

%%   
%  ** Setup model parameters
      S_slow = 4.814 ; % (SKS at 100 degrees from iasp91) ->
      aoi = 11.2262 ; % at 100 km depth, angle of incidence

% %  ** Layer 1 (upper) parameters
%       L1_depth = 0. ; L1_thick = 150. ; L1_dip = 30.0 ; % layer geometry
%       L1_aaz = 0 ; % a-axis azimuth (rel. to down dip direction)
%       L1_faln = 0.3; % fraction aligned
      
%  ** Layer 1 (upper) parameters
      L1_depth = 0. ; L1_thick = 150. ; L1_dip = 0.0 ; % layer geometry
      L1_aaz = 50 ; % a-axis azimuth (rel. to down dip direction)
      L1_faln = 0.8; % fraction aligned     
      
      
%  ** Layer 2 (lower) parameters
       
      % für eine downdip direction von 70° müsste für eine fast axis von
      % 50° relativ zu Nord hier -20° für aaz angegeben werden richtig? Was
      % aber passiert wenn man 2 dipping axes mit unterschiedlichen downdip
      % directions hat???
      
      L2_depth = 250. ; L2_thick = 120. ; L2_dip = 30.0 ; % layer geometry
      L2_aaz = -80.0 ;  % a-axis azimuth (rel. to down dip direction)
      L2_faln = 0.8 ; % fraction aligned 
      
      
   
  %  ** Layer 3 (lower) parameters
      L3_depth = 150. ; L3_thick = 120. ; L3_dip = 0.0 ; % layer geometry
      L3_aaz = 20.0 ;  % a-axis azimuth (rel. to down dip direction)
      L3_faln = 0.3 ; % fraction aligned    
      

     % report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
     %              L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) ;   
   

%  ** imaging parameters
      azi = [min_azi:del_azi:max_azi]; % 0 is down dip direction (perp. to
                        % strike), note that this is the seismic 
                        % backazimuth + 180 degrees. Wave is assumed to be 
                        % polarised in this direction.  
                        
      %!!!!!!!!!!!!!!!!!!!!!!!!!   
      % by -MG- , azi=0 entspreicht immer down-dip direction, davon
      % ausgehend müssen dann die korrekten BAZ-Werte berechnet werden!!!
      %
      % since here azi=0° is downdip direction (fast axis points in direction of downdip), 
      % for an estimated ~90° downdip direction like
      % in Liddell et al. (2017), one have to add 90 (subtract 270)° to the azi values to get the
      % corresponding curves over BAZ
      %!!!!!!!!!!!!!!!!!!!!!!!!!                
                        
      inc = ones(size(azi)).*90 - aoi ; % 90 is vertical


%  ** calculate distances
      [dist1]=distance_in_dipping_layer(L1_dip,aoi,L1_thick,azi) ;
      [dist2]=distance_in_dipping_layer(L2_dip,aoi,L2_thick,azi) ;
      [dist3]=distance_in_dipping_layer(L3_dip,aoi,L3_thick,azi) ;
      
      
%%      
      
%  ** load anisotropy, and generate an isotropic version of it. 
      [Cani,rh] = MS_elasticDB('olivine') ;
      [Ciso] = MS_decomp(MS_axes(Cani)) ;
      Cani = MS_rot3(Cani,90,0,0) ; % orientation for dry upper mantle

%  ** generate layer elasticities:
%        This is a Voigt-Reuss-Hill average of the appropriately rotated olivine
%        tensor and its isotropic equivalent.
      [L1_C,~] = MS_VRH([L1_faln 1-L1_faln],...
         MS_rot3(Cani,0,-L1_dip,L1_aaz,'order',[3 2 1]),rh, Ciso, rh) ;

      [L2_C,~] = MS_VRH([L2_faln 2-L2_faln],...
         MS_rot3(Cani,0,-L2_dip,L2_aaz,'order',[3 2 1]),rh, Ciso, rh) ;

           
     [L3_C,~] = MS_VRH([L3_faln 3-L3_faln],...
         MS_rot3(Cani,0,-L3_dip,L3_aaz,'order',[3 2 1]),rh, Ciso, rh) ;
     
     
%  ** interrogate elasticities to generate splitting parameters for each layer.
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L1_C, rh, inc, azi ) ;   
      fast1 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag1 = dist1./vs2' - dist1./vs1' ;
            
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L2_C, rh, inc, azi ) ;   
      fast2 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag2 = dist2./vs2' - dist2./vs1' ;

      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L3_C, rh, inc, azi ) ;   
      fast3 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag3 = dist3./vs2' - dist3./vs1' ;
      
      
      

      save('fast1','fast1')
      save('tlag1','tlag1')
      
      save('fast2','fast2')
      save('tlag2','tlag2')
      
      
%  ** calculate the effective splitting between 2 layers    
      fast_eff = zeros(size(azi)) ;
      tlag_eff = fast_eff ;
      for i = 1:length(azi)
         [fast_eff(i),tlag_eff(i)] = ...
            MS_effective_splitting_N(0.125,azi(i), ...
            [fast3(i) fast2(i) fast1(i)],[tlag3(i) tlag2(i) tlag1(i)], 'mode', e_split_mode);
        %[fast2(i) fast1(i)],[1 -0.5], 'mode', e_split_mode);
      end
      

      
      save('fast_eff','fast_eff')
      

%  ** make a figure of the results
      scrsz = get(0,'ScreenSize');
      figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;

%  ** lag times
      subplot(2,1,2)
      plot(azi,tlag1,'r--') ; hold on
      plot(azi,tlag2,'g--') ;
      plot(azi,tlag_eff,'k-','LineWidth',1.5)
      %axis([min(azi) max(azi) 0 4])
      xlabel('Polarisation (relative to downdip direction)')
      ylabel('Lag times (s)')
      legend('Upper layer','Lower layer','Total') ;
      ylim([0 4])

%  ** fast directions
      subplot(2,1,1)
      plot(azi,fast1,'r--') ; hold on
      plot(azi,fast2,'g--') ;
      plot(azi,fast_eff,'k-','LineWidth',1.5)
      %axis([min(azi) max(azi) -90 90])
      %ylim([-90 90])
      xlabel('Polarisation (relative to downdip direction)')
      ylabel('Fast shear-wave orientation (degree)')      
      legend('Upper layer','Lower layer','Total') ;
      
      
      savepdf(1,['TEST_' num2str(now)])
%savepdf(1,['TEST_'])
      
      
return
% function [dist]=distance_in_dipping_layer(dip,aoi,thick,azi)
% 
% %  ** calculate apparent dip
%       alp = atand(tand(dip) .* sind(azi-90)) ;
%       
% %  ** calculate distances
%       gam = 90 - alp - aoi ;
%       bet = 90 + alp ;
%       dist = thick .* sind(bet) ./ sind(gam) ;
% 
% return

function report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
                      L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) 

%  ** Banner
      fprintf('\n** SKS Shear-wave splitting in two layer anisotropic models.\n\n')

%  ** Report layer details
      fprintf('Layer 1: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
         L1_depth, L1_thick, L1_dip)
      fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
         L1_aaz,L1_faln)
      fprintf('\n')

      fprintf('Layer 2: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
         L2_depth, L2_thick, L2_dip)
      fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
         L2_aaz,L2_faln)

      fprintf('\n** Plotting backazimuthal variation of splitting parameters ...\n')


return