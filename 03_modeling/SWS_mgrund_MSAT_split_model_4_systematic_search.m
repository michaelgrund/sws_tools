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

function [fast_eff,tlag_eff,fast_eff1,tlag_eff1,azi,azi4plot]=SWS_mgrund_MSAT_split_model_4_systematic_search(INCvals,dips,downdipdir,thickness)




close all

%  ** Handle optional arguments and defaults
   %e_split_mode = 'GaussianWavelet'%,'s&s'; % Default mode for effective splitting calc.
   e_split_mode = 's&s'; % Default mode for effective splitting calc.
   min_azi = 0.0;
   max_azi = 360.0;
   del_azi = 1.0;
   iarg = 1 ;
%    while iarg <= (length(varargin))
%        switch lower(varargin{iarg})
%            case 'mode'
%               e_split_mode = lower(varargin{iarg+1}) ; 
%               iarg = iarg + 2;
%            case 'min_azi'
%                min_azi = varargin{iarg+1};
%                iarg = iarg + 2;
%            case 'max_azi'
%                max_azi = varargin{iarg+1};
%                iarg = iarg + 2;
%            case 'del_azi'
%                del_azi = varargin{iarg+1};
%                iarg = iarg + 2;
%            otherwise
%                error(['Unknown option']);
%        end
%    end


  
%==============================================================    
%==============================================================  
% original download setup
% 
%    %  ** Setup model parameters
%       S_slow = 4.814 ; % (SKS at 100 degrees from iasp91) ->
%       aoi = 11.2262 ; % at 100 km depth, angle of incidence
% 
% %  ** Layer 1 (upper) parameters
%       L1_depth = 0. ; L1_thick = 100. ; L1_dip = 0.0 ; % layer geometry
%       L1_aaz = 30.0 ; % a-axis azimuth (rel. to down dip direction)
%       L1_faln = 0.3 ; % fraction aligned
%       
% %  ** Layer 2 (lower) parameters
%       L2_depth = 150. ; L2_thick = 60. ; L2_dip = 30.0 ; % layer geometry
%       L2_aaz = 0.0 ;  % a-axis azimuth (rel. to down dip direction)
%       L2_faln = 0.3 ; % fraction aligned 
%==============================================================     
%==============================================================      
 % own setting tests 

%  ** Setup model parameters
      S_slow = 4.814 ; % (SKS at 100 degrees from iasp91) ->
      aoi =11.2262 ; % at 100 km depth, angle of incidence
      
      aoi=INCvals;

%  ** Layer 1 (upper) parameters

      %downdipdir=230;

      L1_depth = 0. ; L1_thick = thickness ; L1_dip = dips ; % layer geometry
      L1_aaz =0 ; % a-axis azimuth (rel. to down dip direction)
      L1_faln = 0.3; % fraction aligned
      
%  ** Layer 2 (lower) parameters
      L2_depth = 0. ; L2_thick = 0 ; L2_dip = 0.0 ; % layer geometry
      L2_aaz =0 %80 ;  % a-axis azimuth (rel. to down dip direction)
      L2_faln = 0 ; % fraction aligned 

%  ** Layer 3 (lower) parameters
      L3_depth = 150. ; L3_thick = 200 ; L3_dip = 0.0 ; % layer geometry
      L3_aaz =70 %80 ;  % a-axis azimuth (rel. to down dip direction)
      L3_faln = 0.6 ; % fraction aligned 
      
      
%==============================================================     
%============================================================== 



      report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
                   L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) ; 

%  ** imaging parameters
      azi = [min_azi:del_azi:max_azi]-downdipdir-180; % 0 is down dip direction (perp. to
                        % strike), note that this is the seismic 
                        % backazimuth + 180 degrees. Wave is assumed to be 
                        % polarised in this direction.  
                 
                        
      %!!!!!!!!!!!!!!!!!!!!!!!!!   
      % by -MG-
      % since here azi=0° is downdip direction (fast axis points in direction of downdip), 
      % for an estimated ~90° downdip direction like
      % in Liddell et al. (2017), one have to add 90 (subtract 270)° to the azi values to get the
      % corresponding curves over BAZ
      %!!!!!!!!!!!!!!!!!!!!!!!!!                
                        
      inc = ones(size(azi)).*90 - aoi ; % 90 is vertical (set aoi==0)


%  ** calculate distances
      [dist1]=distance_in_dipping_layer(L1_dip,aoi,L1_thick,azi) ;
      [dist2]=distance_in_dipping_layer(L2_dip,aoi,L2_thick,azi) ;
      
 %%     
      
 %fayalite , albite
 
      
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
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L1_C, rh, inc, azi) ;   
      fast1 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag1 = dist1./vs2' - dist1./vs1' ;
      
      
      
      %%% when aio == 0 then all distances are the same
      %tlag1 = 0.3.*ones(length(dist1)); %% MG KAF settings
      
            
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L2_C, rh, inc, azi) ;   
      fast2 = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
      tlag2 = dist2./vs2' - dist2./vs1' ;
      
      %%% when aio == 0 then all distances are the same
      %tlag2 = 0.1.*ones(length(dist2));
      
      
      [ pol, ~, vs3, vs3, ~, ~, ~ ] = MS_phasevels( L3_C, rh, inc, azi ) ; 
      fast3 = MS_unwind_pm_90((azi+pol')); 
      tlag3 = 1.*ones(length(dist2));
      
%  ** calculate the effective splitting between 2 layers    
      fast_eff = zeros(size(azi)) ;
      tlag_eff = fast_eff ;
      for i = 1:length(azi)
         [fast_eff(i),tlag_eff(i)] = ...
            MS_effective_splitting_N(0.125,azi(i), ...
            [fast1(i)],[tlag1(i)], 'mode', e_split_mode);
        
      end
      
      
      
%       for i = 1:length(azi)
%          [fast_eff(i),tlag_eff(i)] = ...
%             MS_effective_splitting_N(0.125,azi(i), ...
%             [fast3(i) fast2(i) fast1(i)],[tlag3(i) tlag2(i) tlag1(i)], 'mode', e_split_mode);
%         
%       end
      
% save parameters for theoretical steroplot, see function <<<
% SWS_modelling_stereoplot_THEO_dipping >>>

%fast_eff=fast_eff+90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is used to make the synthetic stereoplot

fast_eff1=fast_eff;
tlag_eff1=tlag_eff;
azi4plot=azi-180;

save('fast_eff','fast_eff')
save('tlag_eff','tlag_eff')    
save('azi','azi')
save('downdipdir','downdipdir')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  ** make a figure of the results
%       scrsz = get(0,'ScreenSize');
%       figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;
% 
% %  ** lag times
%       subplot(2,1,2)
      %plot(azi,tlag1,'r--') ; hold on
      %plot(azi,tlag2,'g--') ;
      
      % recalculate azi to bazi for correct plotting over bazi between 0
      % and 360 degrees
      azi = [min_azi:del_azi:max_azi];
      
      %plot(azi,tlag_eff,'k-','LineWidth',1.5)
      %axis([min(azi) max(azi) 0 4])
      %xlabel('Polarisation (relative to downdip direction)')
      %ylabel('Lag times (s)')
      %legend('Upper layer','Lower layer','Total') ;
      %ylim([0 4])
      
%       getticks=get(gca,'xticklabel');
%       set(gca,'xticklabel',getticks+90)
      
      

%  ** fast directions
%      subplot(2,1,1)
      %plot(azi,fast1,'r--') ; hold on
      %plot(azi,fast2,'g--') ;
      
      % recalculate azi to bazi for correct plotting over bazi between 0
      % and 360 degrees
      azi = [min_azi:del_azi:max_azi];

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

      % now correclty sort between +-90°
      for ii=1:length(azi)
          if (fast_eff(ii)) > 90
              fast_eff(ii)=fast_eff(ii)-180;
          elseif (fast_eff(ii)) < -90
              fast_eff(ii)=fast_eff(ii)+180;
          else
               fast_eff(ii)=fast_eff(ii);
          end
      end
  
      %plot(azi,fast_eff,'k-','LineWidth',1.5)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output variables that are used to construct the OVERALL model matrix

% - fast_eff
% - tlag_eff   
% - azi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
      
      
%       for ii=1:length(azi)
%           if (fast_eff(ii)+downdipdir) > 90
%               plot(azi(ii),fast_eff(ii)+downdipdir-180,'k-','LineWidth',1.5)
%               hold on
%           elseif (fast_eff(ii)+downdipdir) < -90
%               plot(azi(ii),fast_eff(ii)+downdipdir+180,'k-','LineWidth',1.5)
%               hold on
%           else
%                plot(azi(ii),fast_eff(ii)+downdipdir,'k-','LineWidth',1.5)
%                hold on
%           end
%           hold on
%       end
      
      
%       axis([min(azi) max(azi) -90 90])
%       %set(gca,'xtick',[0:90:360])
%       %set(gca,'xticklabel',{'0','90','180','270','360'})
%       xlabel('Polarisation (relative to downdip direction)')
%       ylabel('Fast shear-wave orientation (degree)')      
%       legend('Upper layer','Lower layer','Total') ;
      
return

function [dist]=distance_in_dipping_layer(dip,aoi,thick,azi)

%  ** calculate apparent dip
      alp = atand(tand(dip) .* sind(azi-90)) ;
      
%  ** calculate distances
      gam = 90 - alp - aoi ;
      bet = 90 + alp ;
      dist = thick .* sind(bet) ./ sind(gam) ;

return

function report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
                      L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) 

% %  ** Banner
%       fprintf('\n** SKS Shear-wave splitting in two layer anisotropic models.\n\n')
% 
% %  ** Report layer details
%       fprintf('Layer 1: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
%          L1_depth, L1_thick, L1_dip)
%       fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
%          L1_aaz,L1_faln)
%       fprintf('\n')
% 
%       fprintf('Layer 2: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
%          L2_depth, L2_thick, L2_dip)
%       fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
%          L2_aaz,L2_faln)
% 
%       fprintf('\n** Plotting backazimuthal variation of splitting parameters ...\n')


return