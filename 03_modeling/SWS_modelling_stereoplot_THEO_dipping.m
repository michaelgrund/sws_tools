function [hndl, marker] = SWS_modelling_stereoplot_THEO_dipping(dips)
%[hndl, marker] = stereoplot(bazi, inc, azim, len)
% plot a stereomap of values, at backazimuth bazi, with inclination inc.
% The direction and length (i.e. delay time) of the marker is determined 
% by azim and len, respectively.
% The optional argument Null is a vector of indices to the values with Nll
% charatesistics. these values are plotted as circles at the correspondig 
% backazimuth and inclination. 
%
% Example:
%  Imagine a station with 10 measurements. 
%  The third, fifth and nineth value are Null measurements:
%  [hndl, marker] = stereoplot(bazi, inc, azim, len, [3 5 9])

%===============================================================================
%===============================================================================

% modified to directly read resultfiles and prepare stereoplots afterwards
% and save them as pdf

% After processing the whole data of a station 
%
% 1)  go to your results folder (folder where result lists of SL measurments are saved)
% 2a) run this function without argmuent SWS_Analysis_BASICS_SLresults_stereoplot() 
%     => if a results list is available, it is loaded and processed completly automatically 
%
% or
%
% 2b) if you first want to have a look into the results struct: load a results file via 
%     function SWS_Analysis_BASICS_SLresults_read and use their ouput
%     arguments RES_split & RES_nulls as input for this 
%     function SWS_Analysis_BASICS_SLresults_stereoplot(RES_split, RES_nulls)

%===============================================================================
%%

% generate models for visualisation

%close all

hold on

steigung=-0;

x=0:100;
y=steigung.*x;
plot(x,y)
hold on

x=0:100;
y=steigung.*x-100;
plot(x,y)

% xlim([0 10])
% ylim([-1000 0])
%%
%close all

L=10;
alpha=200

%L is the length
%angle is alpha
x2=x+(L*cos(alpha));
y2=y+(L*sin(alpha));
plot([x x2],[y y2])






















%%


close all
clc

% normal settings
linew=3;
marks=8;
linewcirc=2;

% discrepant settings
discrepant_col=[204 0 0]./256;
discrepant_linew=4.5;
discrepant_marks=10;
discrepant_linewcirc=2.5;

usecmap=parula(181);

savenamevar=0;

%

% [70;35]	[0,500000000000000;0,750000000000000]
% [55;165]	[1,50000000000000;0,500000000000000]
% [60;20]	[1;0,500000000000000]
% [80;35]	[0,500000000000000;1]
% [130;45]	[1;2]
% [60;30]	[0,750000000000000;0,500000000000000]
% [130;45]	[1,25000000000000;2,25000000000000]
% [85;35]	[0,500000000000000;1]
% [55;160]	[1,25000000000000;0,500000000000000]
% [115;40]	[0,500000000000000;1,25000000000000]
% [105;40]	[0,500000000000000;1,50000000000000]
% [65;20]	[0,750000000000000;0,500000000000000]
% [60;5]	[1,25000000000000;0,500000000000000]
% [125;45]	[0,750000000000000;2]
% [60;25]	[1;0,500000000000000]
% [60;0]	[1,25000000000000;0,500000000000000]
% [60;180]	[1,25000000000000;0,500000000000000]
% [130;45]	[1,50000000000000;2,50000000000000]
% [55;160]	[1,50000000000000;0,500000000000000]
% [60;10]	[1,25000000000000;0,500000000000000]
% [65;30]	[0,500000000000000;0,500000000000000]



% index 1 nominates the lower layer,
% index 2 represents the upper layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best model for KAF based on RMS when only phi is fitted,
% note that dt is scaled down below, otherwise bars overlap with edges and
% are not shown!!!
% 06-2018 MG
    phi1 = 80;
    phi2 = 35;
    dt1  = 4;
    dt2  = 3.75;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best model for KEF based on RMS when only phi is fitted,
% note that dt is scaled down below, otherwise bars overlap with edges and
% are not shown!!!
% 06-2018 MG
    phi1 = 105;
    phi2 = 60;
    dt1  = 3.5;
    dt2  = 3.25;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best model for FIA1 based on RMS when only phi is fitted,
% note that dt is scaled down below, otherwise bars overlap with edges and
% are not shown!!!
% 06-2018 MG
    phi1 = 80;
    phi2 = 30;
    dt1  = 0.4;
    dt2  = 0.25;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use synthetic generated with split_model
% and prepare for plotting here

clear phi0 tlag0 azi bazi

load tlag_eff.mat
load fast_eff.mat
load azi.mat
load downdipdir.mat


% fast_eff=fast_eff+90;
% 
% for ii=1:length(fast_eff)
%     
%     if fast_eff(ii) > 90 
%         
%         fast_eff(ii)=fast_eff(ii)-180;
%         
%     end
% end





 phi0=fast_eff(1:5:360);
 dt0=tlag_eff(1:5:360);
% 
phi0=phi0';
dt0=dt0';



    
    
    
    






%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%     phi1 = BEST_models_sort(1).phi(1);
%     phi2 = BEST_models_sort(1).phi(2);
%     dt1  = BEST_models_sort(1).dt(1);
%     dt2  = BEST_models_sort(1).dt(2);


    period = 8;

    bazi=azi(1:5:180);
    %
    for ii=1:length(bazi)
       
        if bazi(ii) > 360
            
            bazi(ii)=bazi(ii)-360;

        end
        
        
    end
    
    
    %bazi=0:5:179;
    
   %[phi0, dt0]=twolayermodel(bazi, phi1,dt1, phi2, dt2, period);
 
     % if a 2 layer model is used here the pattern is duplicate since there
     % is a 90° periodicity, this is not valid when a single layer that is
     % dipping is used...
    %phi0 = [phi0; phi0];
    %dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';

    %inc =ones(size(bazi))*10; %default 10deg inclination
    inc =ones(size(bazi))*10; %default 10deg inclination
    
    azim=phi0;
    len=dt0/1;  % for JOF scaling normal ==1
 % 

%
%

 % find values perpendicular to incoming pol. direction    
    findnulls=floor(azim+90); 



%


%
% %===============================================================================
% 
% 
% if ~isempty(varargin) && nargin == 3 % use function after manually loading result files in folder 
%     
%     if nargin == 2
% 
%         RES_split=varargin{1};
%         RES_nulls=varargin{2};
% 
%     else
%         error('only two inputs are allowed!')
%     end
%     
% elseif ~isempty(varargin) && nargin == 2 % use function via mgrund_RESULTS_VIS_helper, 
%                                          % varargin have to be the destination path for all created plots
%     
%        RES_dir=varargin{1};
%        qual=2;
%        [RES_split, RES_nulls]=SWS_Analysis_BASICS_SLresults_read(qual);
%     
% elseif isempty(varargin)                % use function to automatically gen plot for one station in corresponding results
%                                         % folder 
%     [RES_split, RES_nulls]=SWS_Analysis_BASICS_SLresults_read();
% 
% end
% 
% 
% if ~exist('fast_col','var')
%     fast_col=1;
% end
% 
% 
% disp(' ')
%


%===============================================================================
% color discrepant pairs if available in color

%colordisc=input('Color discrepant SKS-SKKS pairs? [1]=yes [0]=no');
%colordisc=1;
%===============================================================================










%if ~isempty(RES_split)
%     bazi=[RES_split.baz];
%     inc=[RES_split.inc];
      azim_pre=phi0';
%      azim=azim_pre;
%     len=[RES_split.dtSC];
    
    staname='teststation';
%end

% if ~isempty(RES_nulls)
% %     bazi_nulls=[RES_nulls.baz];
% %     inc_nulls=[RES_nulls.inc];
% %     azim_nulls_pre=[RES_nulls.phiSC];
% %     azim_nulls=azim_nulls_pre;
% %     len_nulls=[RES_nulls.dtSC];
%     
%     staname=RES_nulls.staname;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find discrepant pairs and plot them in red

% for ii=1:length(RES_split)
% 
%     for jj=1:length(RES_nulls)
%         
%         
%         if strcmp(RES_split(ii).date_doy,RES_nulls(jj).date_doy)...
%                 && RES_split(ii).baz==RES_nulls(jj).baz...
%                 && (strcmp(RES_split(ii).quality_manual,'good') || strcmp(RES_split(ii).quality_manual,'fair'))...
%                 && (strcmp(RES_nulls(jj).quality_manual,'good') || strcmp(RES_nulls(jj).quality_manual,'fair'))
%         
%             RES_split(ii).phase_col='r';%[255 0 0]./256;
%             RES_nulls(jj).phase_col='r';%[255 0 0]./256;
%             
%         end
%         
%         
%         
%         
%     
%     end
% 
% end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%===============================================================================
% plot fast axis colored depending on fast axis angle

% disp(' ')
% fast_col=input('Plot fast axis colored?: [1] = yes [0] =no');
% 
% if isempty(fast_col)
%     fast_col=0;
% end

%===============================================================================
% if ~isempty(RES_split)
    m = max(inc);
% elseif ~isempty(RES_nulls)
%     m = max(inc_nulls);
% else
%     hndl=[];
%     marker=[];
%     return
% end











m = round(m/10)*10; %make gridline every 10deg
lim = [-inf m+5];

m1=axesm ('stereo', 'Frame', 'on', 'Grid', 'on' ,'Origin',[90 0],...
    'MlineLocation', 90, 'PlineLocation', 5, 'fLatLimit',lim, 'fLineWidth',1, 'GLinestyle','-', 'GLinewidth',0.4, 'Gcolor','k');
%'Gcolor',[100 100 100]./256

%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mark null area in shaded gray

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% uses function <<< plot_arc >>> to plot wedge
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

colfill=[189,189,189]./256;
% P=plot_arc(0+0.18,1.69*pi,0,0,0.2635,colfill);



startwedge=30
endwedge=80


%P=plot_arc(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,colfill);

set(m1, 'Layer', 'top')



hold on


%






% incnew=14.9;  
% plotm(90-incnew, 30 ,'o','color','g', 'MarkerSize',marks+10,'linewidth',linewcirc);%,'MarkerFaceColor','r');
% hold on  

% incnew=14.9;  
% plotm(90-incnew, 74 ,'o','color','g', 'MarkerSize',marks+10,'linewidth',linewcirc);%,'MarkerFaceColor','r');
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find 'r' in structs to map discrepant pairs

% indicesS=find(strcmp({RES_split.phase_col},'r'));
% indicesN=find(strcmp({RES_nulls.phase_col},'r'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Nulls
% 
% for KK=1:length(RES_nulls)
% 
% if ~isempty(RES_nulls) && fast_col==0
%    marker = plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','r', 'MarkerSize',marks,'linewidth',linewcirc);%,'MarkerFaceColor','r');
%        
%    cmap=flipud(usecmap);
%    colormap(cmap);
% 
% elseif  ~isempty(RES_nulls) && fast_col==1
       cmap=flipud(usecmap);
       colormap(cmap);
% 
%        
%             marker = plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks,'linewidth',linewcirc,'markerFacecolor','w');%,'MarkerFaceColor','r');
% else
%     marker=[];
%     
%     cmap=flipud(usecmap);
%     colormap(cmap);
% end
% 
% end



%



%
%
 
%if ~isempty(RES_split)


    NNull = 1:length(phi0);%     setdiff(1:length(inc), Null);%non-Nulls

    bazi = bazi(:);
    

    
    
    inc  = inc(:);
    len  = len(:);
    azim = azim(:);


    bazi = [bazi(NNull)  bazi(NNull)]';
    inc  = [inc(NNull)   inc(NNull)]';
    len  = [-len(NNull)  len(NNull)]';
    azim = (bazi-[azim(NNull) azim(NNull)]');
%
% else
%     
%     NNull = 1:length(phi0);%     setdiff(1:length(inc), Null);%non-Nulls
% 
%     bazi = bazi_nulls(:);
%     inc  = inc_nulls(:);
%     len  = len_nulls(:);
%     azim = azim_nulls(:);
% 
%     bazi = [bazi(NNull)  bazi(NNull)]';
%     inc  = [inc(NNull)   inc(NNull)]';
%     
%     len  = [-len(NNull)  len(NNull)]';
%     azim = (bazi-[azim(NNull) azim(NNull)]');   
% 
% end


%


%if ~isempty(RES_split)

%scale marker to output size
len=len*2; %one second == 4degrees (2 deg in both directions)


% set some markers manually to avoid overllapping with edge
%     len(:,10)=len(:,10)./1.3;
% len(:,46)=len(:,46)./1.3;

% Marker
hold on


% remove points from plot due to large delay times that would be plotted
% outside of frame
% len(:,2)=len(:,2)./1.1;
% len(:,24)=len(:,24)./1.1;

%
%%% add here value in which direction the axis should dip
bazi=bazi+downdipdir;
%azim=azim+downdipdir;



 [latout, lonout] = reckon(90-inc, bazi, len, azim, 'degrees');

hndl = plotm(latout, lonout, 'Linewidth',1);

hold off
%
% % display fast axis in color
% if fast_col==1
% 
     cmap=flipud(usecmap);
%         
     step_phi=-90:1:90;
% 
%     % phase specific color
%     
     for ii=1:length(hndl) 
%         
%         if ~isempty(RES_split)
                                          %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             azim_rounded=floor(azim_pre(1,ii)+downdipdir); % add here value of downdipdir to get correct color
%         else                                              %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             azim_rounded=floor(azim_nulls_pre(1,ii)); 
%         end
%   azim_rounded



         if azim_rounded < -90 && azim_rounded >=-270
             azim_rounded=azim_rounded+180;
             index=step_phi==azim_rounded;
         elseif azim_rounded > 90 && azim_rounded <=270
             azim_rounded=azim_rounded-180;
             index=step_phi==azim_rounded;
         elseif azim_rounded < -270
             azim_rounded=azim_rounded+360;
             index=step_phi==azim_rounded;
         elseif azim_rounded > 270
             azim_rounded=azim_rounded-360;
             index=step_phi==azim_rounded;
         else
            index=step_phi==azim_rounded;
         end
         
         azim_rounded_ALL(ii)=azim_rounded;




%          if azim_rounded < -90
%              azim_rounded=azim_rounded+180;
%              index=step_phi==azim_rounded;
%          elseif azim_rounded > 90
%              azim_rounded=azim_rounded-180;
%              index=step_phi==azim_rounded;
%          else
%             index=step_phi==azim_rounded;
%          end

         set(hndl(ii),'color',cmap(index,:))
         set(hndl(ii),'linewidth',linew)
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % make discrepant red and thicker
%         if ismember(ii,indicesS) && colordisc==1
%             set(hndl(ii),'color',discrepant_col) 
%             set(hndl(ii),'linewidth',discrepant_linew)
%         end
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     end
% 
     end
     
     
     
     
     
   hold on  

   
%=============================================================================== 
% Plot potential nulls in direction of downdip-dir and perp to it  
% only valid for assuming the fast axis plunges with dip direction

% first plot nulls in downdipdir +-180°
bazi_nulls=[downdipdir downdipdir+180];
inc_nulls=[10 10];   
for KK=1:2
   marker = plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');

end

% second find positions where bar orientation exact +-90° baz
% for dipdir 90, 270° the nulls are not correctly shown
if downdipdir < 90
    first=azim_rounded_ALL+90;
    firstmin=min(first);
 
    second=azim_rounded_ALL-90;
    secondmax=max(second);
elseif downdipdir > 90 && downdipdir < 180
    first=azim_rounded_ALL+90;
    firstmin=max(first);
 
    second=azim_rounded_ALL-90;
    secondmax=min(second);
elseif downdipdir > 270
    first=azim_rounded_ALL-90;
    firstmin=max(first);
 
    second=azim_rounded_ALL+90;
    secondmax=min(second);
    
else
    first=azim_rounded_ALL+90;
    firstmin=max(first);
 
    second=azim_rounded_ALL-90;
    secondmax=min(second);
end


if abs(secondmax-firstmin) < 5
    firstmin=firstmin+180;
end

if abs(secondmax-firstmin) < 5 
    secondmax=secondmax-180;
end

   
bazi_nulls=[firstmin secondmax];
inc_nulls=[10 10];   
for KK=1:2
   marker = plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');

end
%===============================================================================

     
     
     
     
     
%     %%%%%%%%%%%%%%% plot on top potential discrepant NULLs
%     hold on 
%     
%     
%     savenamevar=0;
%     for KK=1:length(RES_nulls)
% 
%         if ~isempty(RES_nulls) && fast_col==1 && strcmp(RES_nulls(KK).phase_col,'r') && colordisc==1
%        
%                 marker = plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color',discrepant_col, 'MarkerSize',discrepant_marks,'linewidth',discrepant_linewcirc,'markerFacecolor','w');%,'MarkerFaceColor','r');
%             savenamevar=1;
%         end
%     
%     end 
    %%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%====================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make legend in lower right corner

%%% ####### 
%%% 1 sec reference
% len=[-1 1]';
% len=len*2; %one second == 4degrees (2 deg in both directions)
% azim=[90 90]';
% bazi=[135 135]';
% inc=[13.05 13.05]';
% azim = (bazi-[azim(1) azim(1)]');
% 
% [latout, lonout] = reckon( 90-inc, bazi, len, azim, 'degrees');
% plotm(latout, lonout, 'Linewidth',12,'color','m');
%  
%%% ####### 
%%% 2 sec reference
% len=[-2 2]';
% len=len*2; %one second == 4degrees (2 deg in both directions)
% azim=[90 90]';
% bazi=[135 135]';
% inc=[3.05 3.05]';
% azim = (bazi-[azim(1) azim(1)]');
% 
%  [latout, lonout] = reckon( 90-inc, bazi, len, azim, 'degrees');
% plotm(latout, lonout, 'Linewidth',12,'color','m');


hold on
col_leg='k';
fontsize_leg=12;

startval=0.25;
plot([startval startval],[.18 .18],'o','linewidth',linewcirc,'markersize',marks,'color',col_leg)
text(startval,.16,'Null', 'HorizontalAlignment', 'center','fontsize',fontsize_leg)

lengthbar=0.0710; % manually adjusted to fit the length to the bars plotted via plotm (see above)
plot([startval-lengthbar/2 startval+lengthbar/2],[.22 .22],'-','linewidth',linew,'color',col_leg)
text(startval,.205,'1 s', 'HorizontalAlignment', 'center','fontsize',fontsize_leg)

lengthbar=2*0.0710;
plot([startval-lengthbar/2 startval+lengthbar/2],[.26 .26],'-','linewidth',linew,'color',col_leg)  
text(startval,.245,'2 s', 'HorizontalAlignment', 'center','fontsize',fontsize_leg)


   cb = colorbar('location','north');
   zlab = get(cb,'xlabel');
   set(zlab,'String','\phi in \circ');
   set(zlab,'fontsize',10);

   caxis([-90 90])
   set(cb,'xtick',[-90:30:90]);
   set(cb,'fontsize',10)
   
   axpos = get(gca,'position');
   cbpos=get(cb,'position');
   
   cbpos(1) = 0.38+cbpos(1);
   cbpos(2) = 0.1+cbpos(2);
   
   cbpos(3) = 0.4*cbpos(3);
   cbpos(4) = 0.4*cbpos(4);
   set(cb,'position',cbpos)
   set(gca,'position',axpos)





%====================================================================
% plot in edge the corresponding two layer model parameters


% lengthbar=0.0790; % manually adjusted to fit the length to the bars plotted via plotm (see above)
% plot([-startval+0.023-lengthbar/2 -startval+lengthbar],[.22 .22],'-','linewidth',linew+15,'color',[0.9922    0.5529    0.2353])
% plot([-startval+0.023-lengthbar/2 -startval+lengthbar],[.2465 .2465],'-','linewidth',linew+15,'color',[0.4549    0.7686    0.4627])
% 
% % index 1 nominates the lower layer,
% % index 2 represents the upper layer
% fontsize_leg=8;
% text(-startval+0.027-lengthbar/2,.22,['\phi = ' num2str(phi2) '\circ, \deltat = ' num2str(dt2) ' s'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')
% text(-startval+0.027-lengthbar/2,.2465,['\phi = ' num2str(phi1) '\circ, \deltat = ' num2str(dt1) ' s'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mark dip angle schematic in lower left corner

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% uses function <<< plot_arc >>> to plot wedge
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


colfill=[210 210 210]./256;%[189,189,189]./256;
startwedge=90;
endwedge=90+dips;

xdirset=0.0;
ydirset=-.26;

P=plot_arc(deg2rad(startwedge-90),deg2rad(endwedge-90),-startval+xdirset-lengthbar/2,ydirset,0.12,colfill);

plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset ydirset],'-','linewidth',linew-1,'color','k')



Anfang = [-startval+xdirset-lengthbar/2 ydirset];
Winkel = dips; % 30°, falls im Bogenmaß unten cos, sin verwenden
Laenge = 0.13;
Ende = Anfang + Laenge*[cosd(Winkel) sind(Winkel)];
  
if dips==0
      plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset Ende(2)],'-','linewidth',linew+35,'color',[255 	165 	079 ]./256)
else
    plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset Ende(2)],'-','linewidth',linew+20,'color',[255 	165 	079 ]./256 )
end
% small white line above to kaschieren ;)
plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset-0.005 ydirset-0.005],'-','linewidth',linew+5,'color',[255.9999999999999 255.9999999999999 255.9999999999999]./256)

plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset ydirset],'-','linewidth',linew-1,'color','k')

% plot number of dip angle
fontsize_leg=14;
text(-startval+0.085-lengthbar/2,-.276,['\Psi = ' num2str(dips) '\circ'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')


% line 2 (center line)
scaleddist=1;
Laenge = 0.12;
Ende = Anfang + Laenge*[cosd(Winkel) sind(Winkel)];

if dips==0
      plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset+0.015 Ende(2)+0.015],'--','linewidth',linew-1,'color','b')
else
      plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset+0.0*scaleddist Ende(2)+0.0*scaleddist],'--','linewidth',linew-1,'color','b')
end




% plot schematic raypath with 11° incidence
xdirset=0.03
ydirset=-.26

Anfang = [-startval+xdirset-lengthbar/2 ydirset];
Winkel = 79; % 30°, falls im Bogenmaß unten cos, sin verwenden
Laenge = 0.165;
Ende = Anfang + Laenge*[cosd(Winkel) sind(Winkel)]


plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset Ende(2)],'-','linewidth',linew-1,'color','r')

Anfang = [-startval+xdirset-lengthbar/2 ydirset];
Winkel = 101; % 30°, falls im Bogenmaß unten cos, sin verwenden
Laenge = 0.165;
Ende = Anfang + Laenge*[cosd(Winkel) sind(Winkel)]

plot([-startval+xdirset-lengthbar/2 Ende(1)],[ydirset Ende(2)],'-','linewidth',linew-1,'color','r')

xdirset=-0.291;
ydirset=-.269;
plot(xdirset,ydirset,'^','markerfacecolor','w','markersize',14,'markeredgecolor','k','linewidth',2)


if dips~=0
%%%% PLOT dip dir as arrow
% arrow to show relation for which the distance was calculated

%===============================================
%%%% PLOT dip dir as arrow

color_arrow=[120 120 120]./256;  
dav1=davinci( 'arrow', 'X',             [-0.07 0.07], ...
                  'Y',             [0 0], ...
                  'ArrowType',   'single', ...
                  'Shaft.Width',       .01, ...
                  'Head.Length',      0.04, ...
                  'Head.Sweep',         0.02, ...
                  'Head.Width',         0.04, ...
                  'Color',            color_arrow, ...
                  'EdgeColor',        'k', ...
                  'LineWidth',        1 );
              
rotate(dav1,[0 0 1],downdipdir-90)  

            
%text(750, 0.83, '20 s','fontsize',14,'color',color_arrow)

end







%%% lower left corner
% colfill=[210 210 210]./256%[189,189,189]./256;
% startwedge=90
% endwedge=90+60
% 
% P=plot_arc(deg2rad(startwedge-90),deg2rad(endwedge-90),-startval+0.03-lengthbar/2,.2,0.08,colfill);
% 
% plot([-startval+0.03-lengthbar/2 -startval+lengthbar/3],[.2 .2],'-','linewidth',linew-1,'color','k')
% 
% 
% Anfang = [-startval+0.03-lengthbar/2 .2];
% Winkel = 60; % 30°, falls im Bogenmaß unten cos, sin verwenden
% Laenge = 0.09;
% Ende = Anfang + Laenge*[cosd(Winkel) sind(Winkel)]
%   
% plot([-startval+0.03-lengthbar/2 Ende(1)],[.2 Ende(2)],'-','linewidth',linew-1,'color','k')
% 
% fontsize_leg=14;
% text(-startval+0.06-lengthbar/2,.225,[num2str(60) '\circ'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')
% 
% 
% 








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%====================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   





   
   
   
   
   
   
   
   
% else
%     set(hndl,'color',[ 0    0.4470    0.7410],'linewidth',2.5)       
% end
% 
% else
%     
%     hndl=[];
% 
% end











% phase specific color

% for ii=1:length(hndl)       
%     set(hndl(ii),'color',RES_split(ii).phase_col)
%     set(hndl(ii),'linewidth',2)
% end


axis tight

axis off
L = min(abs(axis));
str =strvcat([num2str(m/2)  char(186)], [num2str(m)  char(186)]);
% t   = textm(90-[m/2; m], [45; 45], str, 'FontName','FixedWidth' );
text(0 , -L-0.005, 'N','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20);
text(L+0.005, 0,   'E','HorizontalAlignment','Left','verticalAlignment','middle','fontsize',20);

view([0 -90])


%text(-0.2, -L-0.005, 'AAA','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20,'backgroundcolor','w','edgecolor','k','color','w');
%text(-0.2, -L+0.015, 'synthetic','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',14,'backgroundcolor','w','edgecolor','w');





%===============================================================================
%===============================================================================

% save plots

if exist('RES_dir','var') % if given in varargin, first change to results folder where all plots should be stored
    cd(RES_dir)
end
  
colordisc=1

if colordisc==1 && savenamevar==1
    filename=['PLOT_RESULTS_stereo_THEO_DISC_dipping_' num2str(dips) '_' staname];
elseif colordisc==1 && savenamevar==0
    filename=['PLOT_RESULTS_stereo_THEO_dipping_' num2str(dips) '_' staname];
elseif colordisc==0 
    filename=['PLOT_RESULTS_stereo_THEO_dipping_' num2str(dips) '_' staname];
end
    
    
    
    
    
%print ('-dpng', '-painters','-r600', [filename '.jpg']) 
print ('-depsc', '-painters','-r600', [filename '.eps']) 
dir_eps_file=dir([filename '.eps']);
[status,cmdout]=system(['epstopdf ' dir_eps_file.name]);
%[status,cmdout]=system(['pdfcrop ' filename '.pdf ' filename '.pdf']);

%plot2svg('TESTSVG',gcf,'png')






