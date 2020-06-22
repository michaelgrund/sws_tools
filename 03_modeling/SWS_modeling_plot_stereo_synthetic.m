function SWS_modeling_plot_stereo_synthetic(modsall_sort)
%
% This function generates stereoplots showing the splitting pattern
% of the best-fitting model



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

close all
clc

% normal settings
linew=3;
marks=8;
linewcirc=3;

usecmap=parula(181);

% sector plotting
colfill=[219,219,219]./256;
col_wedge=[255.9999999999999,255.9999999999999,255.9999999999999]./256;

% modelnumber to be plotted, 1 is best based on RMSE
usenum=14;


if strcmp(modsall_sort(usenum).mod_type,'dipping')


    tlag_eff=modsall_sort(usenum).dt4plot;
    fast_eff=modsall_sort(usenum).fast4plot;
    azi=modsall_sort(usenum).azi4plot;
    downdipdir=modsall_sort(usenum).downdipdir;
    dips=modsall_sort(usenum).dip;


%==========================================================

    stepsize=5;

    phi0=fast_eff(1:stepsize:360);
    dt0=tlag_eff(1:stepsize:360);
 
    phi0=phi0';
    dt0=dt0';
    
    len=dt0;  % for JOF scaling normal ==1

    bazi=azi(1:stepsize:180);
    
    for ii=1:length(bazi)
       
        if bazi(ii) > 360
            
            bazi(ii)=bazi(ii)-360;

        end
        
        
    end

    
elseif strcmp(modsall_sort(usenum).mod_type, 'two_layers')

 % index 1 nominates the lower layer,
% index 2 represents the upper layer   
    
    
    phi0=modsall_sort(usenum).phi_eff(1:5:360);
    dt0=modsall_sort(usenum).dt_eff(1:5:360);

    phi0=phi0';
    dt0=dt0';


    bazi=0:5:179;
  % [phi0, dt0]=twolayermodel(bazi, phi1,dt1, phi2, dt2, period);

    phi0 = [phi0; phi0];
    dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';
    
    
    bazi=bazi';
    len=dt0/2;  % for JOF scaling normal ==1


end

bazi = [bazi, bazi+180]';  
inc =ones(size(bazi))*10; %default 10deg inclination
%  inc2 =ones(size(bazi))*5; % 5deg inclination
azim=phi0;
azim_pre=phi0';
m = max(inc);
m = round(m/10)*10; %make gridline every 10deg
lim = [-inf m+5];

m1 = axesm ('stereo', 'Frame', 'on', 'Grid', 'on' ,'Origin',[90 0],...
   'MlineLocation', 90, 'PlineLocation', 5, 'fLatLimit',lim, 'fLineWidth',1, ...
   'GLinestyle','-', 'GLinewidth',0.4, 'Gcolor','k');

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mark null area in shaded gray

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% uses function <<< plot_arc >>> to plot wedge, modified
% function is <<< plot_arc3D >>> to also add a layer in 3rd dimension!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if modsall_sort(1).modrange_low < modsall_sort(1).modrange_upp

    if (modsall_sort(1).modrange_low ~=0 || modsall_sort(1).modrange_upp~=360)
    
        % first plot whole range in gray as bottom layer
        startwedge=0;
        endwedge=360;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,colfill);

        % then plot modelled range again on top in white ;)
        startwedge=modsall_sort(1).modrange_low;
        endwedge=modsall_sort(1).modrange_upp;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,col_wedge);

    end 

% this is the case when the higher value is slightly higher than 0 and the other one in the SA region
elseif modsall_sort(1).modrange_low > modsall_sort(1).modrange_upp 

    if (modsall_sort(1).modrange_low ~=0 || modsall_sort(1).modrange_upp~=360)
    
        % first plot whole range in gray as botom layer
        startwedge=0;
        endwedge=360;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,colfill);

        % then plot modelled range again on top in white ;) in two steps
        startwedge=0;
        endwedge=modsall_sort(1).modrange_upp;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,col_wedge);

        startwedge=modsall_sort(1).modrange_low;
        endwedge=360;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,col_wedge);

    end 
end

    
hold on
cmap=flipud(usecmap);
colormap(cmap);

NNull = 1:length(phi0);
bazi = bazi(:);

inc  = inc(:);
len  = len(:);
azim = azim(:);

bazi = [bazi(NNull)  bazi(NNull)]';
inc  = [inc(NNull)   inc(NNull)]';
len  = [-len(NNull)  len(NNull)]';
azim = (bazi-[azim(NNull) azim(NNull)]');
  
len=len*2; %one second == 4degrees (2 deg in both directions)
%len2=len2*2; %one second == 4degrees (2 deg in both directions)

% Marker
hold on

%%% add here value in which direction the axis should dip
if strcmp(modsall_sort(usenum).mod_type,'dipping')
    bazi=bazi+downdipdir;
end

[latout, lonout] = reckon(90-inc, bazi, len, azim, 'degrees');
hndl = plotm(latout, lonout, 'Linewidth',2);
hold off

%=======================================
% plot inc =10°
cmap=flipud(usecmap);   
step_phi=-90:1:90;
    
for ii=1:length(hndl) 
  
    if strcmp(modsall_sort(usenum).mod_type,'dipping')
  
        azim_rounded=floor(azim_pre(1,ii)+downdipdir); 
  
    else

        azim_rounded=floor(azim_pre(1,ii));

    end
    
    
    
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
         
         azim_rounded_all(ii)=azim_rounded;
 
         set(hndl(ii),'color',cmap(index,:))
         set(hndl(ii),'linewidth',linew+2)

end

hold on

%=============================================================================== 
% Plot potential nulls in direction of downdip-dir and perp to it  
% only valid for assuming the fast axis plunges with dip direction

if strcmp(modsall_sort(usenum).mod_type,'dipping')

    % first plot nulls in downdipdir +-180°
    bazi_nulls=[downdipdir downdipdir+180];
    inc_nulls=[10 10];   
    for KK=1:2
        plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks+5,'linewidth',linewcirc,'MarkerFaceColor','w');
    end

    % second find positions where bar orientation exact +-90° baz
    % for dipdir 90, 270° the nulls are not correctly shown
    if downdipdir < 90
       first=azim_rounded_all+90;
        firstmin=min(first);
 
        second=azim_rounded_all-90;
        secondmax=max(second);
    elseif downdipdir > 90 && downdipdir < 180
        first=azim_rounded_all+90;
        firstmin=max(first);
 
        second=azim_rounded_all-90;
        secondmax=min(second);
    elseif downdipdir > 270
        first=azim_rounded_all-90;
        firstmin=max(first);
 
        second=azim_rounded_all+90;
        secondmax=min(second);
    
    elseif downdipdir > 180 && downdipdir < 270

        first=azim_rounded_all(azim_rounded_all < 0 )-90;
        firstmin=max(first);
        second=azim_rounded_all(azim_rounded_all > 0 )-90;
        secondmax=min(second);

    else
        first=azim_rounded_all+90;
        firstmin=max(first);
 
        second=azim_rounded_all-90;
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
        plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks+5,'linewidth',linewcirc,'MarkerFaceColor','w');
    end
%===============================================================================

     
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make legend in lower right corner

hold on
col_leg='k';
fontsize_leg=12;
startval=0.25;

plot([startval startval],[.18 .18],'o','linewidth',linewcirc,'markersize',marks,'color',col_leg)
text(startval,.16,'Null', 'HorizontalAlignment', 'center','fontsize',fontsize_leg)

lengthbar=0.0710; % manually adjusted 
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
  
cbpos(1) = 0.42+cbpos(1);
cbpos(2) = 0.1+cbpos(2);
   
cbpos(3) = 0.4*cbpos(3);
cbpos(4) = 0.4*cbpos(4);
set(cb,'position',cbpos)
set(gca,'position',axpos)

%====================================================================

if strcmp(modsall_sort(usenum).mod_type, 'two_layers')
   
    %plot in edge the corresponding two layer model parameters

    lengthbar=0.07; % manually adjusted to fit the length to the bars plotted via plotm (see above)
    plot([-startval-0.013-lengthbar/2 -startval+lengthbar],[-.215 -.215],'-','linewidth',linew+15,'color',[0.9922    0.5529    0.2353])
    plot([-startval-0.013-lengthbar/2 -startval+lengthbar],[-.2465 -.2465],'-','linewidth',linew+15,'color',[0.4549    0.7686    0.4627])

    % index 1 nominates the lower layer,
    % index 2 represents the upper layer
    phi1=modsall_sort(usenum).phi(2);
    phi2=modsall_sort(usenum).phi(1);
    dt1=modsall_sort(usenum).dt(2);
    dt2=modsall_sort(usenum).dt(1);

    fontsize_leg=8;
    text(-startval-0.01-lengthbar/2,-.215,['\phi = ' num2str(phi1) '\circ, \deltat = ' num2str(dt1) ' s'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')
    text(-startval-0.01-lengthbar/2,-.2465,['\phi = ' num2str(phi2) '\circ, \deltat = ' num2str(dt2) ' s'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')


elseif strcmp(modsall_sort(usenum).mod_type, 'dipping')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %mark dip angle schematic in lower left corner

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %uses function <<< plot_arc >>> to plot wedge
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    colfill=[210 210 210]./256;%[189,189,189]./256;
    startwedge=90;
    endwedge=90+dips;

    xdirset=0.0;
    ydirset=-.26;

    plot_arc(deg2rad(startwedge-90),deg2rad(endwedge-90),-startval+xdirset-lengthbar/2,ydirset,0.12,colfill);

    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset ydirset],'-','linewidth',linew-1,'color','k')

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = dips; % 30°, falls im Bogenmaß unten cos, sin verwenden
    leng = 0.13;
    fin = beg + leng*[cosd(ang) sind(ang)];
  
    if dips==0
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],'-','linewidth',linew+35,'color',[255 	165 	079 ]./256)
    else
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],'-','linewidth',linew+20,'color',[255 	165 	079 ]./256 )
    end


    % plot small blue bars to represent aligned crystals
    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = dips; 
    scaleddist=1;
    leng = 0.12;
    fin = beg + leng*[cosd(ang) sind(ang)];

    if dips==0
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],'--','linewidth',linew-1,'color','b')
    else
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset+0.0*scaleddist fin(2)+0.0*scaleddist],'--','linewidth',linew-1,'color','b')
    end

    % small white and black lines for polishing
    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset-0.006 ydirset-0.006],'-','linewidth',linew+5,'color',[255.9999999999999 255.9999999999999 255.9999999999999]./256)
    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],[ydirset ydirset],'-','linewidth',linew-1,'color','k')

    % plot number of dip angle
    fontsize_leg=14;
    text(-startval+0.085-lengthbar/2,-.276,['\Psi = ' num2str(dips) '\circ'], 'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')

    % plot schematic raypath with 11° incidence
    xdirset=0.03;
    ydirset=-.26;

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = 79;
    leng = 0.165;
    fin = beg + leng*[cosd(ang) sind(ang)];

    plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],'-','linewidth',linew-1,'color','r')

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = 101; 
    leng = 0.165;
    fin = beg + leng*[cosd(ang) sind(ang)];

    plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],'-','linewidth',linew-1,'color','r')

    xdirset=-0.291;
    ydirset=-.271;
    plot(xdirset,ydirset,'^','markerfacecolor','w','markersize',14,'markeredgecolor','k','linewidth',2)
    
    %===============================================
    % PLOT dip dir as arrow

    % adopted from the drawarrow function available here:
    % https://de.mathworks.com/matlabcentral/fileexchange/55181-drawarrow
    
    
    p0 = [-0.07 0];
    p1 = [0.07,0];
    color_arrow = [120 120 120]./256;
 
    W1 = 0.14;   % half width of the arrow head, normalized by length of arrow
    W2 = 0.04;  % half width of the arrow shaft
    L1 = 0.38;   % Length of the arrow head, normalized by length of arrow
    L2 = 0.33;  % Length of the arrow inset

    x0 = p0(1);
    y0 = p0(2);
    x1 = p1(1);
    y1 = p1(2);

    P = [...
        0, (1-L2), (1-L1), 1, (1-L1), (1-L2), 0;
        W2,    W2,     W1, 0,    -W1,    -W2, -W2];

    dx = x1-x0;
    dy = y1-y0;
    Length = sqrt(dx*dx + dy*dy);
    Angle = atan2(-dy,dx);
    P = Length*P;  
    P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P; 
    P = p0(:)*ones(1,7) + P;  
    
    arrw = patch(P(1,:), P(2,:),color_arrow);  
    arrw.EdgeColor = 'k';
    set(arrw,'LineWidth',1.5)
    rotate(arrw,[0 0 1],downdipdir-90)  
  
end


axis tight

axis off
L = min(abs(axis));
str =strvcat([num2str(m/2)  char(186)], [num2str(m)  char(186)]);
% % t   = textm(90-[m/2; m], [45; 45], str, 'FontName','FixedWidth' );
 text(0 , -L-0.005, 'N','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20);
 text(L+0.005, 0,   'E','HorizontalAlignment','Left','verticalAlignment','middle','fontsize',20);

view([0 -90])





%===============================================================================
%===============================================================================

% save plots

if exist('RES_dir','var') % if given in varargin, first change to results folder where all plots should be stored
    cd(RES_dir)
end
  
colordisc=1;

% if colordisc==1 && savenamevar==1
%     filename=['PLOT_modelling_4GMT_DIPPING_' modsall_sort(1).staname '_THEO_stereo_dipdir' num2str(modsall_sort(1).downdipdir)...
%         '_dip' num2str(dips) '_thick' num2str(modsall_sort(1).thick)];
% elseif colordisc==1 && savenamevar==0
%         filename=['PLOT_modelling_4GMT_DIPPING_' modsall_sort(1).staname '_THEO_stereo_dipdir' num2str(modsall_sort(1).downdipdir)...
%         '_dip' num2str(dips) '_thick' num2str(modsall_sort(1).thick)];
% elseif colordisc==0 
%         filename=['PLOT_modelling_4GMT_DIPPING_' modsall_sort(1).staname '_THEO_stereo_dipdir' num2str(modsall_sort(1).downdipdir)...
%         '_dip' num2str(dips) '_thick' num2str(modsall_sort(1).thick)];
% end
    
    
    
    
    
%print ('-dpng', '-painters','-r600', [filename '.jpg']) 
%print ('-depsc', '-painters','-r600', [filename '.eps']) 
%dir_eps_file=dir([filename '.eps']);
%[status,cmdout]=system(['epstopdf ' dir_eps_file.name]);


%savefig('TEST1')


%plot2svg('TESTSVG',gcf,'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%exportgraphics(gcf,'barchart.png','Resolution',600)

saveas(gcf,'Barchart.png')
%print('BarPlot','-dpdf','-painters')









