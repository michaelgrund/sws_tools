function SWS_modeling_plot_stereo_synthetic(modsall_sort,plotnum)
%
% generate stereoplots displaying the (synthetic) splitting pattern
% of the best-fitting model, partly is based on SplitLab function
% stereoplot.m
%
% INPUT:
% modsall_sort: struct with models sorted based on RMSE,
%               see function SWS_modeling_calc_misfit(_GR)
% plotnum: model(s) to be plotted, 1 is best based on RMSE
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

figure()

% plotting settings
linew=3;
marks=8;
linewcirc=2;
usecmap=parula(181);

% sector plotting
colfill=[219,219,219]./256;
col_wedge=[255.9999999999999,255.9999999999999,255.9999999999999]./256;

%==========================================================

if strcmp(modsall_sort(plotnum).mod_type,'dipping')
    
    tlag_eff=modsall_sort(plotnum).dt4plot;
    fast_eff=modsall_sort(plotnum).fast4plot;
    azi=modsall_sort(plotnum).azi4plot;
    downdipdir=modsall_sort(plotnum).downdipdir;
    dips=modsall_sort(plotnum).dip;
    
    stepsize=5;
    phi0=fast_eff(1:stepsize:360);
    dt0=tlag_eff(1:stepsize:360);
    
    phi0=phi0';
    dt0=dt0';
    len=dt0/1; % scale bar lengths if they e.g. range over frame
    bazi=azi(1:stepsize:180);
    
    for ii=1:length(bazi)
        if bazi(ii) > 360
            bazi(ii)=bazi(ii)-360;
        end
    end

elseif strcmp(modsall_sort(plotnum).mod_type, 'two_layers')

    phi0=modsall_sort(plotnum).phi_eff(1:5:360);
    dt0=modsall_sort(plotnum).dt_eff(1:5:360);

    phi0=phi0';
    dt0=dt0';
    bazi=0:5:179;

    phi0 = [phi0; phi0];
    dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';

    bazi=bazi';
    len=dt0/1; % scale bar lengths if they e.g. range over frame

elseif strcmp(modsall_sort(plotnum).mod_type,'single_layer')
    
    phi0=modsall_sort(plotnum).phi_eff(1:5:360);
    dt0=modsall_sort(plotnum).dt_eff(1:5:360);

    phi0=phi0';
    dt0=dt0';
    bazi=0:5:179;

    phi0 = [phi0; phi0];
    dt0  = [dt0; dt0];
    bazi = [bazi, bazi+180]';

    bazi=bazi';
    len=dt0/1; % scale bar lengths if they e.g. range over frame
    
end

%==========================================================

bazi = [bazi, bazi+180]';  
inc = ones(size(bazi))*10; %default 10deg inclination
% inc2 = ones(size(bazi))*5; % 5deg inclination
azim=phi0;
azim_pre=phi0';
m = max(inc);
m = round(m/10)*10; %make gridline every 10deg
lim = [-inf m+5];

axesm ('stereo', 'Frame', 'on', 'Grid', 'on' ,'Origin',[90 0],...
   'MlineLocation', 90, 'PlineLocation', 5, 'fLatLimit',lim, 'fLineWidth',1, ...
   'GLinestyle','-', 'GLinewidth',0.4, 'Gcolor','k');

hold on

%==========================================================
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
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),...
            0,0,0.2635,colfill);

        % then plot modelled range again on top in white ;)
        startwedge=modsall_sort(1).modrange_low;
        endwedge=modsall_sort(1).modrange_upp;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),...
            0,0,0.2635,col_wedge);

    end 

% this is the case when the higher value is slightly higher than 0 and the other one in the SA region
elseif modsall_sort(1).modrange_low > modsall_sort(1).modrange_upp 

    if (modsall_sort(1).modrange_low ~=0 || modsall_sort(1).modrange_upp~=360)
    
        % first plot whole range in gray as botom layer
        startwedge=0;
        endwedge=360;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),...
            0,0,0.2635,colfill);

        % then plot modelled range again on top in white ;) in two steps
        startwedge=0;
        endwedge=modsall_sort(1).modrange_upp;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),...
            0,0,0.2635,col_wedge);

        startwedge=modsall_sort(1).modrange_low;
        endwedge=360;
        plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),...
            0,0,0.2635,col_wedge);

    end 
end

%==========================================================
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

if strcmp(modsall_sort(plotnum).mod_type,'dipping')
    bazi=bazi+downdipdir;
end

% plot the splits as bars
[latout, lonout] = reckon(90-inc, bazi, len, azim, 'degrees');
hndl = plotm(latout, lonout,-0.1);
%==========================================================

step_phi=-90:1:90;
    
for ii=1:length(hndl) 
  
    if strcmp(modsall_sort(plotnum).mod_type,'dipping')
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
         set(hndl(ii),'linewidth',linew)
end

hold on

%==========================================================
% plot potential nulls in direction of downdip-dir and perp to it  
% only valid for assuming the fast axis plunges with dip direction

if strcmp(modsall_sort(plotnum).mod_type, 'dipping')

    % first plot nulls in downdipdir +-180°
    bazi_nulls=[downdipdir downdipdir+180];
    inc_nulls=[10 10];   
    for k=1:2
        plotm(90-inc_nulls(k), bazi_nulls(k), -0.15,'o','color','k', ...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
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
    
    for k=1:2
        plotm(90-inc_nulls(k), bazi_nulls(k), -0.15,'o','color','k',...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
    end  
    
elseif strcmp(modsall_sort(plotnum).mod_type, 'single_layer')
    modsall_sort(plotnum).phi
    
    incs=10; 

    plotm(90-incs, modsall_sort(plotnum).phi, -0.15,'o','color','k',...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
    plotm(90-incs, modsall_sort(plotnum).phi+90, -0.15,'o','color','k',...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
    plotm(90-incs, modsall_sort(plotnum).phi+180, -0.15,'o','color','k',...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
    plotm(90-incs, modsall_sort(plotnum).phi+270, -0.15,'o','color','k',...
            'MarkerSize',marks,'linewidth',linewcirc,'MarkerFaceColor','w');
 
end   

%==========================================================
% make legend in lower right corner, always shown

hold on
col_leg='k';
fontsize_leg=10;
startval=0.25;

plot([startval startval],[.18 .18],'o','linewidth',linewcirc,...
    'markersize',marks,'color',col_leg)
text(startval,.1575,'Null', 'HorizontalAlignment',...
    'center','fontsize',fontsize_leg)

lengthbar=0.0710; % manually adjusted 
plot([startval-lengthbar/2 startval+lengthbar/2],[.22 .22],'-',...
    'linewidth',linew,'color',col_leg)
text(startval,.205,'1 s', 'HorizontalAlignment', 'center',...
    'fontsize',fontsize_leg)

lengthbar=2*0.0710;
plot([startval-lengthbar/2 startval+lengthbar/2],[.26 .26],'-',...
    'linewidth',linew,'color',col_leg)  
text(startval,.245,'2 s', 'HorizontalAlignment', 'center',...
    'fontsize',fontsize_leg)

%==========================================================
% colorbar in upper right corner, always shown

cb = colorbar('location','north');
zlab = get(cb,'xlabel');
set(zlab,'String','\phi in \circ');
set(zlab,'fontsize',10);

caxis([-90 90])
set(cb,'xtick',-90:30:90);
set(cb,'fontsize',10)
   
axpos = get(gca,'position');
cbpos=get(cb,'position');
  
cbpos(1) = 0.42+cbpos(1);
cbpos(2) = 0.1+cbpos(2);
   
cbpos(3) = 0.4*cbpos(3);
cbpos(4) = 0.4*cbpos(4);
set(cb,'position',cbpos)
set(gca,'position',axpos)

%==========================================================
% display station name in lower left corner, always shown

text(-startval,.205,modsall_sort(1).staname, 'HorizontalAlignment', ...
    'center','fontsize',fontsize_leg+8)

%==========================================================

if strcmp(modsall_sort(plotnum).mod_type, 'two_layers')
   
    %plot the corresponding two layer model parameters in upper right
    %corner
    lengthbar=0.075; 
    plot([-startval-0.016-lengthbar/2 -startval+lengthbar],[-.2215 -.2215],...
        '-','linewidth',linew+15,'color',[0.8 0.8 0.8])
    plot([-startval-0.016-lengthbar/2 -startval+lengthbar],[-.2588 -.2588],...
        '-','linewidth',linew+15,'color',[0.9922    0.5529    0.2353])

    % index 1 nominates the lower layer,
    % index 2 represents the upper layer
    phi1=modsall_sort(plotnum).phi(2);
    phi2=modsall_sort(plotnum).phi(1);
    dt1=modsall_sort(plotnum).dt(2);
    dt2=modsall_sort(plotnum).dt(1);

    fontsize_leg=8;
    text(-startval-0.015-lengthbar/2,-.2215,...
        ['\phi = ' num2str(phi1) '\circ, \deltat = ' num2str(dt1) ' s'],...
        'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')
    text(-startval-0.015-lengthbar/2,-.2588,...
        ['\phi = ' num2str(phi2) '\circ, \deltat = ' num2str(dt2) ' s'],...
        'HorizontalAlignment', 'left','fontsize',fontsize_leg,'color','k')

elseif strcmp(modsall_sort(plotnum).mod_type, 'dipping')

    %mark dip angle schematic in lower left corner

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %uses function <<< plot_arc >>> to plot wedge
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    colfill=[210 210 210]./256;
    startwedge=90;
    endwedge=90+dips;

    xdirset=0.0;
    ydirset=-.26;

    plot_arc(deg2rad(startwedge-90),deg2rad(endwedge-90),...
        -startval+xdirset-lengthbar/2,ydirset,0.12,colfill);
    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],...
        [ydirset ydirset],'-','linewidth',linew-1,'color','k')

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = dips; 
    leng = 0.13;
    fin = beg + leng*[cosd(ang) sind(ang)];
  
    if dips==0
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],...
            '-','linewidth',linew+35,'color',[255 	165 	079 ]./256)
    else
        plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],...
            '-','linewidth',linew+20,'color',[255 	165 	079 ]./256 )
    end
    
    % plot small blue bars to represent aligned crystals
    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = dips; 
    scaleddist=1;
    leng = 0.12;
    fin = beg + leng*[cosd(ang) sind(ang)];

    if dips==0
        plot([-startval+xdirset-lengthbar/2 fin(1)],...
            [ydirset fin(2)],'--','linewidth',linew-1,'color','b')
    else
        plot([-startval+xdirset-lengthbar/2 fin(1)],...
            [ydirset+0.0*scaleddist fin(2)+0.0*scaleddist],'--',...
            'linewidth',linew-1,'color','b')
    end

    % add small white and black lines for polishing
    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],...
        [ydirset-0.006 ydirset-0.006],'-','linewidth',linew+5,'color',...
        [255.9999999999999 255.9999999999999 255.9999999999999]./256)
    plot([-startval+xdirset-lengthbar/2 -startval+lengthbar/1.65],...
        [ydirset ydirset],'-','linewidth',linew-1,'color','k')

    % plot number of dip angle
    fontsize_leg=10;
    text(-startval+0.085-lengthbar/2,-.276,...
        ['\Psi = ' num2str(dips) '\circ'], 'HorizontalAlignment', ...
        'left','fontsize',fontsize_leg,'color','k')

    % plot schematic raypath with ~ 10° incidence
    xdirset=0.03;
    ydirset=-.26;

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = 79;
    leng = 0.165;
    fin = beg + leng*[cosd(ang) sind(ang)];

    plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],...
        '-','linewidth',linew-1,'color','r')

    beg = [-startval+xdirset-lengthbar/2 ydirset];
    ang = 101; 
    leng = 0.165;
    fin = beg + leng*[cosd(ang) sind(ang)];

    plot([-startval+xdirset-lengthbar/2 fin(1)],[ydirset fin(2)],...
        '-','linewidth',linew-1,'color','r')

    xdirset=-0.291;
    ydirset=-.271;
    plot(xdirset,ydirset,'^','markerfacecolor','w','markersize',14,...
        'markeredgecolor','k','linewidth',2)
    
    %===============================================
    % plot dip dir as arrow

    % adopted from the drawarrow function available here:
    % https://de.mathworks.com/matlabcentral/fileexchange/55181-drawarrow
    %
    % Matthew Kelly (2020). drawArrow (https://www.mathworks.com/matlabcentral/fileexchange/55181-drawarrow), 
    % MATLAB Central File Exchange. Retrieved June 22, 2020. 

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

%==========================================================

axis tight
axis off

lval = min(abs(axis));
text(0 , -lval - 0.005, 'N','HorizontalAlignment','Center',...
    'VerticalAlignment','Base','fontsize',20);
text(lval + 0.005, 0,   'E','HorizontalAlignment','Left',...
    'verticalAlignment','middle','fontsize',20);

view([0 -90])

%==========================================================
% save plots

filename=['PLOT_modeling_' modsall_sort(1).staname '_stereoplot_syn_' num2str(plotnum)];

vers_out=SWS_modeling_check_matlab_version();

if vers_out == 1
    saveas(gcf,[filename '.png'])
    %%% uncomment if you want plots in pdf format, takes time...
    % print([filename '.pdf'],'-dpdf','-painters') 
else
    exportgraphics(gcf,[filename '.png'],'Resolution',600)
    %%% uncomment if you want plots in pdf format, takes time...
    %print([filename '.pdf'],'-dpdf','-painters') 
end

%==========================================================
