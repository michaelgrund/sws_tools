function SWS_stereoplot_col(fast_col,RES_split,RES_nulls,savedir)
%
% This function is a modified version of SplitLabs < stereoplot.m >
% function (original author: A. Wüstefeld). It allows to plot the
% splitting parameters given via the structs RES_split and RES_nulls in 
% stereoplot view over backazimuth and incidence angle. Additionally, it
% allows to color-code the individual bars with respect to the measured 
% fast axis using the SC method.
%
% The function is primarily  modified to read the full data set published by
% Grund & Ritter (2020). Functionality with other data sets is not
% guaranteed. For the required mat-file structure see the data available
% from: 
%
%  https://publikationen.bibliothek.kit.edu/1000091427
%
% 2019-04-10 -MG-
% ORCID: https://orcid.org/0000-0001-8759-2018
% GitHub: https://github.com/michaelgrund/sws_tools
%
% see also function: SWS_read_evstruct
%===============================================================================

close all
clc

% basic settings
linew=4; % linewidth split bars
marks=9; % markersize nulls
linewcirc=2.5; % width of null circles
nocolorfast=[.3 .3 .3]; % color used if fast axis are not colored with respect to phi
nullcol=[0.8500 0.3250 0.0980];
usecmap=parula(181); % define colormap

%===============================================================================

if ~isempty(RES_split)
    bazi=[RES_split.BAZ];
    inc=[RES_split.incline];
 
    for ii=1:length(RES_split)
        azim_pre(ii)=RES_split(ii).phiSC(2); % only take measured value, errors are not included here
    end

    azim=azim_pre;
    len=[RES_split.dtSC];
    staname=RES_split.staname;
end

if ~isempty(RES_nulls)
    bazi_nulls=[RES_nulls.BAZ];
    inc_nulls=[RES_nulls.incline];
    
    for ii=1:length(RES_nulls)
        azim_nulls_pre(ii)=RES_nulls(ii).phiSC(2); % only take measured value, errors are not included here
    end
    
    azim_nulls=azim_nulls_pre;
    len_nulls=[RES_nulls.dtSC];
    staname=RES_nulls.staname;
end

%===============================================================================

if ~isempty(RES_split)
    m = max(inc);
elseif ~isempty(RES_nulls)
    m = max(inc_nulls);
else
    hndl=[];
    return
end

m = round(m/10)*10; %make gridline every 10 deg
lim = [-inf m+5];

axesm ('stereo', 'Frame', 'on', 'Grid', 'on' ,'Origin',[90 0],...
    'MlineLocation', 90, 'PlineLocation', 5, 'fLatLimit',lim, 'fLineWidth',1, 'GLinestyle','-', 'GLinewidth',0.4, 'Gcolor','k');

% PLOT Nulls
for KK=1:length(RES_nulls)

    if ~isempty(RES_nulls) && fast_col==0
        plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color',nullcol, 'MarkerSize',marks,'linewidth',linewcirc);
    elseif  ~isempty(RES_nulls) && fast_col==1
       cmap=flipud(usecmap);
       colormap(cmap);
       plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks,'linewidth',linewcirc,'markerFacecolor','w');
    else
       cmap=flipud(usecmap);
       colormap(cmap);
    end
end

%===============================================================================
 
if ~isempty(RES_split)
    NNull = 1:length(RES_split);
    bazi = bazi(:);
    inc  = inc(:);
    len  = len(:);
    azim = azim(:);
    bazi = [bazi(NNull)  bazi(NNull)]';
    inc  = [inc(NNull)   inc(NNull)]';
    len  = [-len(NNull)  len(NNull)]';
    azim = (bazi-[azim(NNull) azim(NNull)]');
else
    NNull = 1:length(RES_nulls);
    bazi = bazi_nulls(:);
    inc  = inc_nulls(:);
    len  = len_nulls(:);
    azim = azim_nulls(:);
    bazi = [bazi(NNull)  bazi(NNull)]';
    inc  = [inc(NNull)   inc(NNull)]';
    len  = [-len(NNull)  len(NNull)]';
    azim = (bazi-[azim(NNull) azim(NNull)]');   
end

% PLOT Splits
if ~isempty(RES_split)
    %scale marker to output size
    len=len*2; %one second == 4degrees (2 deg in both directions)

    hold on
    [latout, lonout] = reckon( 90-inc, bazi, len, azim, 'degrees');
    hndl = plotm(latout, lonout, 'Linewidth',1);

    % display fast axis in color with respect to phiSC
    if fast_col==1

    cmap=flipud(usecmap); 
    step_phi=-90:1:90;

    for ii=1:length(hndl) 
        
        if ~isempty(RES_split)
            azim_rounded=floor(azim_pre(1,ii));
        else
            azim_rounded=floor(azim_nulls_pre(1,ii)); 
        end
            
        index=step_phi==azim_rounded;
        set(hndl(ii),'color',cmap(index,:))
        set(hndl(ii),'linewidth',linew)
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%====================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make legend in lower right corner

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%====================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% make colorbar

colormap(flipud(usecmap));

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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%====================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    else % display fast axis only in a single color
        set(hndl,'color',nocolorfast,'linewidth',2.5)       
    end

else
    
    if fast_col==1
    
    % plot legend and colorbar also for NULL-stations
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
   
   end
   
end

axis tight
axis off
L=min(abs(axis));
text(0 , -L-0.005, 'N','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20);
text(L+0.005, 0,   'E','HorizontalAlignment','Left','verticalAlignment','middle','fontsize',20);

view([0 -90])
text(-0.2, -L-0.005, staname,'HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20,'backgroundcolor','w','edgecolor','k');

%===============================================================================

% save plots

if fast_col==1 
     filename=['PLOT_RESULTS_stereo_' staname '_fastcol'];
else 
     filename=['PLOT_RESULTS_stereo_' staname]; 
end

print ('-dpdf', '-painters','-r600', [savedir '/' filename '.pdf']) 
 
print ('-depsc', '-painters','-r600', [savedir '/' filename '_cut.eps']) 

% epstopdf installed?
[status,cmdout]=system('which epstopdf');

if status == 0
    dir_eps_file=dir([savedir '/' filename '_cut.eps']);
    [status,cmdout]=system(['epstopdf ' savedir '/' dir_eps_file.name]);
end

end

%===============================================================================
% EOF
