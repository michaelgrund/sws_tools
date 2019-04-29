function SWS_Analysis_BASICS_stereoplot(colmap)
%===============================================================================
% Function to directly read single (SplitLab, SL) and multi-event (StackSplit, SS) 
% resultfiles, prepare stereoplots and save them as publication-ready pdf, png, jpg etc.
% Modified after the original < stereoplot > function of SplitLab.
%
% Main author: Michael Grund (michael.grund@kit.edu)
%
% created: 2018-06-11 -MG-
%===============================================================================
% How to:
%
% after processing the whole data of a station with SL (and SS)
%
% 1)  go to your results folder (folder where result lists of SL and SS measurments are saved)
%     required: splitresults_*.txt, splitresultsNULL_*.txt
%     optional: *_stackresults.mat
%
% 2) run this function SWS_Analysis_BASICS_stereoplot() 
%     => if result lists are available, they are loaded and processed completely automatically 
%
% If bars should be color-coded with respect to the fast axis, insert the
% colormap of your choice as input, e.g.:
%    
%     SWS_Analysis_BASICS_stereoplot('lajolla') 
%     SWS_Analysis_BASICS_stereoplot('viridis')
%     SWS_Analysis_BASICS_stereoplot('winter')
%
% If no argument is given, the default colormap parula (flipped) is used. 
%
%     SWS_Analysis_BASICS_stereoplot() or SWS_Analysis_BASICS_stereoplot
%
% For plotting without color-coding use
%     
%     SWS_Analysis_BASICS_stereoplot('none')
%
%===============================================================================
% If you downloaded the corresponding content the following colormaps are
% supported:
% 
% 1) standard MATLAB cmamps: parula, winter, summer, copper....
%
% 2) MatPlotLib 2.0 Colormaps: Perceptually Uniform and Beautiful 
%    (see: https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-2-0-colormaps-perceptually-uniform-and-beautiful)
%
% 3) Scientific colour-maps by F. Crameri (Zenodo. http://doi.org/10.5281/zenodo.1243862)
%    (see: http://www.fabiocrameri.ch/colourmaps.php)
%
% To plot colored backgrounds or limited sectors the function < plot_arc3D > is 
% attached which is based on < plot_arc > from Matt Fig, see also CHANGEABLE SETTINGS below
%
% Radial axis annotations (incidence angle 5, 10 and 15 degrees) can be plotted in 4 different quadrants: 
%     NE, SE, SW, NW
%===============================================================================
%===============================================================================

close all
clc

% CHANGEABLE SETTINGS 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% plot background (or not considered sector) in 
% specific color or white (function <<< plot_arc3D >>> required)
% default is white between 0 and 360 degrees
lowlim=0;     % adjust for your needs
upplim=360;   % adjust for your needs 

% Examples: 
% 1) only show sector between 20 and 120 degrees in white, rest is gray
%        lowlim=20;     
%        upplim=120;   
%        colfill=[190,190,190]./256;
% 2) shade full background in lighter gray
%        lowlim=0;     
%        upplim=0.001;   
%        colfill=[220,220,220]./256;

% color of radial annotation
annotcol='k';
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% output settings
outform='pdf'; % file format: pdf (default), png, jpg...
outres='600'; % resolution in dpi: 600 (default), 300, 100...

% plot settings
linew=4;
marks=9;
linewcirc=2.5;

% no phi color-coding
splitcol=[.3 .3 .3];
stackcol=[0 0.4470 0.7410];
nullcol=[0.8500 0.3250 0.0980];

% Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
cramericmap={'devon.mat','lajolla.mat','bamako.mat','davos.mat',...
        'bilbao.mat','nuuk.mat','oslo.mat','grayC.mat',...
        'hawaii.mat','lapaz.mat','tokyo.mat','buda.mat',...
        'acton.mat','turku.mat','imola.mat','broc.mat',...
        'cork.mat','vik.mat','lisbon.mat','tofino.mat',...
        'berlin.mat','batlow.mat','roma.mat','oleron.mat'};

% MatPlotLib 2.0 Colormaps: Perceptually Uniform and Beautiful  
mplcmap={'viridis','magma','inferno','plasma'};

% read SL single event results
[RES_split, RES_nulls]=SWS_Analysis_BASICS_read_SLresults();

%===============================================================================    

% define colormap
if ~exist('colmap','var')
    usecmap=parula(181); 
    usecmap=flipud(usecmap); % default
    colmap='parula_flip';
    fast_col=1;
elseif strcmp('none',colmap)
    fast_col=0;
else

    fast_col=1;
    % check for Scientific colour-maps of F. Crameri on your system
    checkcmaps=strfind(cramericmap,colmap);
    idx=find(not(cellfun('isempty',checkcmaps)), 1);

    % check for matplotlib colormaps on your system
    checkcmaps2=strfind(mplcmap,colmap);
    idx2=find(not(cellfun('isempty',checkcmaps2)), 1);

    if ~isempty(idx)
        loadmap=load([colmap '.mat']);
        getfname=fieldnames(loadmap);
        usecmap=loadmap.(getfname{1});
    end

    if ~isempty(idx2)
        usecmap=colormap([colmap '(181)']);
    end

    if isempty(idx) && isempty(idx2)
        if exist(colmap,'file') % check for other cmaps implemented in MATLAB (see MATLAB help: colormap) 
            usecmap=colormap([colmap '(181)']);
        else
            error('Colormap not available!')
        end
    end
end
%===============================================================================
% read stacking results if available and make query

plotmulti=0;

if ~exist('plot_stacks','var')

    dir_res_stack=dir('*_stackresults.mat');

    if ~isempty(dir_res_stack)
        RES_STACK=SWS_Analysis_BASICS_read_SSresults(dir_res_stack,1);
        
        disp(' ')
        plotmulti=input('Plot stacking_results (if available)? [1]=yes [0]=no');
    else
        RES_STACK=[];
    end

end
%===============================================================================

if ~isempty(RES_split)
    bazi=[RES_split.baz];
    inc=[RES_split.inc];
    azim_pre=[RES_split.phiSC];
    azim=azim_pre;
    len=[RES_split.dtSC];
    staname=RES_split.staname;
end

if ~isempty(RES_nulls)
    bazi_nulls=[RES_nulls.baz];
    inc_nulls=[RES_nulls.inc];
    azim_nulls_pre=[RES_nulls.phiSC];
    azim_nulls=azim_nulls_pre;
    len_nulls=[RES_nulls.dtSC];
    staname=RES_nulls.staname;
end

if ~isempty(RES_STACK)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate mean incidence for each stack
    for ii=1:length(RES_STACK)

        for jj=1:length(RES_STACK(ii).used_phases)
            incALL(jj)=RES_STACK(ii).used_phases(jj).results.incline;
        end
        
        incALLmean(ii)=mean(incALL);
        clear incALL
    end   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    bazi_stack=[RES_STACK.meanbaz];
    inc_stack=[incALLmean];
    azim_STACK_pre=[RES_STACK.phiSTACK];
    azim_stack=azim_STACK_pre;
    len_stack=[RES_STACK.dtSTACK];

end

%===============================================================================
% plot annotation on circles

disp(' ')
plotannot=input('Annotate radial scale ? [1]=NE [2]=SE [3]=SW [4]=NW [0]=no');

%===============================================================================

f1=figure;

if ~isempty(RES_split)
    m = max(inc);
elseif ~isempty(RES_nulls)
    m = max(inc_nulls);
end

m = round(m/10)*10; %make gridline every 10deg
lim = [-inf m+5];

axesm ('stereo', 'Frame', 'on', 'Grid', 'on' ,'Origin',[90 0],...
    'MlineLocation', 90, 'PlineLocation', 5, 'fLatLimit',lim, 'fLineWidth',1, 'GLinestyle','-', 'GLinewidth',0.4, 'Gcolor','k');


%===============================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mark null area in shaded gray

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% uses function <<< plot_arc >>> to plot wedge, modified
% function is <<< plot_arc3D >>> also add a layer in 3rd dimension!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if exist('plot_arc3D','file')

    if lowlim < upplim

        if (lowlim ~=0 || upplim ~=360)

            % first plot whole range in gray as bottom layer
            startwedge=0;
            endwedge=360;
            plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,colfill);

            % then plot modelled range again on top in white ;)
            startwedge=lowlim;
            endwedge=upplim;
            plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,[255.99999999 255.99999999 255.99999999]./256);
        end

    elseif lowlim > upplim % this is the case when the higher value is slightly higher than 0 and the other one in the SA region

        if (lowlim ~=0 || upplim~=360)

            % first plot whole range in gray as bottom layer
            startwedge=0;
            endwedge=360;
            plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,colfill);

            % then plot modelled range again on top in white ;) in two steps
            startwedge=0;
            endwedge=upplim;
            plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,[255.99999999 255.99999999 255.99999999]./256);

            startwedge=lowlim;
            endwedge=360;
            plot_arc3D(deg2rad(startwedge-90),deg2rad(endwedge-90),0,0,0.2635,[255.99999999 255.99999999 255.99999999]./256);
        end
    end
end

%===============================================================================
   
hold on

% Nulls
for KK=1:length(RES_nulls)
    if ~isempty(RES_nulls) && fast_col==0
        plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color',nullcol, 'MarkerSize',marks,'linewidth',linewcirc);
    elseif  ~isempty(RES_nulls) && fast_col==1
        cmap=usecmap;
        colormap(cmap);
        plotm(90-inc_nulls(KK), bazi_nulls(KK) ,'o','color','k', 'MarkerSize',marks,'linewidth',linewcirc,'markerFacecolor','w');
        cmap=usecmap;
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

if ~isempty(RES_STACK)

    NNull = 1:length(RES_STACK);% 

    bazi_stack = bazi_stack(:);
    inc_stack  = inc_stack(:);
    len_stack  = len_stack(:);
    azim_stack = azim_stack(:);

    bazi_stack = [bazi_stack(NNull)  bazi_stack(NNull)]';
    inc_stack  = [inc_stack(NNull)   inc_stack(NNull)]';
    len_stack  = [-len_stack(NNull)  len_stack(NNull)]';
    azim_stack = (bazi_stack-[azim_stack(NNull) azim_stack(NNull)]');

end

%===============================================================================

if ~isempty(RES_split)

%scale marker to output size
len=len*2; %one second == 4degrees (2 deg in both directions)

if ~isempty(RES_STACK)
    len_stack=len_stack*2;
end

% singles
[latout, lonout] = reckon( 90-inc, bazi, len, azim, 'degrees');
hndl = plotm(latout, lonout, 'Linewidth',1);

% stacks
if exist('RES_STACK','var') && plotmulti==1
    [latoutstack, lonoutstack] = reckon( 90-inc_stack, bazi_stack, len_stack, azim_stack, 'degrees');
    hndlstack = plotm(latoutstack, lonoutstack, 'Linewidth',1);
end

hold off

% display fast axis in color
if fast_col==1

    cmap=usecmap;  
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

    if exist('RES_STACK','var') && plotmulti==1
    
        linewstack=linew;
    
        for ii=1:length(hndlstack) 

            if ~isempty(RES_STACK)
                azim_rounded=floor(azim_STACK_pre(1,ii));
            end
            
            index=step_phi==azim_rounded;
            set(hndlstack(ii),'color',cmap(index,:))
            set(hndlstack(ii),'linewidth',linewstack)

        end

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

colormap(usecmap);

   cb = colorbar('location','north');
   zlab = get(cb,'xlabel');
   set(zlab,'String','\phi in \circ');
   set(zlab,'fontsize',10);

   caxis([-90 90])
   set(cb,'xtick',[-90:30:90]);
   set(cb,'fontsize',10)
   
%    axpos = get(gca,'position');
%    cbpos=get(cb,'position');
   
%    cbpos(1) = 0.38+cbpos(1);
%    cbpos(2) = 0.1+cbpos(2);
%    
%    cbpos(3) = 0.4*cbpos(3);
%    cbpos(4) = 0.4*cbpos(4);
%    set(cb,'position',cbpos)
%    set(gca,'position',axpos)
   
   % use predefined size values
   set(cb,'position',[0.6157    0.9482    0.2271    0.0238])
   set(gca,'position',[0.1300    0.1100    0.7750    0.8150])

else
    set(hndl,'color',splitcol,'linewidth',2.5)
    set(hndlstack,'color',stackcol,'linewidth',2.5)
end

else
    
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
   
%    axpos = get(gca,'position');
%    cbpos=get(cb,'position');
   
%    cbpos(1) = 0.38+cbpos(1);
%    cbpos(2) = 0.1+cbpos(2);
%    
%    cbpos(3) = 0.4*cbpos(3);
%    cbpos(4) = 0.4*cbpos(4);
%    
%    set(cb,'position',cbpos)
%    set(gca,'position',axpos)
   
   % use predefined size values
   set(cb,'position',[0.6157    0.9482    0.2271    0.0238])
   set(gca,'position',[0.1300    0.1100    0.7750    0.8150])

end

axis tight
axis off
L = min(abs(axis));
text(0 , -L-0.005, 'N','HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20);
text(L+0.005, 0,   'E','HorizontalAlignment','Left','verticalAlignment','middle','fontsize',20);

view([0 -90])
text(-0.2, -L-0.005, staname,'HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20,'backgroundcolor','w','edgecolor','k');

%====================================================================
% annotate radial axis?
hold on

if isempty(plotannot) % default
    plotannot=2;
end

if plotannot==1
   text(0.047,-0.047,'5\circ','fontsize',12,'color',annotcol) 
   text(0.105,-0.105,'10\circ','fontsize',12,'color',annotcol) 
   text(0.168,-0.168,'15\circ','fontsize',12,'color',annotcol)
end

if plotannot==2
   text(0.047,0.047,'5\circ','fontsize',12,'color',annotcol) 
   text(0.105,0.105,'10\circ','fontsize',12,'color',annotcol) 
   text(0.168,0.168,'15\circ','fontsize',12,'color',annotcol)
end

if plotannot==3
   text(-0.055,0.055,'5\circ','fontsize',12,'color',annotcol) 
   text(-0.117,0.117,'10\circ','fontsize',12,'color',annotcol) 
   text(-0.18,0.18,'15\circ','fontsize',12,'color',annotcol)
end

if plotannot==4
   text(-0.055,-0.055,'5\circ','fontsize',12,'color',annotcol) 
   text(-0.117,-0.117,'10\circ','fontsize',12,'color',annotcol) 
   text(-0.18,-0.18,'15\circ','fontsize',12,'color',annotcol)
end

%===============================================================================
%===============================================================================
% set papersize (crop white space at edges) and save plot

filename=['PLOT_RESULTS_stereo_' staname '_' colmap];


%........................................
% predefined and optimized to give nearly identical sizes on most systems
set(f1, 'PaperPosition', [-3.4655 -1.1548 20.3046 15.2284]);
set(f1,'PaperSize',[14 14]); %set paper size 

%........................................
% if saved figures do not fulfill your expectations you can adjust the size
% and scaling here by uncomment the following lines and modify the individual values:

% figsize=get(f1, 'PaperPosition');
% set(f1, 'PaperPosition', [figsize(1)-4.1 figsize(2)-7.5  figsize(3) figsize(4)]);

%........................................
print (['-d' outform], '-painters',['-r' outres], [filename '.' outform]) 

%........................................
% % alternatively on Linux systems print as eps and then convert to pdf
% 
% print ('-depsc', '-painters','-r600', [filename '.eps']) 
% dir_eps_file=dir([filename '.eps']);
% [status,cmdout]=system(['epstopdf ' dir_eps_file.name]);

end