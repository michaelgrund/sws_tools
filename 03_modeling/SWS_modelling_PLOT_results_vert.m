function SWS_modelling_PLOT_results_vert(BAZ,BEST_models_sort,plot_mod_max,...
    meas_BAZ_floor,meas_phiSC,meas_dtSC,meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,color_SS_bf_1,color_SS_bf_2max,linewidth_SS,colorsfill,colorsedge,...
    mymarkersize,mymarkersize_null,linewidth_symbols,linewidth_symbols_null,color_face_null,...
    color_edge_null,myfontsize,fontsize_subletters,modrange_col,modrange_edcol,meas_BAZ_floor4plot,meas_phiSC4plot,meas_dtSC4plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

% use selfdefined step graycolormap depending on how many models should be
% plotted, -1 because best model is plotted in red
%mycolormap = mgrund_custommap(plot_mod_max-1);

% subplot 1 phi
s1=subplot(2,1,1);

% plot model range only if not full range is used
if modrange_low~=0 && modrange_upp~=360 || modrange_low == 0 && modrange_upp~=360 ;
    
    % frist make whole background gray
    xdir=[0 0 360 360];
    ydir=[-90 90 90 -90];
    
    p1=patch(xdir,ydir,'white');
    set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)
    hold on

    xdir=[modrange_low modrange_low modrange_upp modrange_upp];
    ydir=[-90 90 90 -90];
    p2=patch(xdir,ydir,'white');
    set(p2,'edgecolor','w')
    hold on
end

% plot best models #2 -#max
for ii=2:plot_mod_max % 570:600
     plot(BAZ,BEST_models_sort(ii).phi_eff,'linewidth',linewidth_SS,'color',color_SS_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,BEST_models_sort(1).phi_eff,'linewidth',linewidth_SS,'color',color_SS_bf_1)
hold on

% %=========================
% %workaround to plot "axis" on top after real axis was set to bottom layer
% %to allow null circles plotted on top ;)
% hold on
% a = axis; %// gives xmin xmax ymin ymax
% cx = get(gca,'Xcolor'); %// color of x axis
% % plot([a(1) a(2)], [a(3)*0.9 a(3)*0.9], 'color', cx)
% plot([a(1) a(2)], [a(4)*0.9 a(4)*0.9], 'color', cx)
% 
% cy = get(gca,'Ycolor'); %// color of y axis
% plot([a(1) a(1)], [a(3) a(4)], 'color', cy)
% plot([a(2)*0.9 a(2)*0.9], [a(3) a(4)], 'color', cy)
% %=========================

  
% plot measured values
for FF=1:length(meas_phiSC4plot)
h1(FF)=errorbar(meas_BAZ_floor4plot(FF),meas_phiSC4plot(FF,2),abs(meas_phiSC4plot(FF,2)-meas_phiSC4plot(FF,1)),...
    abs(meas_phiSC4plot(FF,2)-meas_phiSC4plot(FF,3)),'ok','markerfacecolor',colorsfill(meas_phiSC4plot(FF,4),:),'markersize',...
    mymarkersize-1,'LineWidth',linewidth_symbols,'markeredgecolor',colorsedge(meas_phiSC4plot(FF,4),:),'color',colorsedge(meas_phiSC4plot(FF,4),:));
end

% set errorbar width
for FF=1:length(meas_phiSC4plot)
    errorbar_tick(h1(FF),6,'units') 
end


for ii=1:length(meas_phiSC_null(:,2))
% plot measured null values                               
plot(meas_BAZ_floor_null(ii),meas_phiSC_null(ii,2),'ok','markerfacecolor',color_face_null,...
    'markersize',mymarkersize_null-2,'LineWidth',linewidth_symbols_null,'markeredgecolor',color_edge_null);
end

ylabel('\phi in \circ','fontsize',myfontsize)
set(gca,'xticklabel',[])

text(0.012,0.97,'\bf(a)\rm' , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fontsize_subletters,'backgroundcolor','w','edgecolor','k');


set(gca,'fontsize',myfontsize)
set(gca,'xtick',-90:45:360)
set(gca,'ytick',-90:45:90)
set(gca, ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on')

set(gca,'layer','top')
box on

xlim([0 360]) 
ylim([-90 90])
 
set(gca,'TickDir','out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% plot RMS vlaue 
text(0.012,0.12,['\bfRMS_{tot} = ' num2str(BEST_models_sort(1).RMS,'%4.2f') ', RMS_{\phi} = ' num2str(BEST_models_sort(1).RMS_phi,'%4.2f') '\circ, RMS_{\deltat} = ' num2str(BEST_models_sort(1).RMS_dt,'%4.2f') ' s'], ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fontsize_subletters,'backgroundcolor','w','edgecolor','k','color',color_SS_bf_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%=============================================================================================
%=============================================================================================
%=============================================================================================

%  xx=meas_BAZ_floor';  
%  yy=[meas_phiSC(:,2)]';
% 
% 
%   MIN_XY=[0,-90];
%   MAX_XY=[360,90];
%   rangeX=0.2;
%   rangeY=0.2;
%   
%   ngrid=2^8
% 
%   %colormap(viridis(230)) 
% 
% %   %clear all  
% %   % generate synthetic data
% %   data=[xx', yy';];
% %   % call the routine, which has been saved in the current directory 
% %   [bandwidth,density,X,Y]=kde2d(data,ngrid,MIN_XY,MAX_XY,rangeX,rangeY);
% %   % plot the data and the density estimate
% %   pcolor(X,Y,density), hold on
% %   %plot(data(:,1),data(:,2),'r.','MarkerSize',5)
% %  shading interp;
% 
% 
% % if show_densplot==1
% % 
% % % include density plot
% % 
% % Border = 15;
% % Sigma = 10;
% % stepsize = 1;
% % 
% % xx=meas_BAZ_floor';  
% % yy=[meas_phiSC(:,2)]';
% % 
% % dd = [xx' yy'];
% % nn = length(xx);
% % 
% % xrange = [0 360];
% % yrange = [-90 90];
% % 
% % % setup coordinate grid
% % [XX YY] = meshgrid(xrange(1):stepsize:xrange(2), yrange(1):stepsize:yrange(2));
% % YY = flipud(YY);
% % 
% % %Parzen parameters and function handle
% % pf1 = @(C1,C2) (1/nn)*(1/((2*pi)*Sigma^2)).*...
% %          exp(-( (C1(1)-C2(1))^2+ (C1(2)-C2(2))^2)/(2*Sigma^2));
% % 
% % PPDF1 = zeros(size(XX));    
% % 
% % %Populate coordinate surface
% % [R C] = size(PPDF1);
% % NN = length(dd);
% % for c=1:C
% %    for r=1:R 
% %        for d=1:nn 
% %             PPDF1(r,c) = PPDF1(r,c) + ...
% %                 pf1([XX(1,c) YY(r,1)],[dd(d,1) dd(d,2)]); 
% %        end
% %    end
% % end
% % 
% % colormap(parula(230))
% % 
% % %Normalize data
% % m1 = max(PPDF1(:));
% % PPDF1 = PPDF1 / m1;
% % 
% % %Set up visualization
% % %set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',20)
% % %stem3(D(:,1),D(:,2),zeros(N,1),'b.','markersize',20,'markerfacecolor','red','markeredgecolor','red');
% % %hold on;
% % 
% % grid off
% % 
% % %Add PDF estimates to figure
% % sssss1 = pcolor(XX,YY,PPDF1);
% % 
% % shading interp;
% % %alpha(s1,0.1);
% % 
% % end


%=============================================================================================
%=============================================================================================
%=============================================================================================

%............................................................

% subplot 2 dt
s2=subplot(2,1,2);

% plot model range only if not full range is used
if modrange_low~=0 && modrange_upp~=360 || modrange_low == 0 && modrange_upp~=360 ;
    
    xdir=[0 0 360 360];
    ydir=[0 4 4 0];
    p1=patch(xdir,ydir,'white');
    set(p1,'facecolor',modrange_col,'edgecolor',modrange_edcol)
    hold on
    
    
    xdir=[modrange_low modrange_low modrange_upp modrange_upp];
    ydir=[0 4 4 0];
    p2=patch(xdir,ydir,'white');
    set(p2,'edgecolor','white')
    hold on
end

% plot best models #2 -#max
for ii=2:plot_mod_max 
     plot(BAZ,BEST_models_sort(ii).dt_eff,'linewidth',linewidth_SS,'color',color_SS_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,BEST_models_sort(1).dt_eff,'linewidth',linewidth_SS,'color',color_SS_bf_1)
hold on

% plot measured values
for FF=1:length(meas_phiSC4plot)
h2(FF)=errorbar(meas_BAZ_floor4plot(FF),meas_dtSC4plot(FF,2),abs(meas_dtSC4plot(FF,2)-meas_dtSC4plot(FF,1)),...
    abs(meas_dtSC4plot(FF,2)-meas_dtSC4plot(FF,3)),'ok','markerfacecolor',colorsfill(meas_phiSC4plot(FF,4),:),'markersize',...
    mymarkersize-1,'LineWidth',linewidth_symbols,'markeredgecolor',colorsedge(meas_phiSC4plot(FF,4),:),'color',colorsedge(meas_phiSC4plot(FF,4),:));
end

% set errorbar width
for FF=1:length(meas_phiSC4plot)
    errorbar_tick(h2(FF),6,'units') 
end



%=========================
% workaround to plot "axis" on top after real axis was set to bottom layer
% to allow null circles plotted on top ;)
hold on
a = axis; %// gives xmin xmax ymin ymax
cx = get(gca,'Xcolor'); %// color of x axis
plot([a(1) a(2)], [a(3) a(3)], 'color', cx)
plot([a(1) a(2)], [a(4) a(4)], 'color', cx)

cy = get(gca,'Ycolor'); %// color of y axis
plot([a(1) a(1)], [a(3) a(4)], 'color', cy)
plot([a(2)*0.9 a(2)*0.9], [a(3) a(4)], 'color', cy)
%=========================

% plot measured null values, all delay times of nulls are set to zero (0)!!!   
for ii=1:length(meas_dtSC_null(:,2))
hnulls=plot(meas_BAZ_floor_null(ii),[0],'ok','markerfacecolor',...
    color_face_null,'markersize',mymarkersize_null-2,'LineWidth',linewidth_symbols_null,...
    'markeredgecolor',color_edge_null);
end
   
xlabel('Backazimuth in \circ','fontsize',myfontsize)
ylabel('\deltat in s','fontsize',myfontsize)


text(0.012,0.97,'\bf(b)\rm' , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fontsize_subletters,'backgroundcolor','w','edgecolor','k');


set(gca,'fontsize',myfontsize)
set(gca,'xtick',0:45:360)
set(gca,'ytick',0:1:4)
set(gca, ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on')

set(gca,'layer','bottom')
box on

xlim([0 360])
ylim([0 4])

set(gca,'TickDir','out');
%............................................................

pos=get(s1,'Position');
set(s1,'Position',[pos(1) pos(2)+0 pos(3)*0.9 pos(4)]);    


pos=get(s2,'Position');
set(s2,'Position',[pos(1) pos(2)+0.1 pos(3)*0.9 pos(4)]); 

%............................................................



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE colored density plot

% to revise!!!




% if show_densplot==1
% %%
% % include density plot
% 
% Border = 0.5;
% Sigma = 10;
% stepsize = 1;
% 
% xx=meas_BAZ_floor';  
% yy=[meas_dtSC(:,2)]'
% 
% length(xx)
% 
% %   rangeX=0.1;
% %   rangeY=0.005;
% 
%     rangeX=0.2;
%   rangeY=0.2;
%   
% 
% %%
% 
% 
%   MIN_XY=[0,0];
%   MAX_XY=[360,4];
% 
%   %colormap(parula(230))
% 
% %   %clear all  
% %   % generate synthetic data
% %   data=[xx', yy';];
% %   % call the routine, which has been saved in the current directory 
% %   [bandwidth,density,X,Y]=kde2d(data,ngrid,MIN_XY,MAX_XY,rangeX,rangeY);
% %   % plot the data and the density estimate
% %   pcolor(X,Y,density), hold on
% %   %plot(data(:,1),data(:,2),'r.','MarkerSize',5)
% %  shading interp;
% 
% 
%  
% 
% %[bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% % dd = [xx' yy'];
% % nn = length(xx);
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % xrange = [min(xx)-Border max(xx)+Border];
% % % yrange = [min(yy)-Border max(yy)+Border];
% % 
% % xrange = [0 360];
% % yrange = [0 4];
% % 
% % % setup coordinate grid
% % [XX YY] = meshgrid(xrange(1):stepsize:xrange(2), yrange(1):stepsize:yrange(2));
% % YY = flipud(YY);
% % 
% % %Parzen parameters and function handle
% % pf1 = @(C1,C2) (1/nn)*(1/((2*pi)*Sigma^2)).*...
% %          exp(-( (C1(1)-C2(1))^2+ (C1(2)-C2(2))^17)/(2*Sigma^2));
% % 
% % PPDF1 = zeros(size(XX));    
% % 
% % %Populate coordinate surface
% % [R C] = size(PPDF1);
% % NN = length(dd);
% % for c=1:C
% %    for r=1:R 
% %        for d=1:nn 
% %             PPDF1(r,c) = PPDF1(r,c) + ...
% %                 pf1([XX(1,c) YY(r,1)],[dd(d,1) dd(d,2)]); 
% %        end
% %    end
% % end
% % 
% % colormap(parula(230))
% % 
% % %Normalize data
% % m1 = max(PPDF1(:));
% % PPDF1 = PPDF1 / m1;
% % 
% % %Set up visualization
% % %set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',20)
% % %stem3(D(:,1),D(:,2),zeros(N,1),'b.','markersize',20,'markerfacecolor','red','markeredgecolor','red');
% % %hold on;
% % 
% % grid off
% % 
% % %Add PDF estimates to figure
% % sssss1 = pcolor(XX,YY,PPDF1);
% % axis([xrange(1) xrange(2) yrange(1) yrange(2)])
% % 
% % 
% % shading interp;
% % %alpha(s1,0.1);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
