function SWS_modeling_plot_results(BAZ,models_sort,plot_mod_max,...
    meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,colmod_bf_1,colmod_bf_2max,lw_mod,colfill,...
    coledge,ms,ms_null,lw_symbols,lw_symbols_null,color_face_null,...
    col_edge_null,fs,fs_RMSE,modrange_col,modrange_edcol,...
    meas_BAZ_floor4plot,meas_phiSC4plot,meas_dtSC4plot)

% plot modeling results

%================================================================
figure(1)

%================================================================
% subplot 1 phi
s1=subplot(2,1,1);

% plot model range only if not full range is used
if modrange_low~=0 && modrange_upp~=360 || modrange_low == 0 && modrange_upp~=360
    
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
     plot(BAZ,models_sort(ii).phi_eff,'linewidth',lw_mod,...
         'color',colmod_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,models_sort(1).phi_eff,'linewidth',lw_mod,'color',colmod_bf_1)
hold on

% plot measured values
for FF=1:length(meas_phiSC4plot)
    h1(FF)=errorbar(meas_BAZ_floor4plot(FF),meas_phiSC4plot(FF,2),...
        abs(meas_phiSC4plot(FF,2)-meas_phiSC4plot(FF,1)),...
        abs(meas_phiSC4plot(FF,2)-meas_phiSC4plot(FF,3)),...
        'ok','markerfacecolor',colfill(meas_phiSC4plot(FF,4),:),...
        'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
        coledge(meas_phiSC4plot(FF,4),:),'color',...
        coledge(meas_phiSC4plot(FF,4),:));
end

for ii=1:length(meas_phiSC_null(:,2))
    % plot measured null values                               
    plot(meas_BAZ_floor_null(ii),meas_phiSC_null(ii,2),...
        'ok','markerfacecolor',color_face_null,'markersize',ms_null-2,...
        'LineWidth',lw_symbols_null,'markeredgecolor',col_edge_null);
end

ylabel('\phi in \circ','fontsize',fs)
set(gca,'xticklabel',[])
set(gca,'fontsize',fs)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot RMS vlaue 
text(0.018,0.18,['\bfRMSE_{tot} = ' num2str(models_sort(1).RMS,'%4.2f') ', RMSE_{\phi} = ' num2str(models_sort(1).RMS_phi,'%4.2f') '\circ, RMSE_{\deltat} = ' num2str(models_sort(1).RMS_dt,'%4.2f') ' s'], ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fs_RMSE,'backgroundcolor','w','edgecolor','k','color',colmod_bf_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================
% subplot 2 dt
s2=subplot(2,1,2);

% plot model range only if not full range is used
if modrange_low~=0 && modrange_upp~=360 || modrange_low == 0 && modrange_upp~=360
    
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
     plot(BAZ,models_sort(ii).dt_eff,'linewidth',lw_mod,'color',colmod_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,models_sort(1).dt_eff,'linewidth',lw_mod,'color',colmod_bf_1)
hold on

% plot measured values
for FF=1:length(meas_phiSC4plot)
    h2(FF)=errorbar(meas_BAZ_floor4plot(FF),meas_dtSC4plot(FF,2),...
        abs(meas_dtSC4plot(FF,2)-meas_dtSC4plot(FF,1)),...
        abs(meas_dtSC4plot(FF,2)-meas_dtSC4plot(FF,3)),...
        'ok','markerfacecolor',colfill(meas_phiSC4plot(FF,4),:),...
        'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
        coledge(meas_phiSC4plot(FF,4),:),...
        'color',coledge(meas_phiSC4plot(FF,4),:));
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

% plot measured null values, all delay times of nulls are set to zero (0)   
for ii=1:length(meas_dtSC_null(:,2))
    hnulls=plot(meas_BAZ_floor_null(ii),[0],'ok','markerfacecolor',...
    color_face_null,'markersize',ms_null-2,'LineWidth',lw_symbols_null,...
    'markeredgecolor',col_edge_null);
end
   
xlabel('Backazimuth in \circ','fontsize',fs)
ylabel('\deltat in s','fontsize',fs)

set(gca,'fontsize',fs)
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

%================================================================
% adjust figure size and position

pos=get(s1,'Position');
set(s1,'Position',[pos(1) pos(2)+0 pos(3)*0.9 pos(4)]);    

pos=get(s2,'Position');
set(s2,'Position',[pos(1) pos(2)+0.1 pos(3)*0.9 pos(4)]); 

%================================================================