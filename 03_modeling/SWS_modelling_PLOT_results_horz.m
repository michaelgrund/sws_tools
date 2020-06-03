function SWS_modelling_PLOT_results_horz()


figure(1);

% subplot 1 phi
s1=subplot(1,2,1);


% plot model range

xdir=[modrange_low modrange_low modrange_upp modrange_upp];
ydir=[-90 90 90 -90];
patch(xdir,ydir,'red')
hold on



% plot best models #2 -#max
for ii=2:plot_mod_max % 570:600
     plot(BAZ,BEST_models_sort(ii).phi_eff,'linewidth',linewidth_SS,'color',color_SS_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,BEST_models_sort(1).phi_eff,'linewidth',linewidth_SS,'color',color_SS_bf_1)
hold on
  
% plot measured values

h1=errorbar(meas_BAZ_floor,[meas_phiSC(:,2)],[abs(meas_phiSC(:,2)-meas_phiSC(:,1))],...
    [abs(meas_phiSC(:,2)-meas_phiSC(:,3))],'dk','markerfacecolor',color_face,'markersize',...
    mymarkersize-1,'LineWidth',linewidth_symbols,'markeredgecolor',colorsedge(meas_phiSC(:,4),:),'color',color_error);

errorbar_tick(h1,2,'units') 
hold on

for ii=1:length(meas_phiSC_null(:,2))
% plot measured null values                               
plot(meas_BAZ_floor_null(ii),meas_phiSC_null(ii,2),'ok','markerfacecolor',color_face_null,...
    'markersize',mymarkersize_null-2,'LineWidth',linewidth_symbols_null,'markeredgecolor',color_edge_null);
end

xlabel('Backazimuth in \circ','fontsize',myfontsize)
ylabel('\phi in \circ','fontsize',myfontsize)

text(0.012,0.97,'\bf(a)\rm' , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fonsize_subletters,'backgroundcolor','w','edgecolor','k');


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
 
%............................................................

% subplot 2 dt
s2=subplot(1,2,2);
 
% plot best models #2 -#max
for ii=2:plot_mod_max 
     plot(BAZ,BEST_models_sort(ii).dt_eff,'linewidth',linewidth_SS,'color',color_SS_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,BEST_models_sort(1).dt_eff,'linewidth',linewidth_SS,'color',color_SS_bf_1)
hold on

% plot measured values
h2=errorbar(meas_BAZ_floor,[meas_dtSC(:,2)],[abs(meas_dtSC(:,2)-meas_dtSC(:,1))],...
    [abs(meas_dtSC(:,2)-meas_dtSC(:,3))],'dk','markerfacecolor',color_face,'markersize',...
    mymarkersize-1,'LineWidth',linewidth_symbols,'markeredgecolor',color_edge,'color',color_error);

errorbar_tick(h2,2,'units') 
      
% plot measured null values   
for ii=1:length(meas_dtSC_null(:,2))
plot(meas_BAZ_floor_null(ii),[meas_dtSC_null(ii,2)],'ok','markerfacecolor',...
    color_face_null,'markersize',mymarkersize_null-2,'LineWidth',linewidth_symbols_null,...
    'markeredgecolor',color_edge_null);
end
   

xlabel('Backazimuth in \circ','fontsize',myfontsize)
ylabel('\deltat in s','fontsize',myfontsize)


text(0.012,0.97,'\bf(b)\rm' , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fonsize_subletters,'backgroundcolor','w','edgecolor','k');


set(gca,'fontsize',myfontsize)
set(gca,'xtick',0:45:360)
set(gca,'ytick',0:1:4)
set(gca, ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on')


xlim([0 360])
ylim([0 4])
 
%............................................................

pos=get(s1,'Position');
set(s1,'Position',[pos(1) pos(2) pos(3)*1.2 pos(4)*0.4]);    

pos=get(s2,'Position');
set(s2,'Position',[pos(1)+0.02 pos(2) pos(3)*1.2 pos(4)*0.4]); 

%............................................................    



