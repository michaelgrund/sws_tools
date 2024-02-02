function SWS_modeling_plot_results(BAZ,models_sort,plot_mod_max,...
    meas_BAZ_floor_null,meas_phiSC_null,meas_dtSC_null,...
    modrange_low,modrange_upp,colmod_bf_1,colmod_bf_2max,lw_mod,colfill,...
    coledge,ms,ms_null,lw_symbols,lw_symbols_null,color_face_null,...
    col_edge_null,fs,fs_RMSE,modrange_col,modrange_edcol,...
    meas_BAZ_floor,meas_phiSC,meas_dtSC,staname_split,nameend)

% plot shear-wave splitting modeling results
%
% LICENSE
%
% Copyright (C) 2020  Michael Grund, Karlsruhe Institute of Technology (KIT).
% ORCID: https://orcid.org/0000-0001-8759-2018
% GitHub: https://github.com/michaelgrund/sws_tools
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

%================================================================
% check your matlab version
vers_out=SWS_modeling_check_matlab_version();

%================================================================
% plot models together with measured data
figure()

%================================================================
% subplot 1 phi
s1=subplot(2,1,1);

% plot model range only if not full range is used
if (modrange_low~=0 && modrange_upp~=360) || ...
   (modrange_low==0 && modrange_upp~=360) || ...
   (modrange_low~=0 && modrange_upp==360)
    
    % first make whole background gray
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

% plot best models #2 - #max
for ii=2:plot_mod_max % 570:600
     plot(BAZ,models_sort(ii).phi_eff,'linewidth',lw_mod,...
         'color',colmod_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,models_sort(1).phi_eff,'linewidth',lw_mod,'color',colmod_bf_1)
hold on

% plot measured values

sizeSC = size(meas_phiSC);

for FF=1:sizeSC
    h1(FF)=errorbar(meas_BAZ_floor(FF),meas_phiSC(FF,2),...
        abs(meas_phiSC(FF,2)-meas_phiSC(FF,1)),...
        abs(meas_phiSC(FF,2)-meas_phiSC(FF,3)),...
        'ok','markerfacecolor',colfill(meas_phiSC(FF,4),:),...
        'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
        coledge(meas_phiSC(FF,4),:),'color',...
        coledge(meas_phiSC(FF,4),:));
end


if ~isempty(meas_phiSC_null)
    for ii=1:length(meas_phiSC_null(:,2))
        % plot measured null values
        plot(meas_BAZ_floor_null(ii),meas_phiSC_null(ii,2),...
            'ok','markerfacecolor',color_face_null,'markersize',ms_null-2,...
            'LineWidth',lw_symbols_null,'markeredgecolor',col_edge_null);
    end
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
% plot RMSE value
text(0.018,0.165,['\bfRMSE_{tot} = ' num2str(models_sort(1).RMSE,'%4.2f') ', RMSE_{\phi} = ' num2str(models_sort(1).RMSE_phi,'%4.2f') '\circ, RMSE_{\deltat} = ' num2str(models_sort(1).RMSE_dt,'%4.2f') ' s'], ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fs_RMSE,'backgroundcolor','w','edgecolor','k','color',colmod_bf_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================
% subplot 2 dt
s2=subplot(2,1,2);

% plot model range only if not full range is used
if (modrange_low~=0 && modrange_upp~=360) || ...
   (modrange_low==0 && modrange_upp~=360) || ...
   (modrange_low~=0 && modrange_upp==360)
    
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

% plot best models #2 - #max
for ii=2:plot_mod_max 
     plot(BAZ,models_sort(ii).dt_eff,'linewidth',lw_mod,'color',colmod_bf_2max) 
     hold on
end

% plot best model #1
plot(BAZ,models_sort(1).dt_eff,'linewidth',lw_mod,'color',colmod_bf_1)
hold on

% plot measured values
for FF=1:sizeSC
    h2(FF)=errorbar(meas_BAZ_floor(FF),meas_dtSC(FF,2),...
        abs(meas_dtSC(FF,2)-meas_dtSC(FF,1)),...
        abs(meas_dtSC(FF,2)-meas_dtSC(FF,3)),...
        'ok','markerfacecolor',colfill(meas_phiSC(FF,4),:),...
        'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
        coledge(meas_phiSC(FF,4),:),...
        'color',coledge(meas_phiSC(FF,4),:));
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
if ~isempty(meas_dtSC_null)
    for ii=1:length(meas_dtSC_null(:,2))
        plot(meas_BAZ_floor_null(ii),[0],'ok','markerfacecolor',...
            color_face_null,'markersize',ms_null-2,'LineWidth',lw_symbols_null,...
            'markeredgecolor',col_edge_null);
    end
end

%=========================
% create legend
% predefined colors for singles and multis
singlcol=[66 91 169]./256;

% splits and stacks in struct?
wcols=unique(colfill(meas_phiSC(:,4),:),'rows');
swcols=size(wcols);


% plot measured null values, all delay times of nulls are set to zero (0) 
if ~isempty(meas_dtSC_null)
    for ii=1:length(meas_dtSC_null(:,2))
        plot(meas_BAZ_floor_null(ii),[0],'ok','markerfacecolor',...
            color_face_null,'markersize',ms_null-2,'LineWidth',lw_symbols_null,...
            'markeredgecolor',col_edge_null);
    end
end



if swcols(1)==2 % both, singles and multis in dataset
    if isequal(wcols(1,:),singlcol)
        l1 = plot(-10,-10,'ok','markerfacecolor',wcols(1,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]);
        l2 = plot(-10,-10,'ok','markerfacecolor',wcols(2,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]); 
        
        handvec = [l1(1),l2(1)];
        strvec = {'split (s)','split (m)'};

    else
        l1 = plot(-10,-10,'ok','markerfacecolor',wcols(2,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]);
        l2 = plot(-10,-10,'ok','markerfacecolor',wcols(1,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]); 
        
        handvec = [l2(1),l1(1)];
        strvec = {'split (s)','split (m)'};
 
    end

else % only singles or multis in dataset
    if isequal(wcols(1,:),singlcol)
    
        l1 = plot(-10,-10,'ok','markerfacecolor',wcols(1,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]);
        
        handvec = [l1(1)];
        strvec = {'split (s)'};
 
    else
        l1 = plot(-10,-10,'ok','markerfacecolor',wcols(2,:),...
            'markersize',ms-1,'LineWidth',lw_symbols,'markeredgecolor',...
            [0 0 0]);
        
        handvec = [l1(1)];
        strvec = {'split (m)'};
    end

end

% add nulls to legend if in dataset
if ~isempty(meas_dtSC_null)
    n1 = plot(-10,-10,'ok','markerfacecolor',...
        color_face_null,'markersize',ms_null-2,'LineWidth',lw_symbols_null,...
        'markeredgecolor',col_edge_null);
    handvec(end+1)=n1;
    strvec{end+1}='null';

end

legend(handvec,strvec);

%=========================

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

% save file
filename=['PLOT_modeling_' staname_split '_' num2str(plot_mod_max) '_best_models_'...
    nameend '_range_' num2str(modrange_low) '_' num2str(modrange_upp)];

if vers_out == 1 %%% requires R2020a or later
    exportgraphics(gcf,[filename '.pdf'],'ContentType','vector')
else %%% if your version is below 2020a use
    print ('-dpdf', '-painters','-r600', [filename '.pdf'])
end

%%% if you work on a Linux machine, you can also try:
%print ('-depsc', '-painters','-r600', [filename '.eps'])
%dir_eps_file=dir([filename '.eps']);
%[status,cmdout]=system(['epstopdf ' dir_eps_file.name]);
%%
%================================================================
%================================================================
% plot overview and statistics of modeling
figure()

s11=subplot(2,1,1);
modtype1 = {models_sort.mod_type};
idx1 = find(strcmp(modtype1, 'dipping'));
idx2 = find(strcmp(modtype1, 'two_layers'));
idx3 = find(strcmp(modtype1, 'single_layer'));

models_sort(1).color = 0;

for kk = 1:length(idx1)
   models_sort(idx1(kk)).color = [0 0.4470 0.7410];
end

for kk = 1:length(idx2)
   models_sort(idx2(kk)).color = [0.6350 0.0780 0.1840];
end

for kk = 1:length(idx3)
   models_sort(idx3(kk)).color = [0.9290 0.6940 0.1250];
end

for ii = 1:length(models_sort)
    plot(ii, models_sort(ii).RMSE, 'color', models_sort(ii).color, ...
        'marker', 'o', ...
        'markerfacecolor',models_sort(ii).color,...
        'markersize', 10)
    hold on
end

% legend
ll1=plot(-10, -10, 'o', 'markeredgecolor', [0 0.4470 0.7410], ...
        'markerfacecolor',[0 0.4470 0.7410],...
        'markersize', 10);

hold on

ll2=plot(-10, -10, 'o', 'markeredgecolor', [0.6350 0.0780 0.1840], ...
        'markerfacecolor',[0.6350 0.0780 0.1840],...
        'markersize', 10);
    
hold on    
    
ll3=plot(-10, -10, 'o', 'markeredgecolor', [0.9290 0.6940 0.1250],...
        'markerfacecolor',[0.9290 0.6940 0.1250],...
        'markersize', 10);

legend([ll1,ll2,ll3],{'dipping layer','two-layers','single-layer'})           

xlim([0,length(models_sort)])
ylim([0,max([models_sort.RMSE])])
xlabel('worst \leftarrow sorted models \rightarrow best','fontsize',fs)
ylabel('RMSE_{tot}','fontsize',fs)


if vers_out == 1
    sgtitle([num2str(length(models_sort)) ' best models (based on RMSE)'],'fontsize',fs)
else
    title([num2str(length(models_sort)) ' best models (based on RMSE)'],'fontsize',fs)
end

set(gca,'XDir','reverse','fontsize',fs); 
set(gca, ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on')

grid on
set(gca,'TickDir','out');
set(gca,'fontsize',fs)

% add bar plot to show distribution
s22=subplot(2,1,2);

b1=bar(1, length(idx1));
hold on
b2=bar(2, length(idx2));
hold on
b3=bar (3, length(idx3));

set(b1,'FaceColor',[0 0.4470 0.7410]);
set(b2,'FaceColor',[0.6350 0.0780 0.1840]);
set(b3,'FaceColor',[0.9290 0.6940 0.1250]);

ylabel('model count','fontsize',fs)
set(gca,'xticklabel',[])
set(gca,'xtick',[])
set(gca,'TickDir','out');
set(gca,'fontsize',fs)

%================================================================
% save file
filename=['PLOT_modeling_' staname_split '_' num2str(plot_mod_max) '_modeldstr_'...
    nameend '_range_' num2str(modrange_low) '_' num2str(modrange_upp)];

if vers_out == 1 %%% requires R2020a or later
    exportgraphics(gcf,[filename '.pdf'],'ContentType','vector')
else %%% if your version is below 2020a use
    print ('-dpdf', '-painters','-r600', [filename '.pdf'])
end

%================================================================




