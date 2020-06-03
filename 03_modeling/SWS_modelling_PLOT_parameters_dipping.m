function SWS_modelling_PLOT_parameters_dipping(dipdirs,dips,thicks,staname_split)
% make roseplots for parameters from modelling
% and count the determined dip angles and down-dip direction






%clear all
clc
close all




%fillcol='r';
%edgecol=[178 34 34]./256;

coldef=viridis(10);
fillcol=coldef(6,:);
edgecol=coldef(1,:);

fontsizeall=16;

figure();

subplot(1,2,1)

clear theta

%theta = deg2rad([0 78  80 78  80  78 80 78  80 78  80  78 80 78 40 40 40 45 45 45 45 45 45  23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 ])

%recalculate dip dirs

dipdirs=dipdirs-90;

theta = deg2rad(dipdirs);


% all roseplots have a fixed radial axis up to 20, since this corresponds
% to the maximum number of the best plotted models, therefore all plots are
% dircetly comparable

theta_fake = deg2rad(ones(1,20)*10);

% whole circle 0:360 is dividided into 36 bins (each 10°)
n = 36;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hax = axes();
hax = subplot(2,2,1)

rose(hax,theta_fake,n); % generate axis length in radial dir with fake_theta
hold on
rose(hax,theta,n); % Use this to initialise polar axes

cla
[theta,rho] = rose(hax,theta,n); % Get the histogram edges
theta(end+1) = theta(1); % Wrap around for easy interation
rho(end+1) = rho(1); 
hold on;
for j = 1:floor(length(theta)/4)
    k = @(j) 4*(j-1)+1; % Change of iterator
    h = polar(theta(k(j):k(j)+3),rho(k(j):k(j)+3));
    %set(h,'color',colours(j,:)); % Set the color
    set(h,'color',edgecol); % Set the colo
    [x,y] = pol2cart(theta(k(j):k(j)+3),rho(k(j):k(j)+3));
    hh = patch(x,y,'');
    %set(h,'FaceColor',colours(j,:));
    set(hh,'FaceColor',fillcol);
    uistack(hh,'down');
end

grid on 
axis equal;

% val1=get(hax,'XLim');
% val2=get(hax,'YLim');
% set(hax,'XLim',[0-0.0125 val1(2)-0.0125]);
% set(hax,'YLim',[0-0.0125 val2(2)-0.0125]);

set(hax,'ydir','reverse');

th = findall(hax,'Type','text');
lh = findall(hax,'Type','line');

% Define what to keep
legit = {'0','30','60','90','120','150','180','210','240','270','300','330','360'};

ZZ=1;
for ii = 1:length(th),
    set(th(ii),'FontSize',fontsizeall)
    if ii > 16 % place radial scale 
        testget(ZZ,:)=get(th(ii),'position');
        maxval=max(theta);
        set(th(ii),'position',[-5-5/maxval testget(ZZ,2)-0.4 testget(ZZ,3)])
        ZZ=ZZ+1;
        
    else  

      % Take the others and set them to empty string
      idx = ismember(get(th(ii),'string'),legit);
      if idx ~=0
         
        %===================  
        % reorder labels to get 0° to north  
        currstr=get(th(ii),'string');
        currstrnum=str2num(currstr);
        
        if currstrnum-270 <= -270 || currstrnum-270 <= -180
            currstr=num2str(currstrnum-270+360);
        elseif currstrnum-270 <= -30 || (currstrnum-270 > -180 && currstrnum-270 < 0)
            currstr=num2str(currstrnum-270+360);
            
        else
            currstr=num2str(currstrnum-270);
        end
         %===================  
        
        set(th(ii),'string',[currstr '\circ' ])  
        testget(ZZ,:)=get(th(ii),'position');
        set(th(ii),'position',[testget(ZZ,1)*1.07 testget(ZZ,2)*1.05 testget(ZZ,3)])
      else
        %set(th(ii),'string',[legit2{ii} '\circ' ]) 
        testget(ZZ,:)=get(th(ii),'position');
      end
        
    end
    
end

% for ii = 5:16
%     set(th(ii),'string',[legit2{ii}])   
% end





for ii = 1:length(lh),
     set(lh(ii),'linewidth',1.5)
     set(lh(ii),'linestyle','-')
     set(lh(ii),'color','k')
end

hline = findobj(gca,'Type','line');
set(hline,'LineWidth',2)
set(hline,'color',edgecol)

pos1=get(hax,'position');
set(hax,'position',[pos1(1)-0.03 pos1(2)-0.105 pos1(3)*1.3 pos1(4)*1.3])

text(-25,-26,'\bf(a)','fontsize',23)
text(25,-26,'\bf(b)','fontsize',23)
text(-25,28,'\bf(c)','fontsize',23)

text(24.3, 30, staname_split,'HorizontalAlignment','Center','VerticalAlignment','Base','fontsize',20,'backgroundcolor','w','edgecolor','k');


text(-25,-27,'                  dip direction ','fontsize',18)
text(25,-27,'                            dip angle','fontsize',18)
  
%view([-90 90])


hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hax = axes();
hax2 = subplot(2,2,2);

clear theta

%theta = deg2rad([80 78 80 78 80 78  80 78  80  7 8 78 80 78  80 78  80   80 78  80 78 80 78  80 78  80  78 80 78  80 78  80  78 80 78 40 40 40 45 45 45 45 45 45  23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 ])

theta = deg2rad(dips);

% all roseplots have a fixed radial axis up to 20, since this corresponds
% to the maximum number of the best plotted models, therefore all plots are
% dircetly comparable
theta_fake2 = deg2rad(ones(1,20)*10);

% whole circle 0:360 is dividided into 36 bins (each 10°)
n = 36;

rose(hax2,theta_fake2,n); % generate axis length in radial dir with fake_theta
hold on
rose(hax2,theta,n); % Use this to initialise polar axes

cla
[theta,rho] = rose(hax2,theta,n); % Get the histogram edges
theta(end+1) = theta(1); % Wrap around for easy interation
rho(end+1) = rho(1); 
hold on;
for j = 1:floor(length(theta)/4)
    k = @(j) 4*(j-1)+1; % Change of iterator
    h = polar(theta(k(j):k(j)+3),rho(k(j):k(j)+3));
    %set(h,'color',colours(j,:)); % Set the color
    set(h,'color',edgecol); % Set the colo
    [x,y] = pol2cart(theta(k(j):k(j)+3),rho(k(j):k(j)+3));
    hh = patch(x,y,'');
    %set(h,'FaceColor',colours(j,:));
    set(hh,'FaceColor',fillcol);
    uistack(hh,'down');
end

grid on 
axis equal;

val1=get(hax2,'XLim');
val2=get(hax2,'YLim');
set(hax2,'XLim',[0-0.0505 val1(2)+.1]);
set(hax2,'YLim',[0-0.0505 val2(2)-0.0125]);

set(hax2,'ydir','reverse');

th = findall(hax2,'Type','text');
lh = findall(hax2,'Type','line');

% Define what to keep
legit = {'0','30','60','90'};

ZZ=1;
for ii = 1:length(th),
    set(th(ii),'FontSize',fontsizeall)
    if ii > 16 % place radial scale to y axis
        testget(ZZ,:)=get(th(ii),'position');
        maxval=max(theta);
        set(th(ii),'position',[-2.2-5/maxval testget(ZZ,2)+0.7 testget(ZZ,3)])
        ZZ=ZZ+1;
        
    else  

      % Take the others and set them to empty string
      idx = ismember(get(th(ii),'string'),legit);
      if idx ~=0
        set(th(ii),'string',[get(th(ii),'string') '\circ' ])
        testget(ZZ,:)=get(th(ii),'position');
        set(th(ii),'position',[testget(ZZ,1)*0.98 testget(ZZ,2)*0.98 testget(ZZ,3)])
      else
        set(th(ii),'string',' ')
      end
        
    end
    
end

for ii = 1:length(lh),
     set(lh(ii),'linewidth',1.5)
     set(lh(ii),'linestyle','-')
     set(lh(ii),'color','k')
end

hline = findobj(gca,'Type','line');
set(hline,'LineWidth',2)
set(hline,'color',edgecol)

hold on


pos2=get(hax2,'position');
%set(hax2,'position',[pos2(1)-0.022 pos2(2)-0.2914 pos2(3)*1.15 pos2(4)*1.15])
set(hax2,'position',[pos2(1)-0.022 pos2(2)-0.132 pos2(3)*1.3 pos2(4)*1.3])


%hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hax3=subplot(2,2,[3,4]);

thicks2comp=0:25:500;

thick_count(1:21)=0;

for ii=1:length(thicks)
    for jj=1:length(thicks2comp)
        if thicks(ii)==thicks2comp(jj)
            thick_count(jj)=thick_count(jj)+1;
        end
    end
end
%



colormap(viridis)
imagesc(thick_count)
set(gca,'xtick',1:2:21)
set(gca,'TickLength',[0 0])
set(gca,'xticklabel',{0:50:500})
set(gca,'fontsize',fontsizeall)
xlabel('thickness in km','fontsize',fontsizeall)
set(gca,'yticklabel',[])

pos3=get(hax3,'position');
set(hax3,'position',[pos3(1) pos3(2)+0.2 pos3(3)*0.95 pos3(4)*0.2])

for ii=1:20
    
   hold on
   plot([ii+0.5 ii+0.5],[-3 3],'-w') 
    
end

cb=colorbar('northoutside');
set(cb,'fontsize',fontsizeall)
poscb=get(cb, 'Position');
set(cb,'Position',[poscb(1)+0.5 poscb(2)-0.03 poscb(3)*0.3 poscb(4)*1.1]);
xlabel(cb','count');

savename=['PLOT_modelling_DIPPING_' staname_split '_dip_rose'];
savepdf_dipping(1,savename)
system(['pdfcrop ' savename '.pdf ' savename '.pdf'])












