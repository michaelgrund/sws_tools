function SWS_histogram(phisall,dtsall)
%
% This function generates a histogram showing the distribution of the fast
% axes and delay times included in the full data set published by
% Grund & Ritter (2019). 
%
% 2019-05-02 -MG- (michael.grund@kit.edu)
%
% see also function: SWS_read_evstruct
%===============================================================================

col_fill=[0.2407    0.2978    0.5401];
fontsize_main=10;
lwsmooth=1.5;
lwcol='r';
lwstyle='-';
colglobmean=[0 205 102]./256;
fonsize_subletters=9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHI (fast axis)

binsphi=-90:5:90;

s1=subplot(1,2,1);
hist(phisall,binsphi);

%====================
% plot smoothed envelope
C=histc(phisall,binsphi);
windowSize = 4; 
b=(1/windowSize)*ones(1,windowSize);
a=1.1;

B=filter(b,a,C);
hold on
plot(binsphi,B,'color',lwcol,'linewidth',lwsmooth,'linestyle',lwstyle)
%====================

set(gca,'TickDir','out','linewidth',1.2,'fontsize',fontsize_main,'TickLength',[.01 .01])

h=findobj(gca,'Type','patch');
set(h,'FaceColor',col_fill,'EdgeColor','k');

xlim([-90 90])
set(gca,'xtick',-90:45:90)

text(0.018,0.979,['\bfall regions\rm'] , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fonsize_subletters,'backgroundcolor','w','edgecolor','k');

ylim([0 112])
xlabel('\phi in \circ')
ylabel('count')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DT (delay time)

bins2=0:0.1:2.5;

s2=subplot(1,2,2);
hist(dtsall,bins2);
a=hist(dtsall,bins2);

% global mean
hold on
plot([1 1],[0 max(a)+30],'linestyle','--','linewidth',2,'color',colglobmean);

xlim([0 3])
ylim([0 250])

%====================
% plot smoothed envelope
C=histc(dtsall,bins2);
windowSize = 4; 
b=(1/windowSize)*ones(1,windowSize);
a=1.1;

B=filter(b,a,C);
hold on
plot(bins2,B,'color',lwcol,'linewidth',lwsmooth,'linestyle',lwstyle)
%====================

set(gca,'TickDir','out','linewidth',1.2,'fontsize',fontsize_main,'TickLength',[.01 .01])

h=findobj(gca,'Type','patch');
set(h,'FaceColor',col_fill,'EdgeColor','k');

text(1.1,230,'ACR','color',colglobmean)

text(0.788,0.979,['N = ' num2str(length(dtsall))] , ...    
'Units', 'normalized', ...   
'HorizontalAlignment', 'left', ...
'VerticalAlignment', 'top','fontsize',fonsize_subletters,'backgroundcolor','w','edgecolor','k');

xlabel('\deltat in s')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust positions

xpos=0.05;
stretch=1;
stretchy=0.6;

pos=get(s1,'position');
set(s1,'position',[pos(1) pos(2) pos(3)*stretch pos(4)*stretchy])

pos=get(s2,'position');
set(s2,'position',[pos(1)-xpos pos(2) pos(3)*stretch pos(4)*stretchy])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save plots

filename='PLOT_RESULTS_histo'; 

print ('-dpdf', '-painters','-r600', [filename '.pdf']) 
 
print ('-depsc', '-painters','-r600', [filename '_cut.eps']) 
dir_eps_file=dir([filename '_cut.eps']);
[status,cmdout]=system(['epstopdf ' dir_eps_file.name]);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
