function PLOT_SCHEMA_2_layers()
% plot 2 layer models in 3D with synthetic raypaths
%
% 2018-03-14 -MG-
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% from:
% https://de.mathworks.com/matlabcentral/answers/342134-how-can-i-find-strike-and-dip-of-a-plane-from-equation-of-a-plane

finddipdir=1; % ~ 0°, however only for 1° a stable solution is found
finddip=30;


close all

%[controls; controls theta; controls dip]

%n = [1 2 3]';
n2 = [0 0 1]';

n_e = [0 0 1]'; % normal vector of earths surface pointing upwards in z direction

increment=1;

ZZ=1;
for ii=1:increment:100
    
    for jj=1:increment:100
        
        for zz=1:increment:100

        n = [ii jj zz]';  
            
        normal = n/norm(n); 
        
        dip = acosd(dot(normal,n_e));
        theta = atan2d(normal(1), normal(2)) + 90;
        downdipdir=theta-90;
        
        %round(dip)==finddip && round(downdipdir)==finddipdir
        if round(dip)==finddip && round(downdipdir)==finddipdir
        
            finalVALS(ZZ,:)=[ii jj zz];
            
            ZZ=ZZ+1;


        end
        
        end
        
    end
    
end
%%

% normal = n/norm(n); 
% n_e = [0 0 1]'; % normal vector of earths surface pointing upwards in z direction
% dip = acosd(dot(normal,n_e));
% theta = atan2d(normal(1), normal(2)) + 90;
% downdipdir=theta-90;



close all

n=finalVALS(1,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if selected values are working!!!
normal = n/norm(n); 
dip = acosd(dot(normal,n_e)) 
theta = atan2d(normal(1), normal(2)) + 90;
downdipdir=theta-90
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f1=figure(1); clf; hold on;
hV=drawVector(n, {'hallo'});% the normal

% lower layer
DIPdir1=0;
zval=-42;
laylowcol=[238,44,44]./256;

[hs2,v2]=drawPlane(n2, zval,laylowcol,DIPdir1);   % plot first the lower plane to get correct overstacks, the plane shifted by "5"
drawPlane(n2, zval-5,laylowcol,DIPdir1);
% upper layer thru origin
DIPdir2=30;
layuppcol=[16,78,139]./256;

[hs,v]=drawPlane(n,0,layuppcol,DIPdir2);               % unshifted plane, comes through the origin
drawPlane(n, -450,layuppcol,DIPdir2); 

%remove normal vector from plot!
set(hV.p,'visible','off')
set(hV.l,'visible','off')
set(hV.t,'visible','off')

XXX=v(:,:,1);
YYY=v(:,:,2);
ZZZ=v(:,:,3);

XXX2=v2(:,:,1);
YYY2=v2(:,:,2);
ZZZ2=v2(:,:,3);

%=============================================%=============================================
%=============================================%=============================================
%%% plot fast axis on layers

%%%%%% UPPER LAYER

% generate grid to get minimum and maximum values for later plotting the
% fast axis manually, as to be done again depending on the fast axis on the
% plane

%plot3(xq,yq,vq,'ow','markerfacecolor','r','markersize',4);
%
% m1=mesh(xq,yq,vq)
% 
%         direction=[0 0 1];
%         rotate(m1,direction,0)
       
%%%%%%%%%%%
% cases
% 1) in direction of 30° dip, so dipdir==fast dir==30°
[xq,yq] = meshgrid(-5:10:5, -40:80:40);
vq = griddata(XXX,YYY,ZZZ+0.5,xq,yq);

barcol=[255,127,0]./256;
l1=surf(xq,yq,vq,'FaceColor', barcol, 'EdgeColor', 'k','LineWidth', 1);

% 2) in direction of XX° dip, so dipdir==fast dir==XX°

%%%%%%%%%%%

% now both , plane and fast axis have a direction dipping to east
direction=[0 0 1];

%rotate to dip dir
rotate(l1,direction,360-DIPdir2)

%=============================================

%%%%%% LOWER LAYER

%%%%%%%%%%%
% cases
% 1) in direction of 0° (==N)
[xq,yq] = meshgrid(-5:10:5, -40:80:40);
vq = griddata(XXX2,YYY2,ZZZ2+0.5,xq,yq);

l2=surf(xq,yq,vq,'FaceColor', 'w', 'EdgeColor', 'k','LineWidth', 1);
% 2) in direction of XX° dip, so dipdir==fast dir==XX°

%%%%%%%%%%%

% now both , plane and fast axis have a direction dipping to east
direction=[0 0 1];

%rotate to dip dir
rotate(l2,direction,360)


%=============================================%=============================================
%=============================================%=============================================

% add dip direction 
ll1=line([0 0],[0 75],[0 0],'linestyle','--','color','k');
rotate(ll1,direction,360-DIPdir2)

% plot dip direction angle arc
ang([0 0],37.5,[0 degtorad(DIPdir2)],'k-')

% plot dip direction angle arc
ang([0 0],37.4,[0 degtorad(DIPdir2)],'k-')
ang([0 0],37.3,[0 degtorad(DIPdir2)],'k-')
ang([0 0],37.2,[0 degtorad(DIPdir2)],'k-')
ang([0 0],37.1,[0 degtorad(DIPdir2)],'k-')

% plot phi symbol
text(3,30,'\phi','fontsize',26,'color','w')


% plot dip angle arc (is manually so placed that it is hided beneath the
% dashed line when looking from top ;))
h=ang3D([0 0],31.6,[0 degtorad(finddip)],'k-');
rotate(h,direction,-DIPdir2)
h=ang3D([0 0],31.5,[0 degtorad(finddip)],'k-');
rotate(h,direction,-DIPdir2)
h=ang3D([0 0],31.4,[0 degtorad(finddip)],'k-');
rotate(h,direction,-DIPdir2)
h=ang3D([0 0],31.3,[0 degtorad(finddip)],'k-');
rotate(h,direction,-DIPdir2)
h=ang3D([0 0],31.2,[0 degtorad(finddip)],'k-');
rotate(h,direction,-DIPdir2)


%=============================================%=============================================
%=============================================%=============================================
%
%%% plot raypaths in circular order

[h,xunit,yunit,zunit] = circle3(0,0,-65, 30,'k');

staN=0;
staE=0;
staZ=50;

for ii=1:length(xunit)
    hold on
    plot3([xunit(ii) staN],[yunit(ii) staE], [zunit(ii) staZ],'k')
end


% PLOT station triangle on top of raypaths
stacol=[0,150,130]./256; % KIT grün ;)
plot3(0,0,staZ+2,'marker','^','markersize',20,'markerfacecolor',stacol,...
    'markeredgecolor','k','LineWidth',5)
plot3(0,0,staZ+2,'marker','^','markersize',19.5,'markerfacecolor',stacol,...
    'markeredgecolor','k','LineWidth',5)
plot3(0,0,staZ+2,'marker','^','markersize',19,'markerfacecolor',stacol,...
    'markeredgecolor','k','LineWidth',5)
plot3(0,0,staZ+2,'marker','^','markersize',18.5,'markerfacecolor',stacol,...
    'markeredgecolor','k','LineWidth',5)

% lims=100;
%zlim([-80 40])
%xlim([-lims lims])
%ylim([-lims lims])

%=============================================%=============================================
%=============================================%=============================================
% 3D cones around the x-y-z-axis

maxnXY=77.3;

r = 1.7015;
h = 6;

m = h/r;

hold on 

[R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));
X = -m*R;
Y = R .* cos(A);
Z = R .* sin(A);

XP=maxnXY;
YP=0;

hold on
surf(X+min(min(XP))+0.25,Y+min(min(YP))+0.02,Z)
colormap([0  0 0])


%%%%%%%%%%%%%%%%%%%
% y axis 

[R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));

X = R .* cos(A);
Y= -m*R;
Z = R .* sin(A);

XP=0;
YP=maxnXY;

hold on
hold on
surf(X+min(min(XP))+0.02,Y+min(min(YP))+0.25,Z)
colormap([0  0 0])

%%%%%%%%%%%%%%%%%%%%%%%%
% Z

[R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));

X = R .* cos(A);
Y = R .* sin(A);
Z = -m*R;

maxn=maxnXY;

XP=0;
YP=0;

hold on
surf(X+min(min(XP))+0.02,Y+min(min(YP))+0.02,Z+0.25+maxn)
colormap([0 0 0])
%box on

set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])

set(gca,'ticklength',[0 0])



%light('Position',[55 55 55],'Style','local')
%light('Position',[-55 -55 55],'Style','local')
%light('Position',[-55 -55 0],'Style','local')

set(gca,'fontsize',14)

%=============================================%=============================================
%=============================================%=============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% necessary for box line plotting
valxlim=xlim();
valylim=ylim();
valxlim=zlim();

%=========================
% plot in different view directions

clear az ele

az(1)=42;
ele(1)=19;

az(end+1)=106;
ele(end+1)=13;

az(end+1)=0;
ele(end+1)=90;

az(end+1)=62;
ele(end+1)=-2;

for ii=3%1:length(az)
    
    view(az(ii),ele(ii))
    
    % manually add lines for "open box look", depending on the view direction
    if az(ii)==42 && ele(ii) == 19

    % this is only valid for view(42,19)
        plot3([valxlim(2) valxlim(2)],[valxlim(2) valxlim(2)], [-valxlim(2) valxlim(2)],'k')
        plot3([valxlim(2) -valxlim(2)],[valxlim(2) valxlim(2)], [valxlim(2) valxlim(2)],'k')
        plot3([-valxlim(2) -valxlim(2)],[valxlim(2) -valxlim(2)], [valxlim(2) valxlim(2)],'k')
    
    elseif az(ii)==106 && ele(ii) == 13
   
    elseif az(ii)==0 && ele(ii) == 90
        
        box on
    
    else
    
    end

    filename=['PLOT_3D_2layers_view_az' num2str(az(ii)) '_ele' num2str(ele(ii))];
    print ('-dpng','-r800', [filename '.png']) 
    [status,cmdout]=system(['convert ' filename '.png -trim ' filename '_trim.png']);

end

%=========================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%print(gcf, '-djpg', 'exFIG.jpg'); 
% 
% export_fig exFIG -png -transparent
% 
% print ('-dpng','-r450', [filename '.png']) 
 
% set(f1, 'Renderer','opengl') 
 

%  print ('-dpdf','-r600', [filename '.pdf']) 
 
%savepdf(1,'TESTsavpdf')


%plot2svg(filename,f1,'png')








