function P = plot_arc(a,b,h,k,r,colfill)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians, 
% b is end of arc in radians, 
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig

t = linspace(a,b);
x = r*cos(t) + h;
y = r*sin(t) + k;
x = [x h x(1)];
y = [y k y(1)];
P = fill(x,y,'r');
set(P,'facecolor',colfill,'edgecolor','white','facealpha',1)
axis([h-r-1 h+r+1 k-r-1 k+r+1]) 
axis tight;


if ~nargout
    clear P
end
