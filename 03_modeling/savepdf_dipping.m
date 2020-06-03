function savepdf_dipping(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% savepdf(fignumber,filename,image)
% - filename as string without '.pdf'
% - fignumnber: number of figure you want to plot
% -in case your plot has a lot of images use the optional argument image
%  and set it to 1!! -> otherwise the resolution is really bad 
%   
%##########################################################################################################
%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
%#%													 %#
%#%   Copyright (c) 2009 by the KaSP-Team.								 %#
%#%   This file is part of the Karlsruhe Seismology Processing (KaSP) Toolbox for MATLAB!		 %#
%#%													 %#
%#%   The KaSP toolbox is free software under the terms of the GNU General Public License!		 %#
%#%													 %#
%#%   Please see the copyright/licensce notice distributed together with the KaSP Toolbox!		 %#
%#%													 %#
%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
%##########################################################################################################
%
% 2010-JAN
%
% Tobias Baumann, Karlsruhe Institute of Technology, Geophysical Institute
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 if nargin < 2
     error('there are too less arguments')
 end
 
 if nargin == 2
     fignumber = varargin{1};
     filename  = varargin{2};
     image     = 0;
 end

 if nargin == 3
     fignumber = varargin{1};
     filename  = varargin{2};
     image     = varargin{3};
 end

 if nargin > 3
     error('there are too many arguments')
 end
 

 
set(0,'CurrentFigure',fignumber);
set(gcf, 'PaperOrientation','landscape')
set(gcf, 'PaperType','A4');
papersize = get(gcf, 'PaperSize');

width=29;
height=20;
% width=30;
% height=22;

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height-2)/2;

myfiguresize = [left, bottom+2, width, height];
set(gcf, 'PaperPosition', myfiguresize);
if image == 1
    print ('-dpdf', '-painters','-r600', [filename '.pdf'])
else
    print('-dpdf',[filename '.pdf']);
end
return