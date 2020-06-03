function [dist]=distance_in_dipping_layer(dip,aoi,thick,azi)

%  ** calculate apparent dip
      alp = atand(tand(dip) .* sind(azi-90)) ;
      
%  ** calculate distances
      gam = 90 - alp - aoi ;
      bet = 90 + alp ;
      dist = thick .* sind(bet) ./ sind(gam) ;

return