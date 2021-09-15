
function T = transformation_func(theta , d , a , alpha) 
   T = [cosd(theta)  -sind(theta)*cosd(alpha) sind(theta)*sind(alpha) a*cosd(theta) ;
       sind(theta) cosd(theta)*cosd(alpha) -cosd(theta)*sind(alpha) a*sind(theta);
       0 sind(alpha) cosd(alpha) d;
       0 0 0 1]
  % X= T01 .* T12 .* T23 .* T34 .* T45 .* T56 .* T67


end
