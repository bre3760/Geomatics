clc
clear all
pos_sat = load('pos_sat.txt');
% latitude and longitude of receiver station
lat = 45 + 3/60 + 48/3600;
lon = 7 + 39/60 + 41/3600;  
h = 0;
wgs84 = wgs84Ellipsoid('meter');
[x,y,z] = geodetic2ecef(wgs84,lat,lon,h);
D = zeros(size(pos_sat,1),4);
for r =1:size(pos_sat,1)
    roh = sqrt( (pos_sat(r,2) - x)^2 + (pos_sat(r,3) - y)^2 + (pos_sat(r,4) - z)^2 );
    D1= (pos_sat(r,2)-x)/roh;
    D2= (pos_sat(r,3)-y)/roh;
    D3= (pos_sat(r,4)-z)/roh;
    D4 = -1;
    D(r,1) = D1;
    D(r,2) = D2;
    D(r,3) = D3;
    D(r,4) = D4;
end
Qxx = inv(transpose(D)*D)
R = [-sind(lon) cosd(lon) 0;
    -sind(lat)*cosd(lon) -sind(lat)*sind(lon) cosd(lat);
    cosd(lat)*cosd(lon) cosd(lat)*sind(lon) sind(lat) ];
Qxx_star = Qxx(1:3,1:3);
Quu = R*Qxx_star*R'
xdop = (Quu(1,1))
ydop = (Quu(2,2))
vdop = (Quu(3,3))
TDOP = (Qxx(4,4))
HDOP = sqrt(xdop+ydop)
PDOP = sqrt(xdop + ydop + vdop)
GDOP = sqrt(PDOP^2 + TDOP^2)