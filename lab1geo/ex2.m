Xsat = 15487292.829 ;
Ysat = 6543538.932 ;
Zsat = 20727274.429 ;
lat = 45 + 3/60 + 48.114/3600;
lon = 7 + 39/60 + 40.605/3600;

[x,y,z] = geodetic2ecef(wgs84Ellipsoid, lat,lon,0);

dx = Xsat - x;
dy = Ysat - y;
dz = Zsat - z;

R = [-sind(lon) cosd(lon) 0;
    -sind(lat)*cosd(lon) -sind(lat)*sind(lon) cosd(lat);
    cosd(lat)*cosd(lon) cosd(lat)*sind(lon) sind(lat) ];

%local coordinates

c = R*[dx;dy;dz];
e = c(1);
n = c(2);
u = c(3);
elev = atand(e/n);
azi = atand(u/sqrt(n^2 + e^2));





