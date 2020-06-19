clear all
close all
clc



% p1_lat = 45 + 3/60 + 45.717/3600;
% p1_long = 7 + 47/60 + 26.292/3600;
% p1_latr = deg2rad(p1_lat);
% p1_longr = deg2rad(p1_long);
% 
% p1_cm = 9;
% p1_utm_zone = 32;

p1_lat = 38 + 332/60 + 34.649/3600;
p1_long = 16 + 50/60 + 06.493/3600;
p1_latr = deg2rad(p1_lat);
p1_longr = deg2rad(p1_long);

p1_cm = 15;
p1_utm_zone = 33;


refsys = ["wgs84"];
table = [ 6378137 ,1/298.257223563];
a = 0;
e2 = 0;
e_prime_2 = 0;
for i = 1:size(table,1)
    a = table(i,1);
    alpha = table(i,2);
    c(i) = a*(1-alpha);
    e_prime_2(i) = (a^2- c(i)^2)/c(i)^2;
    e2(i) = (a^2- c(i)^2)/a^2;
end

%Radius of polar curvature
Rp = a^2/c;

% then we evaluate lambda prime which is the difference between
% your longitude and the value of the sample meridian 

lambda_prime = p1_long - (p1_cm); % this is in degrees
lambda_prime_r = deg2rad(lambda_prime);
% lambda_prime = p2_longr - p2_cm;

% longitude = lambda
% latitude = phi
% v1 and v

v1 = sqrt(1+ e_prime_2 * (cos(p1_latr))^2);
% v1 = sqrt(1+ e_prime_2 * (cos(p2_latr))^2);

% define chsi
squiggle = atan((tan(p1_latr))/(cos(v1*lambda_prime_r)));
% squiggle = atan((tan(p2_latr))/(cos(v1*lambda_prime)));

v = sqrt(1+ e_prime_2 * (cos(squiggle))^2);

A1 = 1 - e2/4 - (3*e2^2)/(64) - (5*e2^3)/(256);
A2 = (3*e2)/8 + (3*e2^ 2)/32 + (45*e2^3)/1024;
A4 = (15*e2^2)/256 + (45*e2^3)/1024;
A6 = (35*e2^3)/3072;

% to then estimate east and north we have to multiply 
% x and y by the contraction modulus and only for the east we have
% to add the false east due to the origin 

x = Rp * asinh((cos(squiggle)*tan(lambda_prime_r))/(v));
y = a*(A1*squiggle -A2*sin(2*squiggle) + A4*sin(4*squiggle) - A6*sin(6*squiggle) );
% what is the false east 500 km 
mc = 0.9996;


east = x*mc +500000
north = y*mc






