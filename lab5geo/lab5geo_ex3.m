
clear all
clc
east = 470139.66;
north = 5031468.37;

mc = 0.9996;

refsys = ["wgs84"];
table = [ 6378137 ,1/298.257223563];

for i = 1:size(table,1)
    a = table(i,1);
    alpha = table(i,2);
    c(i) = a*(1-alpha);
    e_prime_2(i) = (a^2- c(i)^2)/c(i)^2;
    e2(i) = (a^2- c(i)^2)/a^2;
end

E1 = (a-c)/(a+c);

B2 = (3*E1)/2 - (27*E1^3)/32;
B4 = (21*E1^2)/16 - (55*E1^4)/32;
B6 = (151*E1^3)/96;
B8 = (1097*E1)/512;

A1 = 1 - e2/4 - (3*e2^2)/(64) -(5*e2^3)/(256);
y = north/mc;
x =(east-500000)/mc;


teta = y/(a*A1);

squiggle = teta + B2*sin(2*teta) + B4 * sin(4*teta) + B6*sin(6*teta) + B8*sin(8*teta) ; 
v = sqrt(1+ e_prime_2 * (cos(squiggle))^2);
Rp = a^2/c;
% these two are in radians for now --> convert to degrees
lambda_prime = atan((v*sinh(x/Rp))/(cos(squiggle)));
lambda_prime_r = rad2deg(lambda_prime);

phi = atan(tan(squiggle)*cos(v*lambda_prime_r));

lambda_mc = 9;
% be careful that to have the correct longitude we have to add the central meridian value 
lambda_true = lambda_prime + lambda_mc;

phi_deg = rad2deg(phi)
lam_deg = (lambda_true)

phi_dms=degrees2dms(phi_deg)
lam_dms=degrees2dms(lam_deg)
