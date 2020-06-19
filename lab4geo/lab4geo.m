% % table has (a , alpha)
clear all 
clc
% refsys = ["delambre";"everest";"bessel";"fisher";"clarke";"helmert";"hayford";"krassovsky";"wgs84" ];
% table = [ 6376985 ,1/308.6;
%           63777276 ,1/300.8;
%           6377397 ,1/299.1528128;
%          6378338 ,1/288.5;
%           6378249 ,1/293.5;
%           6378140 ,1/298.3;
%            6378388, 1/297;
%           6378245, 1/298.3;
%           6378137 ,1/298.257223563];
% 
% 
% for i = 1:size(table,1)
%     a = table(i,1)
%     alpha = table(i,2)
%     c(i) = a*(1-alpha)
%     e_prime_2(i) = (a^2- c(i)^2)/c(i)^2;
%     e2(i) = (a^2- c(i)^2)/a^2;
% end
% 
% 
% % ex 2
% 
% p1_lat = 44 + 45/60 + 01.03930/3600 %phi
% p1_lon = 7 +24/60 +  29.20335/3600 %lambda
% p1_h = 322.4909
% 
% p2_lat = 44 + 47/60 + 10.90505/3600
% p2_lon = 7 +30/60 +  26.53939/3600
% p2_h = 305.7367
% 
% e2_wgs = e2(9)
% a_wgs = table(9,1)
% 
% e2_hay = e2(7)
% a_hay = table(7,1)
% 
% w_1_wgs=sqrt(1-e2_wgs * (sin(p1_lat))^2 )
% w_1_hay=sqrt(1-e2_hay * (sin(p1_lat))^2 )
% 
% 
% x_1_wgs = (a_wgs/w_1_wgs +p1_h)*cosd(p1_lat)*cosd(p1_lon)
% y_1_wgs = (a_wgs/w_1_wgs +p1_h)*cosd(p1_lat)*sind(p1_lon)
% z_1_wgs = (a_wgs/w_1_wgs*(1-e2_wgs) +p1_h)*sin(p1_lat)
% 
% x_1_hay = (a_hay/w_1_hay +p1_h)*cosd(p1_lat)*cosd(p1_lon)
% y_1_hay = (a_hay/w_1_hay +p1_h)*cosd(p1_lat)*sind(p1_lon)
% z_1_hay = (a_hay/w_1_hay*(1-e2_hay) +p1_h)*sin(p1_lat)
% 
% 
% w_2_wgs=sqrt(1-e2_wgs * (sin(p2_lat))^2 )
% w_2_hay=sqrt(1-e2_hay * (sin(p2_lat))^2 )
% 
% x_2_wgs = (a_wgs/w_2_wgs +p2_h)*cosd(p2_lat)*cosd(p2_lon)
% y_2_wgs = (a_wgs/w_2_wgs +p2_h)*cosd(p2_lat)*sind(p2_lon)
% z_2_wgs = (a_wgs/w_2_wgs*(1-e2_wgs) +p2_h)*sin(p2_lat)
% 
% x_2_hay = (a_hay/w_2_hay +p2_h)*cosd(p2_lat)*cosd(p2_lon)
% y_2_hay = (a_hay/w_2_hay +p2_h)*cosd(p2_lat)*sind(p2_lon)
% z_2_hay = (a_hay/w_2_hay*(1-e2_hay) +p2_h)*sin(p2_lat)
% 
% %changing the elevation of h1
% h1_2000 = 322.4909 + 2000
% 
% x_1_wgs_2000 = (a_wgs/w_1_wgs +p1_h)*cosd(p1_lat)*cosd(p1_lon)
% y_1_wgs_2000 = (a_wgs/w_1_wgs +p1_h)*cosd(p1_lat)*sind(p1_lon)
% z_1_wgs_2000 = (a_wgs/w_1_wgs*(1-e2_wgs) +p1_h)*sin(p1_lat)
% 
% % can't see a difference by using x,y,z
% 
% 
% % ex 3
% % from x,y,z to lat,lon, h
% 
% p1_x = 4499525.4271;
% p1_y = 585034.1293;
% p1_z = 4467910.3596;
% 
% p2_x = 4495694.2695;
% p2_y = 592457.8605;
% p2_z =4470744.7781;
% 
% p3_x = 4503484.7172;
% p3_y = 578160.7507;
% p3_z = 4465024.3002;
% 
% p4_x = 4498329.3715;
% p4_y = 562840.7651;
% p4_z = 4472537.6125;
% 
% % lambda_1 = atand(p1_y/p1_x); 
% lambda_4 = atand(p3_y/p3_x); 
% 
% converged = 0;
% iter = 0;
% while converged == 0
% %     r = sqrt(p1_x^2 + p1_y^2);
%     r = sqrt(p4_x^2 + p4_y^2);
%     if iter == 0
% %         phi_1 = atand(p1_z/r);
%         phi_1 = atand(p4_z/r);
% 
%     end
%     w = sqrt(1-e2_wgs*sind(phi_1)^2);
%     n = a_wgs/w;
% %     h = p1_x/(cosd(phi_1)*cosd(lambda_1)) - n;
%     h = p4_x/(cosd(phi_1)*cosd(lambda_4)) - n;
% 
%     phi_2 = atand(p2_z/(r*(1-(e2_wgs*n)/(n+h))));
%     if phi_2-phi_1<10^-8
%         converged = 1;
%         break;
%     end
%     phi_1 = phi_2;
%     iter = iter +1;
% end
% 
% phi_1
% h 
% lambda_4
%ex 4 
% use coordinates found in ex 2 
% input ETRF89
%  4499328.330154904 585008.5025092314  4408529.380210702
% output ITRF89
%  4499328.17400  585008.69870 4408529.51350      



% helmert transformation
helm_pts = load("points_helmert.txt");
ozz = [1,0,0];
zoz = [0,1,0];
zzo = [0,0,1];
helm_pts = helm_pts./1000000;
A = zeros(size(helm_pts,1),7);
I = zeros(size(helm_pts,1),1);
for r = 1:size(helm_pts,1)
    if mod(r,3) == 1
        A(r,1:3) = ozz;
        A(r,4:7) =  [helm_pts(r,1),0 ,-helm_pts(r,3),helm_pts(r,2)];
        I(r) = helm_pts(r,4)-helm_pts(r,1);
    elseif mod(r,3) == 2
        A(r,1:3) = zoz;
        A(r,4:7) =  [helm_pts(r,2),-helm_pts(r,3),0,-helm_pts(r,1)];
        I(r) = helm_pts(r,5)-helm_pts(r,2);
    elseif mod(r,3) == 0
        A(r,1:3) = zzo;
        A(r,4:7) =  [helm_pts(r,3),helm_pts(r,2),helm_pts(r,1),0];
        I(r) = helm_pts(r,6)-helm_pts(r,3);
    end
end


x = (transpose(A)*A)\(transpose(A)*I)
 
residuals = A*x - I


