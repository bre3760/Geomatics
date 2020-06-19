clear all
clc
close all
latitudes = [30,37,45,60];
for index = 1:size(latitudes,2)
    lat = latitudes(index);

% long = zeros(6,1);
count = 1;
for i=0:0.5:3 
    long(count) = i;
    longr(count) = deg2rad(long(count));
    count = count+1;
    
end

latr = deg2rad(lat);
for i = 1:7

ml(i) = 0.9996*(1+((longr(i)^2)/(2) )*(cos(latr))^2) ;

end

figure(index)
plot(long,ml,'-o')
title('Modulus of Linear deformation')
grid on

end
