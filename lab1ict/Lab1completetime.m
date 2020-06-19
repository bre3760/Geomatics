% Version including all time instants
clc 
clear all 
close all
load('/Users/brendanpolidori/Desktop/poli/geomatics/labs/lab1/data/DataSet/NominalUERE/dataset_5_20180328T121529.mat')

% %GPSNASS
% [satellite,time] = size(RHO.GPS);
% 
% %---PART 1---
% %Plot the number of satellite for each time istant
% number_sat = zeros(1,3600);
% for i=1:time
%     number_sat(1,i) = length(find(not(isnan(RHO.GPS(:,i)))));
% end 
% x = 1:length(number_sat);
% figure(1);
% scatter(x,number_sat, 'LineWidth', 1);
% xlim( [1 3600]);
% ylim([2 12]);
% title('Number of visible satellites in each time instant');
% xlabel('time instants');
% ylabel('number of satellites');
% grid;
% 
% % mean(number_sat)
% % Plot the pseudoranges measured by each satellites during the hour for the
% % first database using constellation GPS
% figure(2);
% hold all
% for i=1:satellite
% %     y = (not(isnan(RHO.GPS(i,:))));
%     
%     y = RHO.GPS(i,:);
%     plot(x,y, 'LineWidth',2);
%     xlim( [1 3600]);
%     grid;
% end
% title('Pseudoranges for each satellite in each time instant');
% xlabel('time instants');
% ylabel('pseudorange');
% hold off
% 
% %---PART 2---
% 
% %%% matrix of the estimated position at each time instant
% x_vector = zeros(time,4);
% 
% for t=1:time
%     
%     %initial estimated position for each time instant
%     x_estimated = [0 0 0 1];
%     
%     %visible Satellites for each time instant
%     visibleSatellites = find(not(isnan(RHO.GPS(:,t))));
%     
%     %vector that will containt the calculated pseudoranges
%     rho_hat = zeros(1, length(visibleSatellites));
%     %mesaured pseudoranges (extracted from the RHO.GPS visible satellites)
%     rho = zeros(1, length(visibleSatellites));
%     for s=1:length(visibleSatellites)
%       
%         rho(s) = RHO.GPS(visibleSatellites(s),t);
%     end
%     
%     %Initialization H matrix
%     H = zeros(length(visibleSatellites),4);
%     H(:,4) = 1;
%     delta_rho = zeros(length(visibleSatellites),1);
%     
%     for k=1:9
%     
%     %for each satellite we find its coordinates 
%     %calculate the pseudorange given the user estimated position
%     %update the row of the H matrix corresponding to that satellite
%     
%         for s=1:length(visibleSatellites)
%             sat_coord = SAT_POS_ECEF.GPS(visibleSatellites(s)).pos(t,:);
%             rho_hat(s) = sqrt(sum(((sat_coord) - x_estimated(1:3)).^2)) + x_estimated(4);
%             r(s) = sqrt(sum(((sat_coord) - x_estimated(1:3)).^2));
%             H(s,1:3) = ((sat_coord) - x_estimated(1:3))/r(s);
%         end
% 
%         delta_rho = (rho_hat - rho)';
%         pseudoinv = (transpose(H)*H)\transpose(H);
%         delta_pos = pseudoinv * delta_rho; 
%         x_estimated = x_estimated + delta_pos'; %update of user location
%     
%     end
%     
%     x_vector(t,:) = x_estimated;
%    
% end
% 
% % x_vector has some NaN values, i will remove them
% % indexes_to_remove = find(isnan(x_vector(:,1)));
% % for i = 1:size(x_vector,1)
% %     j=1;
% %     if i == indexes_to_remove(j)
% %         x_vector(i)=eps;
% %         j=j+1;
% %     end
% % end
% 
% 
% 
% x_standard_deviation = std(x_vector(:,1))
% y_standard_deviation = std(x_vector(:,2))
% z_standard_deviation = std(x_vector(:,3))
% b_standard_deviation = std(x_vector(:,4))
% 
% avg_user_pos = mean(x_vector,1);
% lla = ecef2lla(avg_user_pos(1:3));
% 
% writeKML_GoogleEarth('Final', lla(1,1), lla(1,2), lla(1,3));
% 
% % % %SKYPLOT OF SATELLITES
% % % %testing for one satellite 
% % % sat_number = 22;
% % % wgs84 = wgs84Ellipsoid;
% % % % Specify the geodetic coordinates of the local origin. 
% % % % In this example, the local origin is the satellite dish. 
% % % % Specify h0 as ellipsoidal height in kilometers.
% % % lla_vect = ecef2lla(x_vector(:,1:3)); 
% % % 
% % % 
% % % lat0 = lla_vect(:,1); %degrees
% % % lon0 = lla_vect(:,2); %degrees
% % % h0 = lla_vect(:,3); %meters
% % % % Specify the ECEF coordinates of the point of interest. In this example,
% % % % the point of interest is the satellite.
% % % x = SAT_POS_ECEF.GPS(sat_number).pos(:,1);
% % % y = SAT_POS_ECEF.GPS(sat_number).pos(:,2);
% % % z = SAT_POS_ECEF.GPS(sat_number).pos(:,3);
% % % % Then, calculate the AER coordinates of the satellite 
% % % % with respect to the satellite dish. In this example,
% % % % slantRange displays in scientific notation.
% % % [az,elev,slantRange] = ecef2aer(x,y,z,lat0,lon0,h0,wgs84);
% % % 
% % % 
% % % 
% % % figure()
% % % polarplot(az,slantRange)
% % 
% % 
% % 
% % 
% %SKYPLOT OF SATELLITES
% 
% wgs84 = wgs84Ellipsoid;
% lla_vect = ecef2lla(x_vector(:,1:3)); 
% 
% lat0 = lla_vect(:,1); %degrees
% lon0 = lla_vect(:,2); %degrees
% h0 = lla_vect(:,3); %meters
% 
% x = zeros(3600,size(visibleSatellites,1));
% y = zeros(3600,size(visibleSatellites,1));
% z = zeros(3600,size(visibleSatellites,1));
% for ind = 1:size(visibleSatellites,1)
%     sat_j = visibleSatellites(ind);
%     x(:,ind) = SAT_POS_ECEF.GPS(sat_j).pos(:,1);
%     y(:,ind) = SAT_POS_ECEF.GPS(sat_j).pos(:,2);
%     z(:,ind) = SAT_POS_ECEF.GPS(sat_j).pos(:,3);
% end
% 
% for ind = 1:size(visibleSatellites,1)
% [az(:,ind),elev(:,ind),slantRange(:,ind)] = ecef2aer(x(:,ind),y(:,ind),z(:,ind),lat0,lon0,h0,wgs84);
%     
%     for i = 1:3600
%         if x(i,ind) == 0
%                 az(i,ind) = 0;
%         end
%     end
% end
% start_time = 2800;
% end_time = 3600;
% 
% figure()
% 
% % polarplot(az(start_time:end_time,:),slantRange(start_time:end_time,:),'LineWidth',2)
% polarplot(az(start_time:end_time,:),elev(start_time:end_time,:),'LineWidth',2)
% lgd = legend(string(visibleSatellites));
% pax = gca;
% angles = 0:30:360;
% pax.ThetaTick = angles;
% labels = {'E','60','30','N','330','300','W','240','210','S','150','120'};
% pax.ThetaTickLabel = labels;
% 
% title(lgd,'visible satellites')












































% with realistic data
clear all
close all
clc
load('/Users/brendanpolidori/Desktop/poli/geomatics/labs/lab1/data/DataSet/RealisticUERE/dataset_5_20180329T161418.mat')
% %GPSNASS
[satellite,time] = size(RHO.GPS);

%---PART 1---
%Plot the number of satellite for each time istant
number_sat = zeros(1,3600);
for i=1:time
    number_sat(1,i) = length(find(not(isnan(RHO.GPS(:,i)))));
end 
x = 1:length(number_sat);
figure(1);
scatter(x,number_sat, 'LineWidth', 1);
xlim( [1 3600]);
ylim([2 12]);
grid;


%Plot the pseudoranges measured by each satellites during the hour for the
%first database using constellation GPS
figure(2);
hold all
for i=1:satellite
    y = RHO.GPS(i,:);
    plot(x,y, 'LineWidth',2);
    xlim( [1 3600]);
    grid;
end
hold off


errorsintime = zeros(size(RHO.GPS,1),3598);
for sat_j = 1:size(RHO.GPS,1)
    errorsintime(sat_j,:) = diff(RHO.GPS(sat_j,:),2);
end
% figure(3);
% plot(linspace(1,3598,3598),errorsintime(12,:),'-o',linspace(1,3598,3598),errorsintime(17,:),'-x')
% lgd = legend('12','17');
% title(lgd,'Satellites')
% title('Errors in time')
% xlabel('time')
% ylabel('m')
varianceoferror = zeros(32,1);
for i = 1: size(errorsintime,1)
    varianceoferror(i) = std(errorsintime(i,:)); 
end




%matrix of the estimated position at each time instant
x_vector = zeros(time,4);

for t=1:time
    
    %initial estimated position for each time instant
    x_estimated = [0 0 0 1];
    
    %visible Satellites for each time instant
    visibleSatellites = find(not(isnan(RHO.GPS(:,t))));
    
    %vector that will containt the calculated pseudoranges
    rho_hat = zeros(1, length(visibleSatellites));
    %mesaured pseudoranges (extracted from the RHO.GPS visible satellites)
    rho = zeros(1, length(visibleSatellites));
    for s=1:length(visibleSatellites)
        rho(s) = RHO.GPS(visibleSatellites(s),t);
    end
    
    %Initialization H matrix
    H = zeros(length(visibleSatellites),4);
    H(:,4) = 1;
    delta_rho = zeros(length(visibleSatellites),1);
    
    R = eye(size(visibleSatellites,1));
    varsatvis = zeros(size(visibleSatellites,1),1);
    
    for i = 1:size(visibleSatellites,1)
        index = visibleSatellites(i);
        varsatvis(i) = varianceoferror(index);
    end
    
    R = R.*varsatvis;
    
    W = inv(R);
    
    
    for k=1:9
    
    %for each satellite we find its coordinates 
    %calculate the pseudorange given the user estimated position
    %update the row of the H matrix corresponding to that satellite
    
        for s=1:length(visibleSatellites)
            sat_coord = SAT_POS_ECEF.GPS(visibleSatellites(s)).pos(t,:);
            rho_hat(s) = sqrt(sum(((sat_coord) - x_estimated(1:3)).^2)) + x_estimated(4);
            r(s) = sqrt(sum(((sat_coord) - x_estimated(1:3)).^2));
            H(s,1:3) = ((sat_coord) - x_estimated(1:3))/r(s);
        end

        delta_rho = (rho_hat - rho)';
        pseudoinv_w = (transpose(H)*W*H)\transpose(H)*W;
        delta_pos = pseudoinv_w*delta_rho;
        x_estimated = x_estimated + delta_pos'; %update of user location
    
    end
    
    x_vector(t,:) = x_estimated;
   
end
x_standard_deviation = std(x_vector(:,1))
y_standard_deviation = std(x_vector(:,2))
z_standard_deviation = std(x_vector(:,3))
b_standard_deviation = std(x_vector(:,4))
avg_user_pos = mean(x_vector,1);
lla = ecef2lla(x_estimated(1:3));

writeKML_GoogleEarth('Final', lla(1,1), lla(1,2), lla(1,3));

