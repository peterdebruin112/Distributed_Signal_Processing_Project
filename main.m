clear all
close all

%% Consider the situation where we randomly place sensors in a square area of 100x100 m 2 to
% measure a certain quantity (for example temperature, density, ...).
rng(5) % Set random seed for repeatability
area = 100; % Meter (area x area)
n = 50; % Number of nodes
lower = 10;
upper = 30;

coordinates = area*rand([n 2]);
measurment = lower + (upper - lower).*rand(n,1);
figure
scatter(coordinates(:, 1), coordinates(:, 2))
%% • Design a proper connected sensor network that covers the area of the plant by using a
% reasonable number of sensors. How many do we need to guarantee a connected sensor
% network? Motivate your choice.
sensorRadius = 20; % Meter
nNeeded = ceil(find_n(sensorRadius)*area);
r = sqrt((2*log(n))/(n))*area;
prob = 1 - 1/n^2;
sensorRadius = r;
distanceSensors = zeros([n n]);
inRadius = zeros([n n]);
neighbors = [];
for i = 1:n
    for j = 1:n
        distanceSensors(i,j) = sqrt((coordinates(i,1) - coordinates(j,1))^2 + (coordinates(i,2) - coordinates(j,2))^2);
        if distanceSensors(i,j) <= sensorRadius && i<j
            inRadius(i,j) = 1;
            neighbors = [neighbors ; i , j];
        end
    end
end

figure
hold on
title(['# Sensors = ' num2str(n), ' Radius = ', num2str(sensorRadius), ' prob = ', num2str(prob)])
scatter(coordinates(:,1), coordinates(:,2))
for k = 1:length(neighbors)

    x_temp1 = coordinates(neighbors(k,1),1); % Get first x coordinate
    y_temp1 = coordinates(neighbors(k,1),2); % Get first y coordinate

    x_temp2 = coordinates(neighbors(k,2),1); % Get second x coordinate
    y_temp2 = coordinates(neighbors(k,2),2); % Get second y coordinate
    % Plot lines between neighboring nodes
    plot([x_temp1, x_temp2] ,[y_temp1, y_temp2], '--', Color='k')
end    
hold off


%% • Suppose the sensor network would like to compute the average value of the 
% measurement data. In addition to a randomised gossip implementation, which will
% serve as a baseline method, implement the average consensus problem using the
% PDMM algorithm. Report the performance in terms of convergence speed and number
% of transmissions and compare this to the results obtained by the randomised gossip
% algorithm.

% Randomized Gossip


%% • Suppose the sensor network would like to compute the median of the measurement
% data. Implement the median consensus problem using the PDMM algorithm.

