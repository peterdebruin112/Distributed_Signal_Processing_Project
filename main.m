clear all
close all

%% Create virtual area with randomly placed sensors
% Consider the situation where we randomly place sensors in a square area of 100x100 m 2 to
% measure a certain quantity (for example temperature, density, ...).

% Set random seed for repeatability
rng(5);

% Define the area, amount of sensors and upper and lower bound of the
% measured values from the sensors in the network.
area = 100; % Meter (area x area)
n = 101;     % Number of nodes
lower = -10;% Lower bound of sensor measurement
upper = 30; % Upper bound of sensor measurement

% Define the coordinates of the sensors and place them in the area.
coordinates = area*rand([n 2]);
measurment = lower + (upper - lower).*rand(n,1);
figure(1)
scatter(coordinates(:, 1), coordinates(:, 2))

%% Construct a proper connected sensor network
% • Design a proper connected sensor network that covers the area of the plant by using a
% reasonable number of sensors. How many do we need to guarantee a connected sensor
% network? Motivate your choice.

% Given the maximal coverage distance between consecutive sensors, compute
% the minimum amount of sensors needed to guarantee a proper connected sensor
% network.
sensorRadius = 10;                          % Define sensor coverage in m
nNeeded = ceil(find_n(sensorRadius)*area);  % 
r = sqrt((2*log(n))/(n))*area;
prob = 1 - 1/n^2;
sensorRadius = r;
distanceSensors = zeros(n);
inRadius = zeros(n);
neighbors = [];
neighbors_pdmm = [];
D = zeros(n);
A = zeros(n);

for i = 1:n
    for j = 1:n
        distanceSensors(i,j) = sqrt((coordinates(i,1) - coordinates(j,1))^2 + (coordinates(i,2) - coordinates(j,2))^2);
        if distanceSensors(i,j) <= sensorRadius && i<j
            inRadius(i,j) = 1;
            neighbors = [neighbors ; i , j]; 
        end
        if distanceSensors(i,j) <= sensorRadius && i~=j
            neighbors_pdmm = [neighbors_pdmm ; i , j]; 
            A(i,j) = 1;
            D(i,i) = D(i,i) + 1;
        end  
    end
end

if nnz(diag(D)) ~= length(diag(D))
    disp("ERROR")
    return
end

figure(2)
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
%selectedNode = ceil(1 + (n - 1)*rand(1));
W_i_j = zeros(n);
%W_Bar = zeros([n*n n]);
% Initialize 
P = generate_p_matrix(inRadius);

W_bar = zeros(n);
for i = 1:n
    for j = i+1:n 
        if P(i,j) > 0
            e_i = zeros([n 1]);
            e_i(i,:) = 1;
            e_j = zeros([n 1]);
            e_j(j,:) = 1;
            W_i_j = eye(n) - (1/2)*(e_i - e_j)*(e_i - e_j)';
            W_bar = W_bar + P(i,j)*W_i_j;
        end
    end
end    

numberEdges = size(neighbors,1);

W_cell = cell(numberEdges,1);

for k = 1:numberEdges
    i = neighbors(k,1);
    j = neighbors(k,2);
    e_i = zeros(n,1); e_i(i) = 1;
    e_j = zeros(n,1); e_j(j) = 1;
    W_cell{k} = eye(n) - 0.5 * (e_i - e_j) * (e_i - e_j)';
end

cvx_begin
    variable p(numberEdges) nonnegative

    W_bar = 0;
    for k = 1:numberEdges
        W_bar = W_bar + p(k) * W_cell{k};
    end    
    W_bar = W_bar / n;
    %eigval = eigs(W_bar);
    minimize (lambda_max(W_bar - (1/n)*ones(n)))
    %minimize (lambda_sum_largest(W_bar,2))
    %minimize (svd)

    for i = 1:n
        idx = (neighbors(:,1) == i) | (neighbors(:,2) == i);
        sum(p(idx)) == 1;
    end    
cvx_end
W_bar = full(W_bar);
P_opt = zeros(n);
for k = 1:numberEdges
    i = neighbors(k,1);
    j = neighbors(k,2);
    P_opt(i,j) = p(k);
    P_opt(j,i) = p(k); % symmetric
end
%%
x = measurment;
K = 1000000;
meanBase = mean(x);
error = zeros([K+1 1]);
transmissions_gossip = 0;
for k = 1:K
    error(k,1) = norm(x - meanBase,2)^2/n;
    index = randi([1 numberEdges], 1);
    % pickedNode = randi([1 n], 1);  %Pick uniformly a random node
    % edge_idx = randsample(n, 1, true, P_opt(pickedNode,:)); % pick edge node according to optimal P_opt
    % i = pickedNode;%neighbors(edge_idx, 1);
    % j = edge_idx;%neighbors(edge_idx, 2);
    % avg = (x(i) + x(j))/2;
    % x(i) = avg;
    % x(j) = avg;
    x = W_cell{index}*x;
    transmissions_gossip = transmissions_gossip + 2;
end    
error(k+1,1) = norm(x - meanBase,2)^2/n;

% figure(3)
% plot(error)
% xticks(0:1000:transmissions_gossip)
% set(gca, 'YScale', 'log')
fprintf('Final error: %f \n', error(end))

%% PDMM
A_pdmm = -triu(A) + tril(A);

x_pdmm = measurment;
a = measurment;
meanBase_pdmm = mean(x_pdmm);
K = 100; % Number of iterations
c = 0.2;
error_pdmm = zeros([K  1]);
Z = zeros(n);
%error_pdmm = [];
updated = zeros(n);
Y = zeros(n);
%d = diag(D); % Degree of every node
transmission_pdmm_mean = 0;
for k = 1:K
    error_pdmm(k,1) = norm(x_pdmm - meanBase_pdmm,2)^2/n;
    for i = 1:n  % Synchronous
        sumNeighbors = 0;
        pdmm_neighbors = find(neighbors_pdmm(:,1) == i);
        d = length(pdmm_neighbors);
        for index = pdmm_neighbors' % Transpose such that it iterates through it
            j = neighbors_pdmm(index, 2);
            if A_pdmm(i,j) == 0
                disp('Error: A_pdmm == 0')
                return
            end
            updated(i,j) = updated(i,j) + 1;
            sumNeighbors = sumNeighbors + A_pdmm(i,j)*Z(i,j);
        end
    % Primal Update equation 
        x_pdmm(i) = (a(i) - sumNeighbors)/ (1+ c*d);
        for index = pdmm_neighbors'
            j = neighbors_pdmm(index, 2);
            Y(i,j) = Z(i,j) + 2*c*A_pdmm(i,j)*x_pdmm(i); %Y update equation
            transmission_pdmm_mean = transmission_pdmm_mean + 1;
            %Z(j,i) = Y(i,j); %Z update equation
        end    
    end   
    % for i = 1:n
    %     pdmm_neighbors = find(neighbors_pdmm(:,1) == i);
    %     for index = pdmm_neighbors'
    %         j = neighbors_pdmm(index, 2);
    %         Y(i,j) = Z(i,j) + 2*c*A_pdmm(i,j)*x_pdmm(i); %Y update equation
    %         transmission_pdmm_mean = transmission_pdmm_mean + 1;
    % 
    %     end
    % %transmission_pdmm_mean = transmission_pdmm_mean + 1;
    % end
    for i = 1:n
        pdmm_neighbors = find(neighbors_pdmm(:,1) == i);
        for index = pdmm_neighbors'
            j = neighbors_pdmm(index, 2);
            Z(j,i) = Y(i,j);
        end
    end 
end

%error_pdmm(k + 1,1) = norm(x_pdmm - meanBase_pdmm,2)^2/n;

% figure(5)
% hold on
% grid on
% plot(error_pdmm)
% 
% legend('Randomised Gossip', 'PDMM')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'linear')
% ylim([10e-15 10e5])
% hold off
fprintf('Final error: %f \n', error_pdmm(end))
%% • Suppose the sensor network would like to compute the median of the measurement
% data. Implement the median consensus problem using the PDMM algorithm.

medianBase_pdmm = median(measurment);

A_pdmm_medi = -triu(A) + tril(A);

x_pdmm_medi = measurment;
a_medi = measurment;

K = 10000; % Number of iterations
c = 0.5;
error_pdmm_medi = zeros([K 1]);
Z_medi = zeros(n);
updated = zeros(n);
Y_medi = zeros(n);
%d = diag(D); % Degree of every node
transmissions_pdmm_medi = 0;
for k = 1:K
    error_pdmm_medi(k,1) = norm(x_pdmm_medi - medianBase_pdmm,2)^2/n;
    for i = 1:n 
        sumNeighbors_medi = 0;

        pdmm_neighbors = find(neighbors_pdmm(:,1) == i);
        d_medi = length(pdmm_neighbors);
        for index = pdmm_neighbors' % Transpose such that it iterates through it
            j = neighbors_pdmm(index, 2);
            if A_pdmm_medi(i,j) == 0
                disp('Error: A_pdmm == 0')
                return
            end
            updated(i,j) = updated(i,j) + 1;
            sumNeighbors_medi = sumNeighbors_medi + A_pdmm_medi(i,j)*Z_medi(i,j);
        end
    % X Update equation
        if (-1 - sumNeighbors_medi)/ (c*d_medi) > a(i)
            x_pdmm_medi(i) = (-1 - sumNeighbors_medi)/ (c*d_medi);
        elseif (1 - sumNeighbors_medi)/ (c*d_medi) < a(i)
            x_pdmm_medi(i) = (1 - sumNeighbors_medi)/ (c*d_medi);
        else
            x_pdmm_medi(i) = a(i);%(0 - sumNeighbors_medi)/ (c*d_medi);
        end
        
        for index = pdmm_neighbors'
            j = neighbors_pdmm(index, 2);
            Y_medi(i,j) = Z_medi(i,j) + 2*c*A_pdmm_medi(i,j)*x_pdmm_medi(i); %Y update equation
            %Z(j,i) = Y(i,j); %Z update equation
            transmissions_pdmm_medi = transmissions_pdmm_medi + 1;
        end    
    end   
    for i = 1:n
        pdmm_neighbors_medi = find(neighbors_pdmm(:,1) == i);
        for index = pdmm_neighbors_medi'
            j = neighbors_pdmm(index, 2);
            Z_medi(j,i) = (1/2)*(Y_medi(i,j) + Z_medi(j,i)) ;
        end
    end    
end

%error_pdmm_medi(k + 1,1) = norm(x_pdmm_medi - medianBase_pdmm,2)^2/n;

% figure(4)
% hold on
% grid on
% plot(error_pdmm_medi)
% legend('PDMM')
% set(gca, 'YScale', 'log')
% ylim([10e-15 10e5])
% hold off
fprintf('Final error: %f \n', error_pdmm_medi(end))



%% Plotting everything
% Plot error of mean consensus problem
transmissions = linspace(0, transmissions_gossip, transmissions_gossip/2 + 1);
transmissions_pdmm = linspace(0, transmission_pdmm_mean, 100);

figure(3)
clf(3)
hold on
semilogy(transmissions, error, 'b', 'LineWidth', 1.5)
semilogy(transmissions_pdmm, error_pdmm, 'r', 'LineWidth', 1.5)


legend({'Randomised gossip', ...
        'PDMM'}, ...
       'Location', 'northeast', ...
       'FontSize', 10);


grid on
title(['Mean Consensus problem Nodes = ' , num2str(n)])
xlabel('transmissions', 'FontSize', 12);
ylabel('$\|x^{(k)} - x^*\|$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log');
ylim([1e-10 1e2]);
hold off

% Plot error median consensu problem

transmissions_pdmm_medi_plot = linspace(0, transmissions_pdmm_medi, 10000);

figure(4)
clf(4)
hold on
%semilogy(transmissions, error, 'b', 'LineWidth', 1.5)
semilogy(transmissions_pdmm_medi_plot, error_pdmm_medi, 'b', 'LineWidth', 1.5)


legend({'PDMM'}, ...
       'Location', 'northeast', ...
       'FontSize', 10);


grid on
title(['Median Consensus problem Nodes = ' , num2str(n)])
xlabel('transmissions', 'FontSize', 12);
ylabel('$\|x^{(k)} - x^*\|$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log');
ylim([1e-10 1e2]);
hold off