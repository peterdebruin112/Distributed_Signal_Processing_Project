clear all
close all

%% Create virtual area with randomly placed sensors
% Design a proper connected sensor network that covers the area of the plant by using a
% reasonable number of sensors. How many do we need to guarantee a connected sensor
% network? Motivate your choice.

% Given the maximal coverage distance between consecutive sensors, compute
% the minimum amount of sensors needed to guarantee a proper connected sensor
% network.
d = 2;                          % The dimensionality of the virtual area
area = 100;                     % Define the length and width of the area
prob = 0.9999;                  % Probability of connected graph
n = round(sqrt(1/(1-prob)));     % Compute amount of required nodes
sensorRadius = nthroot((2*log(n)...
                 /(n)),d)*area; % Compute minimum required connection radius

% Set random seed for repeatability
rng(5);

% Define the upper and lower bound of the
% measured values from the sensors in the network.
lower = -10;% Lower bound of sensor measurement
upper = 30; % Upper bound of sensor measurement

% Define the coordinates of the sensors and place them in the area.
coordinates = area*rand([n d]);                 % Define the locations of the sensors
measurment = lower + (upper - lower).*rand(n,1);% Define the measured value for each sensor
figure                                          % Visualise sensor placement in a scatter plot
scatter(coordinates(:, 1), coordinates(:, 2))

%% Construct a matrix defining the indeces of the neighbours of each node

% Initialise matrices
distanceSensors = zeros(n,n);
adjacency = zeros(n,n);
neighbors = [];

% Define for each node i the neighbors by checking if the distance between
% node i and any other node j is within the computed sensorRadius. If that
% is the case, then the edge is defined in the 2D neighbors matrix on a row.
for i = 1:n
    for j = 1:n
        % Compute the distance between selected node i and potential
        % neigbor node j and store the computed distance in a matrix.
        distanceSensors(i,j) = sqrt((coordinates(i,1) - coordinates(j,1))^2 ...
                        + (coordinates(i,2) - coordinates(j,2))^2);
        % Check whether distance is within the sensorRadius and i<j makes
        % sure no duplicate neighbor edges are added to the neighbors
        % matrix.
        if distanceSensors(i,j) <= sensorRadius
            adjacency(i,j) = 1;
            % Store only one time the node number of the selected node and
            % its neigbor
            if i<j
                neighbors = [neighbors ; i , j];
            end
        end
    end
end

% Visualise the connected network by adding the edges between the
% neighboring nodes.
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

%%%%%%%%%%%%%% RANDOMIZED GOSSIP %%%%%%%%%%%%%%%%%%%%
% Derive total number of edges and initialze matrices
numberEdges = size(neighbors,1);
W_ij = cell(numberEdges,1);
I = eye(n);

for k = 1:numberEdges
    i = neighbors(k,1);
    j = neighbors(k,2);

    % Standard basis vector for node i
    e_i = zeros(n,1); 
    e_i(i) = 1;

    % Standard basis vector for neighbor j
    e_j = zeros(n,1); 
    e_j(j) = 1;

    % Compute the random matrix W_ij
    W_ij{k} = I - 0.5 * (e_i - e_j) * (e_i - e_j)';
end

%% Compute P to find fastest averaging algorithm

inc = cell(n,1);
for k = 1:numberEdges
    inc{neighbors(k,1)} = [inc{neighbors(k,1)}, k];
    inc{neighbors(k,2)} = [inc{neighbors(k,2)}, k];
end

cvx_begin
    variable p(numberEdges) nonnegative     % Probabilities on edges
    variable t                              % Lambda_2

    % Declare W_bar to be an expression holder such that Matlab allows CVX
    % objects to be inserted into numeric arrays.
    % (https://cvxr.com/cvx/doc/basics.html)
    expression W_bar(n,n)

    % Compute expected value W_bar
    W_bar = 0;
    for k = 1:numberEdges
        W_bar = W_bar + p(k) * W_ij{k};
    end    
    W_bar = W_bar / n;
    
    % Define the objective function
    minimize(t)
    % minimize (lambda_max(W_bar - (1/n)*ones(n)))

    % Define the constraints
    % The following relaxation makes sure the second eigenvalue is
    % minimized.
    W_bar-(ones(n)/n)<=t*I;
    % Make sure the probabilities are all non-negative values.
    p >= 0;
    p <= 1;
   
    for i = 1:n
        % Check if the probabilities from all the edges from node i to all
        % neighbors sum to 1.
        sum(p(inc{i})) == 1;
    end    
cvx_end

P_opt = zeros(n);
for k = 1:numberEdges
    i = neighbors(k,1);
    j = neighbors(k,2);
    P_opt(i,j) = p(k);
    P_opt(j,i) = p(k);  % symmetric
end
%% Main Procedure Randomized Gossip
x = measurment;
K = 30000;
it_axis = 0:1:K;
% meanBase = ones([n 1])*(lower + upper)*(0.5);
meanBase = mean(measurment);
error_gossip = zeros([K+1 1]);

for k = 1:K
    % Compute the error
    error_gossip(k,1) = (norm(x - meanBase,2)/n)^2;

    % Pick a node i via a uniform distribution with probability p=1/n
    pickedNode = randi([1 n], 1);
    
    % Select the pickedNode and 'load' the probabilities of the
    % corresponding neibhbors of the pickedNode.
    edges_jdx = inc{pickedNode};
    prob_edges = p(edges_jdx);

    % Draw one of the edges according to the optimal probabilities of
    % selecting the different neighbors.
    edge_k = randsample(edges_jdx, 1, true, prob_edges);

    % Obtain the node j at the other end of the selected edge
    node_i = neighbors(edge_k,1);
    node_j = neighbors(edge_k,2);
    
    % Orient the pair of nodes
    if node_j == i
        node_j = node_i; 
        node_i = i; 
    end
    
    % Compute the average and update x in node i and neighbor j
    avg = (x(node_i) + x(node_j))/2;
    x(node_i) = avg;
    x(node_j) = avg;
end    
error_gossip(k+1,1) = (norm(x - meanBase,2)/n)^2;

%% Compute average through the PDMM

% Initialise the primal, dual and auxiliary variables
x_PDMM = measurment;
z = zeros(N);
error_PDMM = zeros([K+1 1]);
c = 1;

for k = 1:K
    %  Iterate through all the nodes in the network to update the x values
    for i = 1:n
        % Collect the indices (node numbers) of the neighbors of node i
        edges_jdx = inc{I};
    
        % Compute the sum Aij^T*zi|j of all the neighbors of node i. Since the
        % Aij is the adjacency matrix containing ones to select the nodes in z,
        % we can ommit this matrix to sum only all the zi|j's
        sum_Z = sum(z(i,edges_jdx));

        % Update the x values
        d = length(edges_jdx);
        x_PDMM(i) = (measurment(i) - sum_Z)/ (1 + c * d);

        % Update yi|j for all neighbors of node i
        for j = 1:d
            yi|j = z()
        end
    end
end

%% Plot converge of both randomized gossip and the PDMM algorithm

figure
plot(it_axis, error_gossip)
grid on
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('MSE')
% yscale('log')

fprintf('Final error: %f', error_gossip(end))
%% • Suppose the sensor network would like to compute the median of the measurement
% data. Implement the median consensus problem using the PDMM algorithm.

