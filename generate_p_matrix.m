function P = generate_p_matrix(adj)
% Generate the gossip probability matrix P for a given adjacency matrix
% adj: n x n symmetric adjacency matrix (0/1, no self-loops)
% P: n x n symmetric matrix, P(i,j) = probability of selecting (i,j)

n = size(adj,1);
P = zeros(n);

% Find all the indicdes of the edges (i < j)
edges = find(triu(adj,1));
num_edges = length(edges);

if num_edges == 0
    error('No edges in the network.');
end

% Assign equal probability to each edge
prob = 1/num_edges;

% Fill P for each edge with uniform probability 1/n
for idx = 1:num_edges
    [i, j] = ind2sub([n n], edges(idx));
    P(i,j) = prob;
    P(j,i) = prob; % symmetric
end

end