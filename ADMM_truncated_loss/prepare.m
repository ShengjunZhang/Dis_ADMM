%% This is used for generating data and graph

clear; clc; close all;

%% Data

rng('default');
s = rng;

m = 150;        % number of examples per agnet
n = 256;        % number of features per data
N = 20;         % number of agents
M = N * m;      % total number of data samples

% % Get the optimal x

x_star =  2. * sprand(n, 1, 0.065);
S = nnz(x_star);
% fprintf('There are %d non-zero elements in optimal x\n', S);
row = find(x_star ~= 0);
col = ones(S, 1);
offset = ones(S, 1);
Offset = sparse(row, col, offset, n, 1);

x_optimal = x_star - Offset;


H_split = {};
b_split = {};
H = randn(M, n);
H = H * spdiags(1./sqrt(sum(H.^2))', 0, n, n);  % normalize columns
v = sqrt(4) * randn(M, 1);                      % noise
b = H * x_optimal + v;

for i = 1:N
    H_split{i} = H((i-1)*150 + 1: i*150, :);
    b_split{i} = b((i-1)*150 + 1: i*150, :);
end

regulation_max = norm(H' * b, 'inf');
regulation = 0.1 * regulation_max;
x_init = 2 * rand(n, N) - 1;


%% Graph

[Lap, W, Adj] = RingLaplacian(N);

W_aug = sparse(kron(W, eye(n)));
L_aug = sparse(kron(Lap, eye(n)));

[Adj,degree,num_of_edge,A,B,D,Lm,edge_index,eig_Lm,min_eig_Lm,WW,LN,L_hat,eig_L_hat,min_eig_L_hat] = graph_transform(Adj, N, n);

%% Cost function

% rho_cost = 3;
% 
% cf = @(H, x, b, rho_cost, m, regulation) (0.5*rho_cost/m) .* log( 1+(b - H * x).^2./ rho_cost ) + regulation * norm(w, 1);
% gf = @(H, x, b, rho_cost, m, regulation) (1/m) * ( (b - H * x)./((1+(b - H * x)^2/rho_cost).*H) )' + regulation * sign(x);

%% Save data and graph

save opt_data.mat;