%% Proposed Linearized ADMM Algorithm

clear; clc; close all;

load('opt_data.mat');


T = 500;
d = n;
n_agents = N;

rho_cost = 3;
% L_aug=A'*A;
sq_grad_linearized_admm = zeros(T, 1);
xminuxbar_linearized_admm = zeros(T, 1);

xs = x_init(:, :);
gs = zeros(d, n_agents);
vs = zeros(d, n_agents);

% tmp_grad = grad_loss(H, sum(xs, 2) / n_agents, b, rho_cost, M, regulation);

gamma1 = 4000000; % 37.6;
eta=1/gamma1;
alpha = 1000000;%10
beta = 20000;%8
% upd = textprogressbar(T);

for t = 1:T
    fprintf('iteration %d.\n', t);
    % Update x
    temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
    temp_xs = reshape(temp_xs, [d, n_agents]);
    
    for k = 1:n_agents
        H_temp = H_split{k}; 
        b_temp = b_split{k};
        gs(:, k) = grad_loss(H_temp, xs(:, k), b_temp, rho_cost, m, regulation);
    end
    
    
    xs = xs - eta * ( alpha * temp_xs + beta*vs + gs);
    
    % broadcast and receive.
    temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
    temp_xs = reshape(temp_xs, [d, n_agents]);
    
    % Update v
    temp_vs = temp_xs;
    vs = vs + eta * beta *temp_vs;
    
    x_avg = sum(xs, 2) / n_agents;
    
    for k = 1:n_agents
        xminuxbar_linearized_admm(t) = xminuxbar_linearized_admm(t)+(norm(xs(:, k)-x_avg))^2;
        %             fminufstar_linearized_admm(t, p) = fminufstar_linearized_admm(t, p)+loss_func(x_avg, y{k}, a_Re{k}, a_Im{k});
    end
    
    tmp_grad = grad_loss(H, x_avg, b, rho_cost, M, regulation);
    sq_grad_linearized_admm(t) = sum(tmp_grad.^2);
    
end

% save ('linearized_admm_grad.mat', 'sq_grad_linearized_admm');
% save ('linearized_admm_x.mat', 'xminuxbar_linearized_admm');
% % save ('linearized_admm_f.mat', 'fminufstar_linearized_admm');

%% Plot
iterations = T;
plot_linearized_admm = xminuxbar_linearized_admm(1:iterations)/n_agents+sq_grad_linearized_admm(1:iterations)/n_agents^2;
semilogy(1: iterations, plot_linearized_admm, '-', 'LineWidth', 2);