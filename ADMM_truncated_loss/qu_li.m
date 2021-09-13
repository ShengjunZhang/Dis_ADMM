%% Tang Li Algorithm

clear; clc; close all;

load('opt_data.mat');

T = 500;
d = n;
n_agents = N;
rho_cost = 3;

% L_aug=A'*A;
sq_grad_qu_li = zeros(T, 1);
xminuxbar_qu_li = zeros(T, 1);

xs = x_init(:, :);
gs = zeros(d, n_agents);
prev_gs = zeros(d, n_agents);

for k = 1:n_agents
    H_temp = H_split{k};
    b_temp = b_split{k};
    prev_gs(:, k) = grad_loss(H_temp, xs(:, k), b_temp, rho_cost, m, regulation);
end

% for k = 1:n_agents
%     prev_gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
% end

ss = prev_gs;
% tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);

eta1 = 1/40000000;
for t = 1:T
    fprintf('iteration %d.\n', t);
    eta = eta1;
    xs = xs - eta * ss;
    xs = reshape(W_aug * reshape(xs, [d*n_agents, 1]), [d, n_agents]);
    
    
    for k = 1:n_agents
        H_temp = H_split{k}; 
        b_temp = b_split{k};
        gs(:, k) = grad_loss(H_temp, xs(:, k), b_temp, rho_cost, m, regulation);
    end
    
    ss = ss + gs - prev_gs;
    ss = reshape(W_aug * reshape(ss, [d*n_agents, 1]), [d, n_agents]);
    
    prev_gs = gs;
    
    
    x_avg = sum(xs, 2) / n_agents;
    
    for k = 1:n_agents
        xminuxbar_qu_li(t) = xminuxbar_qu_li(t)+(norm(xs(:, k)-x_avg))^2;
    end
    
    tmp_grad = grad_loss(H, x_avg, b, rho_cost, M, regulation);
    sq_grad_qu_li(t) = sum(tmp_grad.^2);
end
    

% save ('qu_li_grad.mat', 'sq_grad_qu_li');
% save ('qu_li_x.mat', 'xminuxbar_qu_li');

% save tang_li.mat

% % Plot
%% Plot
iterations = T;
plot_qu_li = xminuxbar_qu_li(1:iterations)/n_agents+sq_grad_qu_li(1:iterations)/n_agents^2;
semilogy(1: iterations, plot_qu_li, '-', 'LineWidth', 2);

