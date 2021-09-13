%% This code is the proposed nonlinearized admm algorithms,
clear; clc; close all;

load('opt_data.mat');

T = 500;
n_inst = 1;
d = n;
n_agents = N;

sq_grad_non_linearized_admm = zeros(T, n_inst);
xminuxbar_non_linearized_admm = zeros(T, n_inst);
% fminufstar_non_linearized_admm = zeros(T, n_inst);

options = optimset('Display','off','Algorithm', 'levenberg-marquardt');

xs = x_init(:, :);
gs = zeros(d, n_agents);
prev_gs = zeros(d, n_agents);
vs = zeros(d, n_agents);
rho_cost = 3;

% tmp_grad = grad_loss(H, sum(xs, 2) / n_agents, b, rho_cost, M, regulation);
L_small = Lap;

%%
% eta1 = 0.1;
% gamma1 = 40;
% alpha = 10;
% beta = 10;
% u1 = 0.1;
gamma1 = 4000000; % 37.6;
eta=1/gamma1;
alpha = 1000000;%10
beta = 20000;%8
for t = 1:T
    fprintf('iteration %d.\n', t);
%     eta = eta1;
    gamma = gamma1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ADMM
    % Update x
    pre_x = xs;
    %         upd = textprogressbar(n_agents);
    for k = 1 : n_agents
        initial = xs(:, k);
        H_temp = H_split{k}; 
        b_temp = b_split{k};
        tic;
        xs(:, k) = fsolve(@(x) admm_fun(H_temp,x,b_temp,rho_cost,m,regulation,beta,vs(:, k), pre_x(:,k),...
            gamma, alpha, L_small(:, k), xs), initial, options);
        toc;
    end
    % (H, x, b, rho_cost, m, regulation)
    
    % broadcast and receive.
    temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
    temp_xs = reshape(temp_xs, [d, n_agents]);
    
    % Update v
    temp_vs = temp_xs;
    vs = vs + (1/gamma) * beta *temp_vs;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x_avg = sum(xs, 2) / n_agents;
    
    
    for k = 1:n_agents
        xminuxbar_non_linearized_admm(t) = xminuxbar_non_linearized_admm(t)+(norm(xs(:, k)-x_avg))^2;
%         fminufstar_non_linearized_admm(t) = fminufstar_non_linearized_admm(t)+loss_func(H, x_avg, b, rho_cost, m, regulation);
    end
    
    tmp_grad = grad_loss(H, x_avg, b, rho_cost, M, regulation);
    sq_grad_non_linearized_admm(t) = sum(tmp_grad.^2);

end
    
%     fprintf('case %d done\n',p);

% save ('non_linearized_admm_grad.mat', 'sq_grad_non_linearized_admm');
% save ('non_linearized_admm_x.mat', 'xminuxbar_non_linearized_admm');
% save ('non_linearized_admm_f.mat', 'fminufstar_non_linearized_admm');

%%
iterations = T;
plot_non_linearized_admm = xminuxbar_non_linearized_admm(1:iterations)/n_agents+sq_grad_non_linearized_admm(1:iterations)/n_agents^2;
plot(1: iterations, plot_non_linearized_admm, '-', 'LineWidth', 2);
