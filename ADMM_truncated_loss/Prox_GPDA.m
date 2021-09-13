clear; clc; 
close all;

load('opt_data.mat');

% % end

T = 500;
beta = 5;
d = n;
n_agents = N;

[sq_grad_Prox_GPDA, xminuxbar_Prox_GPDA, ~, time_Prox_GPDA] =  Prox_GPDA_beta(x_init, edge_index,A,B,d,n_agents,H, H_split,b, b_split, m,M,regulation, T, beta);
    
% save ('Prox_GPDA_grad.mat', 'sq_grad_Prox_GPDA');
% save ('Prox_GPDA_x.mat', 'xminuxbar_Prox_GPDA');
% % Plot

iterations = T;
Prox_GPDA_plot = sq_grad_Prox_GPDA/n_agents^2+xminuxbar_Prox_GPDA/n_agents;
figure(6);
k = 1:iterations;
plot(k, Prox_GPDA_plot, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
