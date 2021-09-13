clear; clc; close all;

load('opt_data.mat');

T = 500;
beta = 5;
d = n;
n_agents = N;

[sq_grad_Prox_PDA, xminuxbar_Prox_PDA, ~, time_Prox_PDA] =  Prox_PDA_beta(x_init, edge_index,A,B,d,n_agents,H, H_split,b, b_split, m,M,regulation, T, beta);
    
% save ('Prox_PDA_grad.mat', 'sq_grad_Prox_PDA');
% save ('Prox_PDA_x.mat', 'xminuxbar_Prox_PDA');

%% Plot

iterations = T;
Prox_PDA_plot = sq_grad_Prox_PDA/n_agents^2+xminuxbar_Prox_PDA/n_agents;
figure(7);
k = 1:iterations;
plot(k, Prox_PDA_plot, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
