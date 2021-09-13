% clear; clc; 
% close all;

 load('opt_data.mat');
% % 
% % x_true = [1; zeros(d-1, 1)];
% % 
% % y = cell(n_agents, 1);
% % a_Re = cell(n_agents, 1);
% % a_Im = cell(n_agents, 1);
% % 
% % y_all = zeros(m * n_agents, 1);
% % a_Re_all = zeros(d, m * n_agents);
% % a_Im_all = zeros(d, m * n_agents);
% % for k = 1:n_agents
% %     [y{k}, a_Re{k}, a_Im{k}] = measure_magnitude(x_true, m, 1);
% %     y_all((k-1)*m+1:k*m) = y{k};
% %     a_Re_all(:, (k-1)*m+1:k*m) = a_Re{k};
% %     a_Im_all(:, (k-1)*m+1:k*m) = a_Im{k};
% % end

T =3000;
beta=1;
[sq_grad_Prox_PDA, xminuxbar_Prox_PDA, ~, time_Prox_PDA] =  Prox_PDA_beta(x0, edge_index,A,B,d,n_agents,y_all, a_Re_all, a_Im_all, y, a_Re, a_Im,T, beta);
    
% save ('Prox_PDA_grad.mat', 'sq_grad_Prox_PDA');
% save ('Prox_PDA_x.mat', 'xminuxbar_Prox_PDA');
% % Plot

iterations = T;
Prox_PDA_plot = sq_grad_Prox_PDA/n_agents^2+xminuxbar_Prox_PDA/n_agents;
figure(7);
k = 1:iterations;
plot(k, Prox_PDA_plot, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
