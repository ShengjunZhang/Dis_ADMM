%% Plot ADMMs and Qu Li

%    clear; clc; 
%   close all;

load ('Prox_GPDA_grad.mat');
load ('Prox_GPDA_x.mat');

load ('Prox_PDA_grad.mat');
load ('Prox_PDA_x.mat');


load ('non_linearized_admm_grad.mat');
load ('non_linearized_admm_x.mat');


load ('linearized_admm_grad.mat');
load ('linearized_admm_x.mat');


load ('first_order_grad.mat');
load ('first_order_x.mat');


load ('qu_li_grad.mat');
load ('qu_li_x.mat');


n_agents=50;
case_choice = 1; % % this could be choose from 1-6
iterations = 400;

plot_Prox_GPDA = xminuxbar_Prox_GPDA(1:iterations, case_choice)/n_agents+sq_grad_Prox_GPDA(1:iterations, case_choice)/n_agents^2;
plot_Prox_PDA = xminuxbar_Prox_PDA(1:iterations, case_choice)/n_agents+sq_grad_Prox_PDA(1:iterations, case_choice)/n_agents^2;
plot_non_linearized_admm = xminuxbar_non_linearized_admm(1:iterations, case_choice)/n_agents+sq_grad_non_linearized_admm(1:iterations, case_choice)/n_agents^2;
plot_linearized_admm = xminuxbar_linearized_admm(1:iterations, case_choice)/n_agents+sq_grad_linearized_admm(1:iterations, case_choice)/n_agents^2;
plot_first_order = xminuxbar_first_order(1:iterations, case_choice)/n_agents+sq_grad_first_order(1:iterations, case_choice)/n_agents^2;
plot_qu_li = xminuxbar_qu_li(1:iterations, case_choice)/n_agents+sq_grad_qu_li(1:iterations, case_choice)/n_agents^2;
for i=2:iterations
    plot_Prox_GPDA(i)=min(plot_Prox_GPDA(i-1:i));
    plot_Prox_PDA(i)=min(plot_Prox_PDA(i-1:i));
    plot_non_linearized_admm(i)=min(plot_non_linearized_admm(i-1:i));
    plot_linearized_admm(i)=min(plot_linearized_admm(i-1:i));
    plot_first_order(i)=min(plot_first_order(i-1:i));
    plot_qu_li(i)=min(plot_qu_li(i-1:i));
end

% Prox_GPDA_plotgrad = sq_grad_Prox_GPDA(1:iterations, case_choice);
% non_linearized_admm_plotgrad = sq_grad_non_linearized_admm(1:iterations, case_choice);
% linearized_admm_plotgrad = sq_grad_linearized_admm(1:iterations, case_choice);
% first_order_plotgrad = sq_grad_first_order(1:iterations, case_choice);
% qu_li_plotgrad = sq_grad_qu_li(1:iterations, case_choice);




figure(1);

k = 1:iterations;


plot(k, plot_non_linearized_admm, '-', 'LineWidth', 2); hold on;
plot(k, plot_linearized_admm, '-.', 'LineWidth', 2); hold on;
% plot(k, plot_first_order, ':', 'LineWidth', 2); hold on;
plot(k, plot_qu_li, '--', 'LineWidth', 2); hold on;
plot(k, plot_Prox_PDA, '-', 'LineWidth', 2); hold on;
plot(k, plot_Prox_GPDA, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('Number of iterations $T$','Interpreter', 'latex', ...
        'FontSize', 15, 'FontWeight','bold');
ylabel('$\min_{k\in [T]}\{\|\nabla f(\bar{x}_k)\|^ {2}+\frac{1}{n}\sum_{i=1}^{n}\|x_{i,k}-\bar{x}_k\|^ {2}\}$', 'Interpreter','latex', ...
        'FontSize', 15, 'FontWeight','bold');
xticks(0:50:iterations);

legend({'Algorithm 1',...
        'Algorithm 2', ...
       'DGTA',...
       'Prox-PDA', ...
       'Prox-GPDA',}, ...
       'Interpreter', 'latex', 'FontSize', 10, 'FontWeight','bold');   
   
% figure(2);
% 
% k = 1:iterations;
% 
% plot(k, admm_plotgrad, '-', 'LineWidth', 2); hold on;
% plot(k, non_linearized_admm_plotgrad, '-', 'LineWidth', 2); hold on;
% plot(k, linearized_admm_plotgrad, '-.', 'LineWidth', 2); hold on;
% plot(k, first_order_plotgrad, ':', 'LineWidth', 2); hold on;
% plot(k, qu_li_plotgrad, '--', 'LineWidth', 2); hold on;
% set(gca,'FontSize', 10);
% set(gca, 'YScale', 'log');
% xlabel('iteration $k$','Interpreter', 'latex', ...
%         'FontSize', 15, 'FontWeight','bold');
% ylabel('$\| \nabla f (\bar{x})\|^ {2}$', 'Interpreter','latex', ...
%         'FontSize', 15, 'FontWeight','bold');
% xticks(0:200:iterations);
% 
% legend({'Classic ADMM',...
%         'Proposed Non-Linearized ADMM',...
%         'Proposed Linearized ADMM with $\gamma = 40$', ...
%         'Proposed first-order algorithm', ...
%        'First-order algorithm in Qu Li',}, ...
%        'Interpreter', 'latex', 'FontSize', 10, 'FontWeight','bold');
   
   
