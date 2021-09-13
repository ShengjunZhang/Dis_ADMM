%% Plot ADMMs and Qu Li

%    clear; clc; 
%   close all;

load ('admm_grad.mat');
load ('admm_x.mat');


load ('admmga_grad.mat');
load ('admmga_x.mat');


% load ('admmgac_grad.mat');
% load ('admmgac_x.mat');


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
iterations = 3000;

admm_plotx = xminuxbar_admm(1:iterations, case_choice);
admmga_plotx = xminuxbar_admmga(1:iterations, case_choice);
admmgac_plotx = xminuxbar_admmgac(1:iterations, case_choice);
% non_linearized_admm_plotx = xminuxbar_non_linearized_admm(1:iterations, case_choice);
linearized_admm_plotx = xminuxbar_linearized_admm(1:iterations, case_choice);
first_order_plotx = xminuxbar_first_order(1:iterations, case_choice);
qu_li_plotx = xminuxbar_qu_li(1:iterations, case_choice);


admm_plotgrad = sq_grad_admm(1:iterations, case_choice);
admmga_plotgrad = sq_grad_admmga(1:iterations, case_choice);
admmgac_plotgrad = sq_grad_admmgac(1:iterations, case_choice);
% non_linearized_admm_plotgrad = sq_grad_non_linearized_admm(1:iterations, case_choice);
linearized_admm_plotgrad = sq_grad_linearized_admm(1:iterations, case_choice);
first_order_plotgrad = sq_grad_first_order(1:iterations, case_choice);
qu_li_plotgrad = sq_grad_qu_li(1:iterations, case_choice);



figure(5);

k = 1:iterations;

plot(k, admm_plotx/n_agents+admm_plotgrad/n_agents, '-', 'LineWidth', 2); hold on;
plot(k, admmga_plotx/n_agents+admmga_plotgrad/n_agents, '-', 'LineWidth', 2); hold on;
plot(k, admmgac_plotx/n_agents+admmgac_plotgrad/n_agents, '-', 'LineWidth', 2); hold on;
plot(k, linearized_admm_plotx/n_agents+linearized_admm_plotgrad/n_agents, '-.', 'LineWidth', 2); hold on;
plot(k, first_order_plotx/n_agents+first_order_plotgrad/n_agents, ':', 'LineWidth', 2); hold on;
plot(k, qu_li_plotx/n_agents+qu_li_plotgrad/n_agents, '--', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', ...
        'FontSize', 15, 'FontWeight','bold');
ylabel('$\sum_{i=1}^{n}\|x_{i,k}-\bar{x}_k\|^ {2}$', 'Interpreter','latex', ...
        'FontSize', 15, 'FontWeight','bold');
xticks(0:500:iterations);

legend({'ADMM with $\gamma=37.6$',...
        'ADMM with $\gamma=18$',...
        'ADMMc with $\gamma=18$',...
        'L-ADMM', ...
        'DFO-PDA', ...
       'DFO-GT',}, ...
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
   
   
% figure(3);
% 
% k = 1:iterations;
% 
% plot(k, admm_plotf, '-', 'LineWidth', 2); hold on;
% plot(k, non_linearized_admm_plotf, '-', 'LineWidth', 2); hold on;
% plot(k, linearized_admm_plotf, '-.', 'LineWidth', 2); hold on;
% plot(k, first_order_plotf, ':', 'LineWidth', 2); hold on;
% plot(k, qu_li_plotf, '--', 'LineWidth', 2); hold on;
% set(gca,'FontSize', 10);
% set(gca, 'YScale', 'log');
% xlabel('iteration $k$','Interpreter', 'latex', ...
%         'FontSize', 15, 'FontWeight','bold');
% ylabel('$f(\bar{x}_k)-f^*$', 'Interpreter','latex', ...
%         'FontSize', 15, 'FontWeight','bold');
% xticks(0:200:iterations);
% 
% legend({'Classic ADMM',...
%         'Proposed Non-Linearized ADMM',...
%         'Proposed Linearized ADMM with $\gamma = 40$', ...
%         'Proposed first-order algorithm', ...
%        'First-order algorithm in Qu Li',}, ...
%        'Interpreter', 'latex', 'FontSize', 10, 'FontWeight','bold');