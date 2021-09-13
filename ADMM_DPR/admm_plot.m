clear; clc; close all;

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
for i = 2:iterations
    plot_Prox_GPDA(i)=min(plot_Prox_GPDA(i-1:i));
    plot_Prox_PDA(i)=min(plot_Prox_PDA(i-1:i));
    plot_non_linearized_admm(i)=min(plot_non_linearized_admm(i-1:i));
    plot_linearized_admm(i)=min(plot_linearized_admm(i-1:i));
    plot_first_order(i)=min(plot_first_order(i-1:i));
    plot_qu_li(i)=min(plot_qu_li(i-1:i));
end





figure(1);

k = 1:iterations;


plot(k, plot_non_linearized_admm, '-', 'LineWidth', 2); hold on;
plot(k, plot_linearized_admm, '-.', 'LineWidth', 2);
plot(k, plot_qu_li, '--', 'LineWidth', 2);
plot(k, plot_Prox_PDA, '-', 'LineWidth', 2);
plot(k, plot_Prox_GPDA, '-', 'LineWidth', 2);
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
   
   
