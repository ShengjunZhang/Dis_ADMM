clear; close all; clc;
iterations = 500;

load('plot_non_linearized_admm.mat');
load('plot_linearized_admm.mat');
load('plot_qu_li.mat');
load('plot_prox_gpda.mat');
load('plot_prox_pda.mat');

%%
semilogy(1: iterations, plot_non_linearized_admm, '-', 'LineWidth', 2); hold on;
semilogy(1: iterations, plot_linearized_admm, '-.', 'LineWidth', 2);
semilogy(1: iterations, plot_qu_li, '--', 'LineWidth', 2);
semilogy(1: iterations, Prox_PDA_plot, '-', 'LineWidth', 2);
semilogy(1: iterations, Prox_GPDA_plot, '-', 'LineWidth', 2);

xlabel('Number of iterations $T$','Interpreter', 'latex', ...
        'FontSize', 15, 'FontWeight','bold');
ylabel('$\min_{k\in [T]}\{\|\nabla f(\bar{x}_k)\|^ {2}+\frac{1}{n}\sum_{i=1}^{n}\|x_{i,k}-\bar{x}_k\|^ {2}\}$', 'Interpreter','latex', ...
        'FontSize', 15, 'FontWeight','bold');
    
xticks(0:100:iterations);

legend({'Algorithm 1',...
        'Algorithm 2',...
        'DGTA',...
        'Prox-PDA', ...
        'Prox-GPDA'}, ...
       'Interpreter', 'latex', 'FontSize', 10, 'FontWeight','bold');   