% clear; clc; close all;

load('opt_data.mat');

T = 3000;
beta_list = [2, 3, 4, 5, 6, 7];%[5, 9, 10, 11, 15, 20]
store = zeros(T,2,length(beta_list));

for i_beta = 1:length(beta_list)
    beta = beta_list(i_beta);
    fprintf('Beta = %d\n', beta);
    [sq_grad_Prox_GPDA, xminuxbar_Prox_GPDA, ~, time_Prox_GPDA] =  Prox_GPDA_beta(x0, edge_index,A,B,d,n_agents,y_all, a_Re_all, a_Im_all, y, a_Re, a_Im,T, beta);
    store(:,:,i_beta) = [sq_grad_Prox_GPDA, xminuxbar_Prox_GPDA];
end

iterations = T;
k = 1:iterations;
figure(100);
plotStyle = {'-','--','-.','-','--','-.'};

for i = 4:length(beta_list)
    plot(k, store(:,1,i)/n_agents^2 + store(1,:,i)/n_agents, plotStyle{i}, 'LineWidth', 2); 
    legendInfo{i} = ['\beta = ' num2str(beta_list(i))];
    hold on;
end
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
legend(legendInfo,'Location','Best');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$\|\nabla f(\bar{x}_k)\|^ {2}+\frac{1}{n}\sum_{i=1}^{n}\|x_{i,k}-\bar{x}_k\|^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
ylim([1e-35, 1e-0]);
xticks(0:500:iterations);
savefig(sprintf('Prox_GPDA_beta'));


figure(1);
plot(1:T, store(:,1,4)/n_agents^2 + store(1,:,4)/n_agents, '-', 'LineWidth', 2); hold on;
plot(1:T, store(:,1,5)/n_agents^2 + store(1,:,5)/n_agents, '-.', 'LineWidth', 2); hold on;
plot(1:T, store(:,1,6)/n_agents^2 + store(1,:,6)/n_agents, ':', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', ...
        'FontSize', 15, 'FontWeight','bold');
ylabel('$\|\nabla f(\bar{x}_k)\|^ {2}+\frac{1}{n}\sum_{i=1}^{n}\|x_{i,k}-\bar{x}_k\|^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
xticks(0:500:iterations);
legend({'$\beta=5$',...
        '$\beta=6$', ...
        '$\beta=7$',}, ...
       'Interpreter', 'latex', 'FontSize', 10, 'FontWeight','bold');   
