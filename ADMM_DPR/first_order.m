%% Proposed first order Algorithm

%clear; clc; 
% close all;

load('opt_data.mat');
n_inst=1;
T = 3000;
% L_aug=A'*A;
sq_grad_first_order = zeros(T, n_inst);
xminuxbar_first_order = zeros(T, n_inst);

for p = 1:n_inst
    xs = x0;
    gs = zeros(d, n_agents);
    vs = zeros(d, n_agents);
    tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);
    
    eta1 = 0.029;%0.029,0.012
    alpha = 10;%10,8
    beta = 8;%8,6
    upd = textprogressbar(T);
    for t = 1:T
        eta = eta1;


        for k = 1:n_agents
            gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
        end
        
        temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
        temp_xs = reshape(temp_xs, [d, n_agents]);
        
        xs = xs - eta * ( alpha * temp_xs + beta*vs + gs);
        
        temp_vs = temp_xs;
        vs = vs + eta * beta *temp_vs;
        x_avg = sum(xs, 2) / n_agents;
        for k = 1:n_agents
            xminuxbar_first_order(t, p) = xminuxbar_first_order(t, p)+(norm(xs(:, k)-x_avg))^2;
        end
        
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad_first_order(t, p) = sum(tmp_grad.^2);
        upd(t);
    end
    
    fprintf('case %d done\n',p);
end

% save ('first_order_grad.mat', 'sq_grad_first_order');
% save ('first_order_x.mat', 'xminuxbar_first_order');


% % Plot
case_choice = 1;
iterations = T;
first_order_plot = sq_grad_first_order(:, case_choice)+xminuxbar_first_order(:, case_choice);
figure(8);
k = 1:iterations;
plot(k, first_order_plot/n_agents, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
