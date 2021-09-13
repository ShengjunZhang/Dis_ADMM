%% Tang Li Algorithm

% clear; clc; 
% close all;

load('opt_data.mat');
n_inst=1;
T = 3000;
% L_aug=A'*A;
sq_grad_qu_li = zeros(T, n_inst);
xminuxbar_qu_li = zeros(T, n_inst);
for p = 1:n_inst
    xs = x0;
    gs = zeros(d, n_agents);
    prev_gs = zeros(d, n_agents);
        for k = 1:n_agents
            prev_gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
        end
    ss = prev_gs;
    tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);
    
    eta1 = 0.03;%0.016,0.022
    upd = textprogressbar(T);
    for t = 1:T
        eta = eta1;
        xs = xs - eta * ss;
        xs = reshape(W_aug * reshape(xs, [d*n_agents, 1]), [d, n_agents]);
        
        
        for k = 1:n_agents
            gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
        end
        ss = ss + gs - prev_gs;
        ss = reshape(W_aug * reshape(ss, [d*n_agents, 1]), [d, n_agents]);
        
        prev_gs = gs;
        
        
        x_avg = sum(xs, 2) / n_agents;
        
        for k = 1:n_agents
            xminuxbar_qu_li(t, p) = xminuxbar_qu_li(t, p)+(norm(xs(:, k)-x_avg))^2;
       end
        
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad_qu_li(t, p) = sum(tmp_grad.^2);
        upd(t);
    end
    
    fprintf('case %d done\n',p);
end

% save ('qu_li_grad.mat', 'sq_grad_qu_li');
% save ('qu_li_x.mat', 'xminuxbar_qu_li');

% save tang_li.mat

% % Plot

case_choice = 1;
iterations = T;
tang_li_plot = sq_grad_qu_li(:, case_choice)+xminuxbar_qu_li(:, case_choice);
figure(2);
k = 1:iterations;
plot(k, tang_li_plot/n_agents, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');


