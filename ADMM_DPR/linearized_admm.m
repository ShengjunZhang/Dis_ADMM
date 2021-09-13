%% Proposed Linearized ADMM Algorithm

% clear; 
% clc; 
%  close all;

% load('opt_data.mat');
n_inst=1;
T = 3000;
% L_aug=A'*A;
sq_grad_linearized_admm = zeros(T, n_inst);
xminuxbar_linearized_admm = zeros(T, n_inst);
% fminufstar_linearized_admm = zeros(T, n_inst);

for p = 1:n_inst
    xs = x0;
    gs = zeros(d, n_agents);
    vs = zeros(d, n_agents);
    tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);
    

    gamma1 = 37.6; % 37.6;
    eta=1/gamma1;
    alpha = 10;%10
    beta = 8;%8
    upd = textprogressbar(T);
    
    for t = 1:T
        


        
        % Update x  
        temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);%Add by Xinlei
        temp_xs = reshape(temp_xs, [d, n_agents]);%Add by Xinlei
        
        for k = 1:n_agents
            gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
        end            
        xs = xs - eta * ( alpha * temp_xs + beta*vs + gs);
        
        % broadcast and receive.
        temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
        temp_xs = reshape(temp_xs, [d, n_agents]);
        
        % Update v       
        temp_vs = temp_xs;
        vs = vs + eta * beta *temp_vs;
        
        x_avg = sum(xs, 2) / n_agents;
        
        for k = 1:n_agents
            xminuxbar_linearized_admm(t, p) = xminuxbar_linearized_admm(t, p)+(norm(xs(:, k)-x_avg))^2;
%             fminufstar_linearized_admm(t, p) = fminufstar_linearized_admm(t, p)+loss_func(x_avg, y{k}, a_Re{k}, a_Im{k});
        end
        
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad_linearized_admm(t, p) = sum(tmp_grad.^2);
        upd(t);
    end
    
    fprintf('case %d done\n',p);
end

% save ('linearized_admm_grad.mat', 'sq_grad_linearized_admm');
% save ('linearized_admm_x.mat', 'xminuxbar_linearized_admm');


% % Plot
case_choice = 1;
iterations = T;
linearized_admm_plot = sq_grad_linearized_admm(:, case_choice)+xminuxbar_linearized_admm(:, case_choice);
figure(3);
k = 1:iterations;
plot(k, linearized_admm_plot/n_agents, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
