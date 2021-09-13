%% Proposed nonlinearized ADMM Algorithm

clear; clc; 
% close all;

load('opt_data.mat');

n_inst=1;
T = 3000;
it_min=1000;
err=0.000001;

sq_grad_admmgac = zeros(T, n_inst);
xminuxbar_admmgac = zeros(T, n_inst);
% fminufstar_admmgac = zeros(T, n_inst);
options = optimset('Display','off');
for p = 1:n_inst
    xs = x0;
    gs = zeros(d, n_agents);
    vs = zeros(d, n_agents);
    tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);
    
    gamma1 = 80; % 33.3, 50;
    gamma2=37.6;
    theta=0.005;
    alpha = 10;
    beta = 8;
    upd = textprogressbar(T);
    for t = 1:T
        

        gamma =gamma2+ (gamma1-gamma2)*t/T;
        
%         Update x
% %         temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);%Add by Xinlei
% %         temp_xs = reshape(temp_xs, [d, n_agents]);%Add by Xinlei
% %         
% %         for k = 1:n_agents
% %             gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
% %         end            
% %         x_in = xs - (1/gamma) * ( alpha * temp_xs + beta*vs + gs);
%         
%         pre_x = xs;
%         for k = 1:n_agents
% %             initial = zeros(d, 1);
%             initial = pre_x(:, k);
% %             initial = x_in(:, k);
%             xs(:, k) = fsolve(@(x) admm_fun(x,y{k},a_Re{k},a_Im{k},beta,vs(:, k), pre_x(:,k),...
%                                                      gamma, alpha,L_small(:, k), xs),initial,options);
%         end
%         
        
        temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);%Add by Xinlei
        temp_xs = reshape(temp_xs, [d, n_agents]);%Add by Xinlei
        
        for k = 1:n_agents
            gs(:, k) = grad_loss(xs(:, k), y{k}, a_Re{k}, a_Im{k});
        end            
        xs_min = xs - theta * ( alpha * temp_xs + beta*vs + gs);
        xs_min2 = xs;
        for k = 1:n_agents
        tmin=2;
        while tmin<=it_min&&norm(xs_min2(:, k)-xs_min(:, k))>=err 
            xs_min2(:, k)=xs_min(:,k);
              gs(:, k) = grad_loss(xs_min2(:, k), y{k}, a_Re{k}, a_Im{k});              
            xs_min(:,k)=xs_min2(:, k)-theta*(alpha*temp_xs(:,k)+beta*vs(:,k)+gs(:,k)+gamma*(xs_min2(:,k)-xs(:,k)));
            tmin=tmin+1;
        end
        xs(:, k)=xs_min(:, k);
        end
        
        
        
        
        
        
        % broadcast and receive.
        temp_xs = L_aug * reshape(xs, [d*n_agents, 1]);
        temp_xs = reshape(temp_xs, [d, n_agents]);
        
        % Update v
        temp_vs = temp_xs;
        vs = vs + (1/gamma) * beta *temp_vs;
%         vs = vs + (1/beta) * alpha *temp_vs;
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        x_avg = sum(xs, 2) / n_agents;


        for k = 1:n_agents
            xminuxbar_admmgac(t, p) = xminuxbar_admmgac(t, p)+(norm(xs(:, k)-x_avg))^2;
%             fminufstar_admmgac(t, p) = fminufstar_admmgac(t, p)+loss_func(x_avg, y{k}, a_Re{k}, a_Im{k});
        end
        
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad_admmgac(t, p) = sum(tmp_grad.^2);
        upd(t);
    end
    
    fprintf('case %d done\n',p);
end

% save ('admmgac_grad.mat', 'sq_grad_admmgac');
% save ('admmgac_x.mat', 'xminuxbar_admmgac');


%Plot
case_choice = 1;
iterations = T;
admm_plot = sq_grad_admmgac(:, case_choice);
figure(4);
k = 1:iterations;
plot(k, admm_plot/n_agents, '-', 'LineWidth', 2); hold on;
set(gca,'FontSize', 10);
set(gca, 'YScale', 'log');
xlabel('iteration $k$','Interpreter', 'latex', 'FontSize', 15, 'FontWeight','bold');
ylabel('$|| \nabla f (\bar{x})||^ {2}$', 'Interpreter','latex','FontSize', 15, 'FontWeight','bold');
