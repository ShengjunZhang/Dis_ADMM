function [sq_grad, xminuxbar, fminufstar, time] = Prox_PDA_beta(x0, edge_index,A,B,d,n_agents,y_all, a_Re_all, a_Im_all, y, a_Re, a_Im,T, beta)

    fprintf('Starting Prox-PDA\n');
    time = zeros(T, 1);
%     Opt_GPDA = zeros(T-1,1);
    sq_grad = zeros(T, 1); % T is the iteration number
    xminuxbar = zeros(T, 1);
    
    fminufstar = zeros(T, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xs = x0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    upd = textprogressbar(T);
    
     f1 = (1/beta)*inv( A'*A + B'*B );
     f3 = beta*( A'*A + B'*B );
     f2 = beta*B'*B;
     f1=sparse(f1);
     f3=sparse(f3);
     f2=sparse(f2);
     A=sparse(A);
     
    x = zeros(n_agents*d, T);
    x0_temp = reshape(x0, [n_agents*d, 1]);
    x(:, 1) = x0_temp;
     
    mu = zeros((edge_index-1)*d,T);
    options = optimset('Display','off');
    for iter  = 2 : T+1
        tic;
        % calculate the gradient
        gradient = zeros(n_agents*d,1);        
        xs_ = reshape(x(:, iter-1), [d, n_agents]);

        for ii = 1 : n_agents
            gradient( (ii-1)*d+1:ii*d ) =  grad_loss( xs_(:, ii), y{ii}, a_Re{ii}, a_Im{ii} );
        end
        
        
        
        %% update x and mu
        pre_x=f2 * x(:,iter-1);
        pre_mu=A'*mu(:,iter-1);
        x_in =  f1 * ( pre_x- gradient - pre_mu);
        
        for k = 1:n_agents
            pre_f3=f3( (k-1)*d+1:k*d,(k-1)*d+1:k*d);
            initial = x_in((k-1)*d+1:k*d);
% %             initial = x((k-1)*d+1:k*d,iter-1);
% %             initial = zeros(d, 1);
            x((k-1)*d+1:k*d,iter) = fsolve(@(x) Prox_PDA_fun(x,y{k},a_Re{k},a_Im{k},pre_mu((k-1)*d+1:k*d), pre_x((k-1)*d+1:k*d),...
                                                     pre_f3),initial,options);
        end
        
        mu(:,iter) = mu(:,iter-1) + beta*  A *x(:,iter);

        %% Calculate the things that can measure the performance
        x_ = reshape( x(:, iter), [d, n_agents] );
        x_avg = sum(x_, 2)/n_agents;

        for k = 1: n_agents
            xminuxbar(iter-1) = xminuxbar(iter-1)+(norm(x_(:, k)-x_avg))^2;
        end
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad(iter-1) = sum(tmp_grad.^2);
        
        upd(iter);
        t_temp = toc;
        time(iter) = time(iter-1) + t_temp;
        
    end
end
