function [sq_grad, xminuxbar, fminufstar, time] = Prox_GPDA_beta(x_init, edge_index,A,B,d,n_agents,H, H_split,b, b_split, m,M,regulation, T, beta)

    fprintf('Starting Prox-GPDA\n');
    time = zeros(T, 1);
%     Opt_GPDA = zeros(T-1,1);
    sq_grad = zeros(T, 1); % T is the iteration number
    xminuxbar = zeros(T, 1);
    
    fminufstar = zeros(T, 1);
    
    rho_cost = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     xs = x0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    upd = textprogressbar(T);
    
   
     f1 = (1/beta)*inv( A'*A + B'*B );
     f2 = beta*B'*B;
     
     f1=sparse(f1);
     f2=sparse(f2);
     A=sparse(A);
     
    x = zeros(n_agents*d, T);
    x0_temp = reshape(x_init, [n_agents*d, 1]);
    x(:, 1) = x0_temp;
     
    mu = zeros((edge_index-1)*d,T);
    
    for iter  = 2 : T+1
        
        tic;
        % calculate the gradient
        gradient = zeros(n_agents*d,1);        
        xs_ = reshape(x(:, iter-1), [d, n_agents]);

        for ii = 1 : n_agents
            H_temp = H_split{ii};
            b_temp = b_split{ii};
            gradient( (ii-1)*d+1:ii*d ) =  grad_loss(H_temp, xs_(:, ii), b_temp, rho_cost, m, regulation );
        end
       
        
        %% update x and mu
        
        x(:,iter) =  f1 * (f2 * x(:,iter-1) - gradient - A'*mu(:,iter-1));
        mu(:,iter) = mu(:,iter-1) + beta*  A *x(:,iter);


        %% Calculate the things that can measure the performance
        x_ = reshape( x(:, iter), [d, n_agents] );
        x_avg = sum(x_, 2)/n_agents;

        for k = 1: n_agents
            xminuxbar(iter-1) = xminuxbar(iter-1)+(norm(x_(:, k)-x_avg))^2;
        end
        tmp_grad = grad_loss(H, x_avg, b, rho_cost, M, regulation);
        sq_grad(iter-1) = sum(tmp_grad.^2);
        
        upd(iter);
        t_temp = toc;
        time(iter) = time(iter-1) + t_temp;
    end
end
