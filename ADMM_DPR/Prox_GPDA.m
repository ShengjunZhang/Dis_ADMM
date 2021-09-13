function [sq_grad, xminuxbar, fminufstar, time] = Prox_GPDA(x0, edge_index,A,B,d,n_agents,y_all, a_Re_all, a_Im_all, y, a_Re, a_Im,T,beta)

    fprintf('Starting Prox-GPDA\n');
    time = zeros(T, 1);
%     Opt_GPDA = zeros(T-1,1);
    sq_grad = zeros(T, 1); % T is the iteration number
    xminuxbar = zeros(T, 1);
    
    fminufstar = zeros(T, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xs = x0;
    tmp_grad = grad_loss(sum(xs, 2) / n_agents, y_all, a_Re_all, a_Im_all);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    upd = textprogressbar(T);
    

     %beta = 9;
   
     f1 = (1/beta)*inv( A'*A + B'*B );
     f2 = beta*B'*B;%beta is missing
     
     
    x = zeros(n_agents*d, T);
    x0_temp = reshape(x0, [n_agents*d, 1]);
    x(:, 1) = x0_temp;
     
    mu = zeros((edge_index-1)*d,T);
    for iter  = 2 : T
        tic;
        % calculate the gradient
        gradient = zeros(n_agents*d,1);        
        xs_ = reshape(x(:, iter-1), [d, n_agents]);

        for ii = 1 : n_agents
            gradient( (ii-1)*d+1:ii*d ) = gradient( (ii-1)*d+1:ii*d ) + grad_loss( xs_(:, ii), y{ii}, a_Re{ii}, a_Im{ii} );
        end
        
        
        
        %% update x and mu
        
        x(:,iter) =  f1 * (f2 * x(:,iter-1) - gradient - A'*mu(:,iter-1));
        mu(:,iter) = mu(:,iter-1) + beta*  A *x(:,iter);


        %% Calculate the things that can measure the performance
        x_ = reshape( x(:, iter), [d, n_agents] );
        x_avg = sum(x_, 2)/n_agents;

        for k = 1: n_agents
            xminuxbar(iter) = xminuxbar(iter)+(norm(x_(:, k)-x_avg))^2;
        end
        tmp_grad = grad_loss(x_avg, y_all, a_Re_all, a_Im_all);
        sq_grad(iter) = sum(tmp_grad.^2);
        
        upd(iter);
        t_temp = toc;
        time(iter) = time(iter-1) + t_temp;
    end
end
