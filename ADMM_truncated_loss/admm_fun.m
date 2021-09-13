function [fval] = admm_fun(H, x, b, rho_cost,m, regulation, beta, v, pre_x, gamma, alpha, L, full)
    
    H = H';
    x = x';
    [~,n] = size(H);
    temp = x*H-b;
    grad = rho_cost*(temp./(n+temp.^2./(rho_cost/n)))*H' + regulation * sign(x);
    grad = sum(grad', 2)/n;
    
    fval = grad + beta*v + gamma*x - gamma*pre_x + alpha*full*L;
%     fval = (0.5*rho_cost/m) .* log( 1+(b - H * x).^2./ rho_cost ) + regulation * norm(x, 1);
end

% @(H, x, b, rho_cost, m, regulation) (0.5*rho_cost/m) .* log( 1+(b - H * x).^2./ rho_cost ) + regulation * norm(w, 1);

