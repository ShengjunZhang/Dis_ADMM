function [fval] = Prox_PDA_fun(H, x, b, rho_cost,m, regulation, mu, pre_x, f3)

    H = H';
    x = x';
    [~,n] = size(H);
    temp = x*H-b;
    grad = rho_cost*(temp./(n+temp.^2./(rho_cost/n)))*H' + regulation * sign(x);
    grad = sum(grad', 2)/n;
    
    fval = grad +mu +f3.*x -pre_x;
end