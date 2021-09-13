function [grad] = grad_loss(H, x, b, rho_cost, m, regulation)
    H = H';
    x = x';
    [~,n] = size(H);
    temp = x*H-b;
    grad = rho_cost*(temp./(n+temp.^2./(rho_cost/n)))*H' + regulation * sign(x);
    grad = sum(grad', 2)/n;
end