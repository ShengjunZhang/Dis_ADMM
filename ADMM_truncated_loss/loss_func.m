function [fval] = loss_func(H, x, b, rho_cost, m, regulation)
    
    H = H';
    x = x';
    [~,n] = size(H);
    pred = x*H;
    obj =   (0.5*rho_cost).*log(1+(pred - b).^2./rho_cost) + regulation * norm(x, 1);
    fval =  obj / n;
end