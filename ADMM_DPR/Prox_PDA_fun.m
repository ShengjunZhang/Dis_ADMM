function [fval] = Prox_PDA_fun(x, y, a_Re, a_Im, mu, pre_x, f3)

    m = length(y);
    diff = y.^2 - ((a_Re'*x).^2+(a_Im'*x).^2);

    tmp_Re = bsxfun(@times, a_Re', diff);
    tmp_Im = bsxfun(@times, a_Im', diff);

    tmp_mat = - a_Re * tmp_Re - a_Im * tmp_Im;
    grad = tmp_mat * x / m + tmp_mat' * x / m;
    
    fval = grad +mu +f3*x -pre_x;
end