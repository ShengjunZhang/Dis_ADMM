function [fval] = admm_fun(x, y, a_Re, a_Im, beta, v, pre_x, gamma, alpha, L, full)

    m = length(y);
    diff = y.^2 - ((a_Re'*x).^2+(a_Im'*x).^2);

    tmp_Re = bsxfun(@times, a_Re', diff);
    tmp_Im = bsxfun(@times, a_Im', diff);

    tmp_mat = - a_Re * tmp_Re - a_Im * tmp_Im;
    grad = tmp_mat * x / m + tmp_mat' * x / m;
    
    fval = grad + beta*v + gamma*x - gamma*pre_x + alpha*full*L;
end