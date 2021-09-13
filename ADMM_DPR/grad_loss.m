function [grad] = grad_loss(x, y, a_Re, a_Im)

    m = length(y);
    diff = y.^2 - ((a_Re'*x).^2+(a_Im'*x).^2);

    tmp_Re = bsxfun(@times, a_Re', diff);
    tmp_Im = bsxfun(@times, a_Im', diff);

    tmp_mat = - a_Re * tmp_Re - a_Im * tmp_Im;
    grad = tmp_mat * x / m + tmp_mat' * x / m;
    
end