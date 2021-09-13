function [fval] = loss_func(x, y, a_Re, a_Im)
    m = length(y);
    fval = 1/(2*m)*sum((y.^2 - ((a_Re'*x).^2+(a_Im'*x).^2)).^2);
end