function [y, a_Re, a_Im] = measure_magnitude(x, m, sig)
d = length(x);
a_Re = randn(d, m) / sqrt(2);
a_Im = randn(d, m) / sqrt(2);

y = abs((a_Re + 1j * a_Im)' * x) + randn(m, 1) * sig;
end