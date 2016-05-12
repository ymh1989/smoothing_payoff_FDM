function x = thomas(alpha, beta, gamma, f)
n = length(f);
for i=2:n
    mult = alpha(i)/beta(i-1);
    beta(i) = beta(i) - mult* gamma(i-1);
    f(i) = f(i) - mult*f(i-1);
end
x(n) = f(n) / beta(n);
for i = n-1:-1:1
    x(i) = (f(i) - gamma(i) * x(i+1)) / beta(i);
end