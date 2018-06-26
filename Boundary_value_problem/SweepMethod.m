function [m, k] = SweepMethod(N, h, p1, d1)

m = zeros(1, N);
k = zeros(1, N);

m(1) = p1;
k(1) = d1;

for i = 2 : 1 : N
    ai = 1 - h * p(i) / 2;
    bi = 1 + (p(i) * h) / 2;
    ci = 2 - h^2 * q(i);

    m(i )= bi / (ci - ai * m(i - 1));
    k(i) = ((ai * k(i - 1)) - (h^2 * f(i))) / (ci - (ai * m(i - 1)));
end;

end

