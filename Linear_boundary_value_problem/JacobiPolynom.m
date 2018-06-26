% Function of computing the Jacobi polynomial:

function [Pn] = JacobiPolynom(n)
syms x;

if n == 0
  Pn = 1;
else
    if n == 1
        Pn = 2 * x;
    else
        Pn = ((n + 1) * (2 * (n - 2) + 5) * x * JacobiPolynom(n - 1) - (n + 1) * n * JacobiPolynom(n - 2)) / ((n + 2) * n);
    end
end

end
