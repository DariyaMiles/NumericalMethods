% Solution of a linear boundary value problem for a second-order ODE

syms x;
N = 10;
b = 1;
a = -1;

alpha = 3;
xi = zeros(1, N);
f = 2 * (3 * x^2 - alpha) / (x^2 + alpha)^3 + 2 * x / (alpha + 1);
solution_y = 1 / (x^2 + alpha) - 1 / (alpha + 1);

% Nods of the Chebyshev polynomial
for i = 1 : 1 : N
    xi(i) = cos(pi * (2 * i - 1) / (2 * N));
end;

disp('  Nods of the Chebyshev polynomial xi = ');
disp(xi);

% Coordinate functions
wi = sym(zeros(N, 1));
Lwi = sym(zeros(N, 1));

for i = 1 : 1 : N
    wi(i) = (1 - x^2) * JacobiPolynom(i - 1);
    Lwi(i) = diff(wi(i), 2) - (x^2 + alpha) * diff(wi(i)) - 2 * x * wi(i);
end;

disp('  wi =');
disp(wi);
disp('  Lwi = ');
disp(Lwi);

y = CollocationMethod(xi, wi, Lwi, f, N);

disp('  y(-0.5) =');
disp(vpa(subs(y, -0.5)));
disp('  Error at point "-0.5" = ');
disp(vpa(subs(y, -0.5) - subs(solution_y, -0.5)));

disp('  y(0)=');
disp(vpa(subs(y, 0)));
disp('  Error at point "0" = ');
disp(vpa(subs(y, 0) - subs(solution_y, 0)));

disp('  y(0.5) = ');
disp(vpa(subs(y, 0.5)));
disp('  Error at point "0.5" = ');
disp(vpa(subs(y, 0.5) - subs(solution_y, 0.5)));

y = GalerkinMethod(a, b, N, wi, Lwi, f);

disp('y(-0.5) = ');
disp(vpa(subs(y, -0.5)));
disp('Error at point "-0.5" = ');
disp(vpa(subs(y, -0.5) - subs(solution_y, -0.5)));

disp('y(0)=');
disp(vpa(subs(y, 0)));
disp('Error at point "0" = ');
disp(vpa(subs(y, 0) - subs(solution_y, 0)));

disp('y(0.5)=');
disp(vpa(subs(y, 0.5)));
disp('Error at point "0.5" = ');
disp(vpa(subs(y, 0.5) - subs(solution_y, 0.5)));
