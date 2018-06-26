% Numerical solution of the Fredholm integral equation of the 2nd kind

syms x;
syms y;
a = 0;
b = 1;
H = cos(x * y) / 2;
f = 1 + x - x^2;

% The method of replacing the kernel by a degenerate

% The expansion in a Taylor series into a sum of 3 terms:
Ht = taylor(H, x, a, 'Order', 5);
disp('The expansion in a Taylor series into a sum of 3 terms:');
disp(Ht);
alpha = [1 / 2, -x^2 / 4, x^4 / 48];
beta = [1, y^2, y^4];
U3 = ReplacingKernel(alpha, beta, f, a, b, 3);
disp('U(a) = ');
disp(vpa(subs(U3, a)));
disp('U((a + b) / 2) = ');
disp(vpa(subs(U3, (a + b) / 2)));
disp('U(b) = ');
disp(vpa(subs(U3, b)));

% Error
disp('Actual error ar point "a":');
disp(vpa(subs(U3 - int(H * subs(U3, y), y, a, b) - f, 0.01)));
disp('Actual error ar point "(a + b) / 2":');
disp(vpa(subs(U3-int(H*subs(U3, y), y, a, b)-f, (a+b)/2)));
disp('Actual error ar point "b":');
disp(vpa(subs(U3 - int(H * subs(U3, y), y, a, b) - f, b)));

% The expansion in a Taylor series into a sum of 4 terms:
Ht = taylor(H, x, a, 'Order', 7);
disp('The expansion in a Taylor series into a sum of 4 terms:');
disp(Ht);
alpha = [1 / 2, -x^2 / 4, x^4 / 48, -x^6 / 1440];
beta = [1, y^2, y^4, y^6];
U4 = ReplacingKernel(alpha, beta, f, a, b, 4);
disp('U(a) = ');
disp(vpa(subs(U4, a)));
disp('U((a + b) / 2) = ');
disp(vpa(subs(U4, (a + b) / 2)));
disp('U(b) = ');
disp(vpa(subs(U4, b)));

% Error
disp('Actual error ar point "a":');
disp(vpa(subs(U4 - int(H * subs(U4, y), y, a, b) - f, 0.01)));
disp('Actual error ar point "(a + b) / 2":');
disp(vpa(subs(U4 - int(H * subs(U4, y), y, a, b) - f, (a + b) / 2)));
disp('Actual error ar point "b":');
disp(vpa(subs(U4 - int(H * subs(U4, y), y, a, b) - f, b)));

% Mechanical quadratures method
N = 2;
disp('Mechanical quadratures method:');
Um = MechanicalQuadratures(H, f, a, b, N);
disp('N = ');
disp(N);
disp('U(a) = ');
disp(vpa(subs(Um, a)));
disp('U((a + b) / 2) = ');
disp(vpa(subs(Um, (a + b) / 2)));
disp('U(b) = ');
disp(vpa(subs(Um, b)));

N = 2 * N;
Un = MechanicalQuadratures(H, f, a ,b, N);
xk = [a; (a + b) / 2; b];
eps = 0.0001;

% Doubling nodes while the difference is greater than eps
while(max(abs(subs(Un, xk) - subs(Um, xk)), xk)) > eps
    N = 2 * N;
    Um = Un;
    Un = MechanicalQuadratures(H, f, a ,b, N);
end;

disp('ÊNumber of nodes N=');
disp(N);
disp('U(a) = ');
disp(vpa(subs(Un, a)));
disp('U((a + b) / 2) = ');
disp(vpa(subs(Un, (a + b) / 2)));
disp('U(b) = ');
disp(vpa(subs(Un, b)));

disp('Comparing the results with the kernel replacement method (3 terms):');
disp('at point "a":');
disp(abs(vpa(subs(Un, a)) - vpa(subs(U3, a))));
disp('at point "(a + b) / 2":');
disp(abs(vpa(subs(Un, (a + b) / 2)) - vpa(subs(U3, (a + b) / 2))));
disp('at point "b":');
disp(abs(vpa(subs(Un, b)) - vpa(subs(U3, b))));

disp('Comparing the results with the kernel replacement method (4 terms):');
disp('at point "a":');
disp(abs(vpa(subs(Un, a)) - vpa(subs(U4, a))));
disp('at point "(a + b) / 2":');
disp(abs(vpa(subs(Un, (a + b) / 2)) - vpa(subs(U4, (a + b) / 2))));
disp('at point "b":');
disp(abs(vpa(subs(Un, b)) - vpa(subs(U4, b))));
