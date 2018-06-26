% Iterative methods for solving SLAE

N = 3;
A = [2.16983 0.41261 0.23320; 0.41261 3.10770 0.61805; 0.23320 0.61805 4.76490];
b = [2.13221; 6.24929; 5.90552];
E = eye(N);
eps = 0.00001;

% Simple iteration method
% 1. Matrix with diagonal predominance
x = DiagonalPreponderance(A, b, N, E, eps);
disp('Solution:');
disp(x);

% Simple iteration method
% 2. Hermitian positive definite matrix
x = HermitianMatrix(A, b, E, eps);
disp('Solution:');
disp(x);

% Seidel method
x = SeidelMethod(A, b, N, E, eps);
disp('Solution:');
disp(x);



