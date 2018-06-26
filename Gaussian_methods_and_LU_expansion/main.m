% Point method of solving a system of linear equations

A = [0.80464 * 10^-4 -0.50072 4.69392; 0.74464 -0.50072 1.92392; 1.02464*10^-4 1.09972 1.37392];
b = [4.95706; 2.42911; 2.92296];
n = 3;

% Gaussian method
disp('Gaussian method');
GaussianMethod(A, b, n);

% Modification of Gaussian method #1
disp('Modification of Gaussian method #1');
ModificationOfGaussianMethod(A, b, n);

% Modification of Gaussian method #2
disp('Modification of Gaussian method #2');
ModificationOfGaussianMethod2(A, b, n);

% Modification of Gaussian method #3
disp('Modification of Gaussian method #3');
ModificationOfGaussianMethod3(A, b, n);

% LU expansion
disp('LU expansion');
LUexpansion(A, b, n);
