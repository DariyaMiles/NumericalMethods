function [nextX] = HermitianMatrix(A, b, E, eps)
disp('Simple iteration method');
disp('2. Hermitian positive definite matrix:');

% Maximum modulus of an eigenvalue A by a power method
maxEigenvalue = PowerMethod(A, eps);
disp('      Maximum modulus of an eigenvalue A = ');
disp(maxEigenvalue);

% Maximum modulus of an eigenvalue A by scalar product method
maxEigenvalue = ScalarProduct(A, eps);
disp('      Maximum modulus of an eigenvalue A = ');
disp(maxEigenvalue);

% Minimum modulus of an eigenvalue A by scalar product method
B = A - maxEigenvalue * E;
minEigenvalue = maxEigenvalue + ScalarProduct(B, eps);
disp('      Minimum modulus of an eigenvalue A = ');
disp(minEigenvalue);

alfa = 2 / (minEigenvalue + maxEigenvalue);
H = E - alfa * A;
g = alfa * b;

prevX = [1; 0; 0];
nextX = H * prevX + g;

count = 0;

while norm(nextX - prevX) > eps
    prevX = nextX;
    nextX = H * prevX + g;
    count = count + 1;
end; 

disp('Error:');
disp(b - A * nextX);
disp('ÊNumber of steps:');
disp(count);

end

