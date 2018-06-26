function [ nextX ] = HermitianMatrix( A, b, E, eps )
disp('Метод простой итерации');
disp('2. Эрмитова положительно определенная матрица:');

% Максимальное по модулю собственное число А степенным методом
maxEigenvalue=PowerMethod(A, eps);
disp('      Максимальное по модулю собственное число А=');
disp(maxEigenvalue);

% Максимальное по модулю собственное число А методом скалярного произведения
maxEigenvalue = ScalarProduct(A, eps);
disp('      Максимальное по модулю собственное число А=');
disp(maxEigenvalue);

% Минимальное по модулю собственное число А методом скалярного произведения
B = A - maxEigenvalue * E;
minEigenvalue=maxEigenvalue +ScalarProduct(B, eps);
disp('      Минимальное по модулю собственное число А=');
disp(minEigenvalue);

alfa = 2/(minEigenvalue+maxEigenvalue);
H = E - alfa*A;
g = alfa * b;

prevX = [1; 0; 0];
nextX = H * prevX + g;

count = 0;

while norm(nextX - prevX) > eps
    prevX = nextX;
    nextX = H * prevX + g;
    count = count + 1;
end; 

disp('Погрешность:');
disp(b - A * nextX);
disp('Количество шагов:');
disp(count);

end

