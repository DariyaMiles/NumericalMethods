function [ nextX ] = DiagonalPreponderance( A, b, N, E, eps )
disp('Метод простой итерации');
disp('1. Матрица с диагональным преобладанием:');

D = zeros(N);

for i = 1:1:N
    D(i, i) = A(i, i);
end;

g = D\b;
H = E - D\A;

prevX = [1; 0; 0];
nextX = H * prevX + g;

count = 0;

disp('||H||=');
disp(norm(H));

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

