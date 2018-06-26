function [ nextX] = SeidelMethod( A, b, N, E, eps )
disp('Метод Зейделя');

D = zeros(N);

for i = 1:1:N
    D(i, i) = A(i, i);
end;

g = D\b;
H = E - D\A;

prevX = [1; 0; 0];
nextX = zeros(N, 1); 

nextX=SearchX(H, g, N, nextX, prevX);

count = 0;

while norm(nextX - prevX) > eps
    prevX=nextX;
    
    nextX=SearchX(H, g, N, nextX, prevX);
    
    count = count + 1;
end;

disp('Погрешность:');
disp(b - A * nextX);
disp('Количество шагов:');
disp(count);

end

