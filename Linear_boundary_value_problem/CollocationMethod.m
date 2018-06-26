function [ y ] = CollocationMethod( xi, wi, Lwi, f, N)
disp('Метод коллокаций:');
syms x;

%Получаем систему уравнений:
A = zeros(N);
b = zeros(N, 1);

for i=1:1:N
    for j=1:1:N
        A(i, j)=subs(Lwi(j), xi(i));
    end;
    
    b(i,1)=subs(f, xi(i));
end;

ci = A\b;
disp('  ci=');
disp(ci);
disp('  Погрешность= ');
disp(b-A*ci);

%Решение в виде линейной комбинации:
y = 0;

for i=1:1:N
    y = y + ci(i)*wi(i);
end;

disp('  y=');
disp(y);

end

