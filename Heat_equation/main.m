% Решение уравнения теплопроводности с постоянными коэффициентами

syms x;
syms t;

ai = 0;
t0 = 0;
T = 0.02;

f = (-0.25)*exp(-0.25*t)*(1-x);
solution_u = exp(-0.25*t)*(sin(0.5*x)+1-x);
U0 = sin(0.5*x)+1-x; % U0=U(0,x)
U1 = exp(-0.25*t); % U1=U(0,t)
U2 = exp(-0.25*t)*sin(0.5); % U2=U(1,t)

%Явная схема
disp('Явная схема:');

h = 0.1;
tau = 0.004;

%Введем сетку
N = 1/h;
M = T/tau;

xi = zeros(1,N+1);
tj = zeros(1,M+1);
f_ji = zeros(M+1,N+1);
U_ji = zeros(M+1,N+1);
solution_U = zeros(M+1,N+1); 

% Подстановка краевых условий и введение обозначений f_ij, xi, tj
for i=1:1:N+1
    xi(i) = ai + (i-1)*h;
    
    U_ji(1,i) = subs(U0, x, xi(i));
    for j=1:1:M+1
        tj(j) = t0 + (j-1)*tau;
        f_ji(j,i) = subs(subs(f, x, xi(i)), t, tj(j));
        
        % значения приближенного решения на нижнем и верхнем слое
        U_ji(j,1) = subs(U1, t, tj(j));
        U_ji(j,N+1) = subs(U2, tj(j));
        
        solution_U(j,i) = subs(subs(solution_u, x, xi(i)), t, tj(j));
    end;
end;

% значения приближенного решения во внутренних точках
for j=1:1:M
    for i=2:1:N
        U_ji(j+1,i) = (1-2*tau/h^2)*U_ji(j,i) + tau/h^2*(U_ji(j, i-1) + U_ji(j, i+1)) + tau*f_ji(j,i);
    end;
end;

disp('   Фактическая погрешность =');
disp(solution_U - U_ji);

%Неявная схема
disp('Неявная схема:');

h = 0.1;
tau = 0.02;

%Введем сетку
N = 1/h;
M = T/tau;

xi = zeros(1,N+1);
tj = zeros(1,M+1);
f_ji = zeros(M+1,N+1);
U_ji = zeros(M+1,N+1);
solution_U = zeros(M+1,N+1);

% Подстановка краевых условий и введение обозначений f_ij, xi, tj
for i=1:1:N+1
    xi(i) = (i-1)*h;
    
    U_ji(1,i) = subs(U0, x, xi(i)); 
    
    for j=1:1:M+1
        tj(j) = t0 + (j-1)*tau;
        f_ji(j,i) = subs(subs(f, x, xi(i)), t, tj(j));
        
        U_ji(j,1) = subs(U1, tj(j));
        U_ji(j,N+1) = subs(U2, tj(j));
        
        solution_U(j,i) = subs(subs(solution_u, x, xi(i)), t, tj(j));
    end;
end;

%Метод прогонки для трехдиагональной матрицы для слоя 1

% коэффициенты прогонки:
ai =  tau/(2*h^2);
bi = tau/(2*h^2);
ci = 1 + tau/h^2;

F_ji = zeros(1, N+1);
j = 1;

for i=2:1:N
    F_ji(j, i) = (1-tau/h^2)*U_ji(j,i) + tau/(2*h^2)*(U_ji(j,i-1) + U_ji(j,i+1))+tau*f_ji(j,i); 
end;

%Из граничных условий:
p1 = 0; 
d1 = U_ji(j+1, 1);

p2 = 0;
d2 = U_ji(j+1, N+1);

%Первый этап прогонки:
m = zeros(1, N);
k = zeros(1, N);

m(1) = p1;
k(1) = d1;

for i=2:1:N
    m(i) = bi/(ci - ai*m(i-1));
    k(i) = (ai*k(i-1) + F_ji(i))/(ci - ai*m(i-1));
end;

%Второй этап
U_ji(j+1,N+1) = (p2*k(N) + d2)/(1 - p2*m(N));

for i=N:-1:2
    U_ji(j+1, i) = m(i)* U_ji(j+1, i+1) + k(i); 
end;

disp('   Фактическая погрешность =');
disp(solution_U - U_ji);
