% Разностный метод решения краевой задачи для ОДУ 2-го порядка

syms x;
N=10;
a=0;
b=1;
m=5;
alpha=0.1*m;
h=(b-a)/N;

P=x^2+alpha;
Q=2*x;
F=2*(3*x^2-alpha)/((x^2+alpha)^3);

alpha0=1;
alpha1=-2;
A=1/alpha;

beta0=1;
beta1=0;
B=1/(1+alpha);

solution=1/(x^2+alpha);

solution_xi = zeros(1,N+1);

f=zeros(1,N+1);
p=zeros(1,N+1);
q=zeros(1,N+1);

% Основная сетка
xi=zeros(1,N+1);
for i=1:1:(N+1)
    xi(i)=a+(i-1)*h;
    
    solution_xi(i)=subs(solution,x,xi(i));
    
    f(i)=subs(F,x,xi(i));
    p(i)=subs(P,x,xi(i));
    q(i)=subs(Q,x,xi(i));
end;

p1=alpha1/(alpha1-alpha0*h);
d1=-A*h/(alpha1-alpha0*h);
p2=beta1/(beta0*h+beta1);
d2=B*h/(beta0*h+beta1);

disp('Метод прогонки для основной сетки');
[m, k]=SweepMethodCoef( N, h, p, q, f, p1, d1 );

% Обратный ход метода прогонки
y=zeros(1,N+1);
y(N+1)=(p2*k(N)+d2)/(1-p2*m(N));

for i=N:-1:1
    y(i)=m(i)*y(i+1)+k(i);
end;

disp('  Фактическая погрешность в узлах сетки:');
for i=1:1:(N+1)
    disp(vpa(abs(y(i)-solution_xi(i))));
end;

% Сдвинутая сетка
for i=1:1:(N+2)
    xi(i)=a-h/2+(i-1)*h;
    
    solution_xi(i)=subs(solution,x,xi(i));
    
    f(i)=subs(F,x,xi(i));
    p(i)=subs(P,x,xi(i));
    q(i)=subs(Q,x,xi(i));
end;

p1=(alpha0*h+2*alpha1)/(2*alpha1-alpha0*h);
d1=2*A*h/(alpha0*h-2*alpha1);
p2=(2*beta1-beta0*h)/(beta0*h+2*beta1);
d2=2*B*h/(beta0*h+2*beta1);

disp('Метод прогонки для сдвинутой сетки');
% добавляется два узла
[m, k]=SweepMethodCoef( N+2, h, p, q, f, p1, d1 );

y=zeros(1,N+2);
y(N+2)=(p2*k(N)+d2)/(1-p2*m(N));

for i=(N+1):-1:1
    y(i)=m(i)*y(i+1)+k(i);
end;

disp('  Фактическая погрешность в узлах сетки:');
for i=1:1:(N+2)
    disp(vpa(abs(y(i)-solution_xi(i))));
end;
