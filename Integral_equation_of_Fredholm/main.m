%Численное решение интегрального уравнения Фредгольма 2 рода
syms x;
syms y;
a=0;
b=1;
H=cos(x*y)/2;
f=1+x-x^2;

%Метод замены ядра на вырожденное

%Разложение в ряд Тейлора в сумму из 3 слагаемых:
Ht=taylor(H, x, a, 'Order', 5);
disp('Разложение в ряд Тейлора в сумму из 3 слагаемых:');
disp(Ht);
alpha= [1/2, -x^2/4, x^4/48];
beta= [1, y^2, y^4];
U3=ReplacingKernel(alpha, beta, f, a, b, 3);
disp('U(a)=');
disp(vpa(subs(U3, a)));
disp('U((a+b)/2)=');
disp(vpa(subs(U3, (a+b)/2)));
disp('U(b)=');
disp(vpa(subs(U3, b)));

%Error
disp('Фактическая погрешность в точке "a"');
disp(vpa(subs(U3-int(H*subs(U3, y), y, a, b)-f, 0.01)));
disp('Фактическая погрешность в точке "(a+b)/2"');
disp(vpa(subs(U3-int(H*subs(U3, y), y, a, b)-f, (a+b)/2)));
disp('Фактическая погрешность в точке "b"');
disp(vpa(subs(U3-int(H*subs(U3, y), y, a, b)-f, b)));

%Разложение в ряд Тейлора в сумму из 4 слагаемых:
Ht=taylor(H, x, a, 'Order', 7);
disp('Разложение в ряд Тейлора в сумму из 4 слагаемых:');
disp(Ht);
alpha= [1/2, -x^2/4, x^4/48, -x^6/1440];
beta= [1, y^2, y^4, y^6];
U4=ReplacingKernel(alpha, beta, f, a, b, 4);
disp('U(a)=');
disp(vpa(subs(U4, a)));
disp('U((a+b)/2)=');
disp(vpa(subs(U4, (a+b)/2)));
disp('U(b)=');
disp(vpa(subs(U4, b)));

%Error
disp('Фактическая погрешность в точке "a"');
disp(vpa(subs(U4-int(H*subs(U4, y), y, a, b)-f, 0.01)));
disp('Фактическая погрешность в точке "(a+b)/2"');
disp(vpa(subs(U4-int(H*subs(U4, y), y, a, b)-f, (a+b)/2)));
disp('Фактическая погрешность в точке "b"');
disp(vpa(subs(U4-int(H*subs(U4, y), y, a, b)-f, b)));

%Метод механических квадратур
N=2;
disp('Метод механических квадратур:');
Um=MechanicalQuadratures(H, f, a, b, N);
disp('N=');
disp(N);
disp('U(a)=');
disp(vpa(subs(Um, a)));
disp('U((a+b)/2)=');
disp(vpa(subs(Um, (a+b)/2)));
disp('U(b)=');
disp(vpa(subs(Um, b)));

N=2*N;
Un=MechanicalQuadratures(H, f, a ,b, N);
xk=[a; (a+b)/2; b];
eps=0.0001;

%Удваиваем узлы пока разность не будет меньше eps
while(max(abs(subs(Un, xk)-subs(Um, xk)), xk))>eps
    N=2*N;
    Um=Un;
    Un=MechanicalQuadratures(H, f, a ,b, N);
end;

disp('Количество узлов N=');
disp(N);
disp('U(a)=');
disp(vpa(subs(Un, a)));
disp('U((a+b)/2)=');
disp(vpa(subs(Un, (a+b)/2)));
disp('U(b)=');
disp(vpa(subs(Un, b)));

disp('Сравним результаты с методом замены ядра(3 слагаемых):');
disp('в точке "a"');
disp(abs(vpa(subs(Un, a))-vpa(subs(U3, a))));
disp('в точке "(a+b)/2"');
disp(abs(vpa(subs(Un, (a+b)/2))-vpa(subs(U3, (a+b)/2))));
disp('в точке "b"');
disp(abs(vpa(subs(Un, b))-vpa(subs(U3, b))));

disp('Сравним результаты с методом замены ядра(4 слагаемых):');
disp('в точке "a"');
disp(abs(vpa(subs(Un, a))-vpa(subs(U4, a))));
disp('в точке "(a+b)/2"');
disp(abs(vpa(subs(Un, (a+b)/2))-vpa(subs(U4, (a+b)/2))));
disp('в точке "b"');
disp(abs(vpa(subs(Un, b))-vpa(subs(U4, b))));
