function U=MechanicalQuadratures(H, f, a, b, N)
%Метод механических квадратур
syms x;
syms y;

h=(b-a)/N;
 
fk=zeros(N, 1);
A=zeros(N);

for i=1:1:N
    xk=(a+(i-1)*h+a+i*h)/2;
    fk(i)=-subs(f, xk);
    for j=1:1:N
       xj=(a+(j-1)*h+a+j*h)/2;
       A(i, j)=subs(H, {x, y}, {xk, xj})*(a+j*h-(a+(j-1)*h));
    end;
A(i, i)=A(i, i)-1;
end;

Uk=A\fk;
sum=0;

for i=1:1:N
    sum=sum+(a+i*h-(a+(i-1)*h))*subs(H, y, a+i*h)*Uk(i);
end;

U=f+sum;

end