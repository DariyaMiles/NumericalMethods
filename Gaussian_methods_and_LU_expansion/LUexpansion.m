function [x]=LUexpansion(A,b,n)
x=[0; 0; 0];
y=[0; 0; 0];
L=[0 0 0; 0 0 0; 0 0 0];
U=[1 0 0; 0 1 0; 0 0 1];

%Находим L и U
for i=1:1:n
    for j=1:1:n
        if i>=j
            sum1=0;
            for k=1:1:(j-1)
                sum1=sum1+L(i,k)*U(k,j);
            end;
            L(i,j)=A(i,j)-sum1;
        else
            sum2=0;
            for k=1:1:(i-1)
                sum2=sum2+L(i,k)*U(k,j);
            end;
            U(i,j)=(A(i,j)-sum2)/L(i,i);
        end;
    end;
end;

disp('Error:');
disp('A-LU=');
disp(A-L*U);

%Решаем Ly=b
y(1)=b(1)/L(1,1);
for i=2:1:n
    sum=0;
    for j=1:1:(i-1)
     sum=sum+L(i,j)*y(j);
    end;
    y(i)=(b(i)-sum)/L(i,i);
end;

disp('y=');
disp(y);
disp('Error:');
disp('b-Ly=');
disp(b-L*y);

%Решаем Ux=y
Uy=[U y];
for i=n:-1:1
    sum=0;
    for j=(i+1):1:n
        sum=sum+Uy(i,j)*x(j);
    end;
    x(i)=Uy(i,n+1)-sum;
end;

disp('x=');
disp(x);
disp('Error:');
disp('y-Ux=');
disp(y-U*x);
disp('b-Ax=');
disp(b-A*x);
end