function U = ReplacingKernel(alpha, beta, f, a ,b, n)
syms x;
syms y;

fi=zeros(n, 1);
A=zeros(n, n);

for i=1:1:n
fi(i)=-int(beta(i)*subs(f, y), y, a, b);
    for j=1:1:n
        A(i, j)=int(beta(i)*subs(alpha(j), y), y, a ,b);
    end;
A(i, i)=A(i, i)-1;
end;

disp('A=');
disp(A);
disp('fi=');
disp(fi);

disp('c=');
c=A\fi;
disp(c);

sum=0;

for i=1:1:n
    sum=sum+c(i)*alpha(i);
end;
U=f+sum;
end
