%������� ������� �������� ���������
syms x;
N=5;
h= zeros(N,N);

for i = 1 : 1 : N
    for j = 1 : 1 : N 
      h(i,j)= 1/(i+j-1);
    end;
end;
disp('h=');
disp(h);

r = rand(N);
disp('r=');
disp(r);

A = r*h;
disp('A=');
disp(A);

b = rand(N,1);
disp('b=');
disp(b);


%����� �������
disp('����� �������:');
disp('x=');
xi=CramerMethod(A, b, N);
disp(vpa(xi));

%������ ��������� ��������� � ������� � �� ������ �������.
disp('����� ������� � �����������:');
h = 0.1;
b1 = zeros(N, 1);
for i = 1 : 1 : N
    b1(i)= b(i)+ h;
end;

disp('����� ������ b:');
disp(b1);

disp('x=');
xk=CramerMethod(A, b1, N);
disp(vpa(xk));

disp('�����������');
for i = 1 : 1 : N
    disp(vpa(xi(i)- xk(i)));
end;

solution1 = sumabs(xi);
difference = sumabs(xi-xk);
error = vpa(difference/solution1);
disp(error);

%������ �����������

%����� ������
disp('||�||=');
normA = norm(A,1);
disp(normA);

B = inv(A);

disp('||B||=');
normB=norm(B);
disp(normB);

disp('����� ���������������=');
cond1 = vpa(normA * normB);
disp(cond1);

normb=sumabs(b);
normb1=sumabs(b1);

disp('������������� �����������:');
disp(cond1*(normb1/normb)- error );
