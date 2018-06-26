function xi=CramerMethod(A, b, N)
xi= zeros(N, 1);

detA = vpa(det(A));

for i = 1 : 1 : N
    A1= A;
    for j = 1 : 1 : N 
        A1(j,i)=b(j);
    end;
 detA1=vpa(det(A1));
 xi(i)=vpa(detA1/detA);
end;

end