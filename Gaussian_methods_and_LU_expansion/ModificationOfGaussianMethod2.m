function [x]=ModificationOfGaussianMethod2(A,b,n)
Ab=[A b];
x=[0; 0; 0];
tmp=[1 2 3];

for k=1:1:n
     max=Ab(k,k);
     num=k;
     for j=(k+1):1:n
         if abs(Ab(k,j))>abs(max)
             max=Ab(k,j);
             num=j;
         end;
     end;
     
     A1=Ab(:,k);
     Ab(:,k)=Ab(:,num);
     Ab(:,num)=A1;
     
     tmp1=tmp(k);
     tmp(k)=tmp(num);
     tmp(num)=tmp1;
     
     %Gaussian method
     Abkk=Ab(k,k);
     for j=k:1:(n+1)
         if Abkk~=0
             Ab(k,j)=Ab(k,j)/Abkk;
         end;
     end;
    
     for i=(k+1):1:n
         Abik=Ab(i,k);
         for j=1:1:(n+1)
            Ab(i,j)=Ab(i,j)-Ab(k,j)*Abik;
         end;
     end;
end;

for i=n:-1:1
    sum=0;
    for j=(i+1):1:n
        sum=sum+Ab(i,j)*x(tmp(j));
    end;
    x(tmp(i))=Ab(i,n+1)-sum;
end;

disp('x=');
disp(x);
disp('Error');
disp('b-Ax=');
disp(b-A*x);

end