function [x] = ModificationOfGaussianMethod3(A, b, n)
Ab = [A b];
x = [0; 0; 0];
tmp = [1 2 3];

for k = 1 : 1 : n
     max = Ab(k, k);
     num_i = k;
     num_j = k;
     for i = k : 1 : n
         for j = k : 1 : n
           if abs(Ab(i, j)) > abs(max)
                 max = Ab(i, j);
                 num_i = i;
                 num_j = j;
            end;
         end;
     end;
     
     A1 = Ab(k, :);
     Ab(k, :) = Ab(num_i, :);
     Ab(num_i, :) = A1;
     
     A2 = Ab(:, k);
     Ab(:, k) = Ab(:, num_j);
     Ab(:, num_j) = A2;
     
     tmp1 = tmp(k);
     tmp(k) = tmp(num_j);
     tmp(num_j) = tmp1;
     
     % Gaussian method
     Abkk = Ab(k, k);
     for j = k : 1 : (n + 1)
         if Abkk ~= 0
             Ab(k, j) = Ab(k, j) / Abkk;
         end;
     end;
    
     for i = (k + 1) : 1 : n
         Abik = Ab(i, k);
         for j = 1 : 1 : (n + 1)
            Ab(i, j) = Ab(i, j) - Ab(k, j) * Abik;
         end;
     end;
end;

for i = n : -1 : 1
    sum = 0;
    for j = (i + 1) : 1 : n
        sum = sum + Ab(i, j) * x(tmp(j));
    end;
    x(tmp(i)) = Ab(i, n + 1) - sum;
end;

disp('x = ');
disp(x);
disp('Error');
disp('b - Ax = ');
disp(b - A * x);

end
