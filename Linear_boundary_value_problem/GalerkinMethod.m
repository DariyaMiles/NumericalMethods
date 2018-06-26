function [y] = GalerkinMethod(a, b, N, wi, Lwi, f)

disp('Galerkin method:');

% Orthogonality conditions
A = zeros(N);
B = zeros(N, 1);

for i = 1 : 1 : N
    for j = 1 : 1 : N
        A(i, j) = int(wi(i) * Lwi(j), a, b);
    end;
    
    B(i, 1) = int(wi(i) * f, a, b);
end;

disp('A = ');
disp(A);
disp('B = ');
disp(B);
ci = A \ B;
disp('ci = ');
disp(ci);
disp('Error = ');
disp(B - A * ci);

% The solution in the form of a linear combination:
y = 0;

for i = 1 : 1 : N
    y = y + ci(i) * wi(i);
end;

disp('y = ');
disp(y);

end

