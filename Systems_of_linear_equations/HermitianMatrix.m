function [ nextX ] = HermitianMatrix( A, b, E, eps )
disp('����� ������� ��������');
disp('2. �������� ������������ ������������ �������:');

% ������������ �� ������ ����������� ����� � ��������� �������
maxEigenvalue=PowerMethod(A, eps);
disp('      ������������ �� ������ ����������� ����� �=');
disp(maxEigenvalue);

% ������������ �� ������ ����������� ����� � ������� ���������� ������������
maxEigenvalue = ScalarProduct(A, eps);
disp('      ������������ �� ������ ����������� ����� �=');
disp(maxEigenvalue);

% ����������� �� ������ ����������� ����� � ������� ���������� ������������
B = A - maxEigenvalue * E;
minEigenvalue=maxEigenvalue +ScalarProduct(B, eps);
disp('      ����������� �� ������ ����������� ����� �=');
disp(minEigenvalue);

alfa = 2/(minEigenvalue+maxEigenvalue);
H = E - alfa*A;
g = alfa * b;

prevX = [1; 0; 0];
nextX = H * prevX + g;

count = 0;

while norm(nextX - prevX) > eps
    prevX = nextX;
    nextX = H * prevX + g;
    count = count + 1;
end; 

disp('�����������:');
disp(b - A * nextX);
disp('���������� �����:');
disp(count);

end

