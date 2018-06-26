function [nextLambda] =ScalarProduct( A, eps )
disp('    Метод скалярного произведения:');

prevY = [1; 0; 0];
nextY = A * prevY;

prevLambda = sum(nextY.*prevY)/sum(prevY.*prevY);

[nextY, prevY]=SearchComponentsY(nextY, A);

nextLambda = sum(nextY.*prevY) / sum(prevY.*prevY);

count = 1;

while abs(nextLambda-prevLambda) > eps
    count = count + 1;
    
    [nextY, prevY]=SearchComponentsY(nextY, A);
    
    prevLambda = nextLambda;
    nextLambda = sum(nextY.*prevY) / sum(prevY.*prevY);
end

disp('      Количество шагов:');
disp(count);

end

