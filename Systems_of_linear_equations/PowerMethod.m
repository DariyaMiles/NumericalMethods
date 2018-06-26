function [nextLambda] = PowerMethod(A, eps)
disp('    Power method:');

prevY =[ 1; 0; 0];
nextY = A * prevY;

index = 1;

while nextY(index) / prevY(index) == 0
    index = index + 1;
end;

prevLambda = nextY(index) / prevY(index);

[nextY, prevY] = SearchComponentsY(nextY, A);

index = 1;

while nextY(index) / prevY(index) == 0 
    index = index + 1;
end;

nextLambda = nextY(index) / prevY(index);

count = 1;

while abs(nextLambda - prevLambda) > eps
    count = count + 1;
    [nextY, prevY] = SearchComponentsY(nextY, A);
    
    prevLambda = nextLambda;
    
    index = 1;
    
    while nextY(index) / prevY(index) == 0 
        index = index + 1;
    end;
    
    nextLambda = nextY(index) / prevY(index);
end;


disp('      ÊNumber of steps:');
disp(count);

end

