function [ nextX ] = SearchX( H, g, N, nextX, prevX )

for i = 1:1:N
    sum1 = 0;
    sum2 = 0;
    
    for j = 1:1:i-1
        sum1 = sum1 + H(i, j) * nextX(j);
    end

    for j = i:1:N
        sum2 = sum2 + H(i, j)*prevX(j);
    end

    nextX(i) = sum1 + sum2 + g(i);
end

end

