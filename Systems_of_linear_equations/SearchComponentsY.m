function [nextY, prevY] = SearchComponentsY(nextY, A)
% finding maximum modulo element and saving its index
indexY = [abs(nextY(1)); abs(nextY(2)); abs(nextY(3))] == max([abs(nextY(1)); abs(nextY(2)); abs(nextY(3))]);

prevY = nextY / nextY(indexY);
nextY = A * prevY;

end

