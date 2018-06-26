function [ nextY, prevY ] = SearchComponentsY( nextY, A )
% Находим максимальный по модулю элемент и сохраняем его индекс
indexY=[abs(nextY(1)); abs(nextY(2)); abs(nextY(3))]== max([abs(nextY(1)); abs(nextY(2)); abs(nextY(3))]);

prevY = nextY /nextY(indexY);
nextY = A * prevY;

end

