
function [i, j] = indV_to_ij(indV, grid_size)

i = mod(indV - 1, grid_size) + 1;
j = floor((indV-1)/grid_size) + 1;

end