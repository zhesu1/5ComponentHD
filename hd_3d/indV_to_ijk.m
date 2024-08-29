
function [i, j, k] = indV_to_ijk(indV, grid_size)

i = mod(indV - 1, grid_size) + 1;
j = floor(mod(indV - 1, grid_size * grid_size)/grid_size) + 1;
k = floor((indV - 1)/(grid_size * grid_size)) + 1;

end
