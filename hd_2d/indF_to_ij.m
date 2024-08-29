
function [i, j] = indF_to_ij(indF, grid_size)

i = mod(indF-1, grid_size-1) + 1;
j = floor((indF-1)/(grid_size-1)) + 1;

end
