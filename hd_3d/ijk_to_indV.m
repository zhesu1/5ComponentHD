

function indV = ijk_to_indV(i, j, k, grid_size)

indV = i + grid_size * (j-1) + grid_size * grid_size * (k-1);

end