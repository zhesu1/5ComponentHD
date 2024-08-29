

function indE = ijkd_to_indE(i, j, k, d, grid_size)

if d == 0 % x positive
    indE = i + (grid_size - 1) * (j-1) + (grid_size - 1) * grid_size * (k-1);
elseif d == 1 % y positive
    indE = i + grid_size * (j-1) + (grid_size - 1) * grid_size * (k-1) + (grid_size-1) * grid_size * grid_size ;
else % z positive
    indE = i + grid_size * (j-1) + grid_size * grid_size * (k-1) + 2 * (grid_size - 1) * grid_size * grid_size;
end
end
