
% is it necessary to vectorized this one?
function indE = ijd_to_indE(i, j, direction, grid_size)

if direction == 0
    indE = i + (grid_size-1)*(j-1);
else
    indE = i + grid_size*(j-1) + (grid_size-1)*grid_size;
end

end