
function vecfFace_grid = one_forms_to_vecfFace(omega, grid_size)

    numF = (grid_size-1)*(grid_size-1);

    vecfFace_grid = zeros(grid_size-1, grid_size-1, 2);
    for indF=1:numF

        [i, j] = indF_to_ij(indF, grid_size);
        indE0 = ijd_to_indE(i, j, 0, grid_size);
        indE1 = ijd_to_indE(i, j+1, 0, grid_size);
        indE2 = ijd_to_indE(i, j, 1, grid_size);
        indE3 = ijd_to_indE(i+1, j, 1, grid_size);
        
        vecfFace_grid(j, i, 1) = (omega(indE0) + omega(indE1))/2;
        vecfFace_grid(j, i, 2) = (omega(indE2) + omega(indE3))/2;
        
    end

end