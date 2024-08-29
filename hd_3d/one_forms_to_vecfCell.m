
function vecfCell_grid = one_forms_to_vecfCell(omega, grid_size)

numC = (grid_size-1)^3;

vecfCell_grid = zeros(grid_size-1, grid_size-1, grid_size-1, 3);

for indC=1:numC

    [i, j, k] = indC_to_ijk(indC, grid_size);


    indE0 = ijkd_to_indE(i, j, k, 0, grid_size); % x direction
    indE1 = ijkd_to_indE(i, j+1, k, 0, grid_size);
    indE2 = ijkd_to_indE(i, j, k+1, 0, grid_size); 
    indE3 = ijkd_to_indE(i, j+1, k+1, 0, grid_size);

    indE4 = ijkd_to_indE(i, j, k, 1, grid_size); % y direction
    indE5 = ijkd_to_indE(i+1, j, k, 1, grid_size);
    indE6 = ijkd_to_indE(i, j, k+1, 1, grid_size); 
    indE7 = ijkd_to_indE(i+1, j, k+1, 1, grid_size);

    indE8 = ijkd_to_indE(i, j, k, 2, grid_size); % z direction
    indE9 = ijkd_to_indE(i, j+1, k, 2, grid_size);
    indE10 = ijkd_to_indE(i+1, j, k, 2, grid_size); 
    indE11 = ijkd_to_indE(i+1, j+1, k, 2, grid_size);
    

    vecfCell_grid(j, i, k, 1) = (omega(indE0) + omega(indE1) ...
        + omega(indE2) + omega(indE3))/4;
    vecfCell_grid(j, i, k, 2) = (omega(indE4) + omega(indE5) ...
        + omega(indE6) + omega(indE7))/4;
    vecfCell_grid(j, i, k, 3) = (omega(indE8) + omega(indE9) ...
        + omega(indE10) + omega(indE11))/4;

end
end