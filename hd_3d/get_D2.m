
function D2 = get_D2(grid_size)

numF = 3*(grid_size-1)*(grid_size-1)*grid_size;
numC = (grid_size-1)^3;

% D2 = sparse(numC, numF);

rows = zeros(6*numC, 1);
cols = zeros(6*numC, 1);
vals = zeros(6*numC, 1);

ind = 1;
for indC = 1:numC
    
    [i, j, k] = indC_to_ijk(indC, grid_size);
    
    rows(ind:ind+5) = indC;
    
    cols(ind) = ijkd_to_indF(i, j, k, 0, grid_size);
    vals(ind) = -1;

    cols(ind+1) = ijkd_to_indF(i, j, k, 1, grid_size);
    vals(ind+1) = -1;

    cols(ind+2) = ijkd_to_indF(i, j, k, 2, grid_size);
    vals(ind+2) = -1;
    
    cols(ind+3) = ijkd_to_indF(i+1, j, k, 2, grid_size);
    vals(ind+3) = 1;

    cols(ind+4) = ijkd_to_indF(i, j+1, k, 1, grid_size);
    vals(ind+4) = 1;

    cols(ind+5) = ijkd_to_indF(i, j, k+1, 0, grid_size);
    vals(ind+5) = 1;
    ind = ind + 6;

end

D2 = sparse(rows, cols, vals, numC, numF);

end