
function D0 = get_D0(grid_size)

numV = grid_size^3;
numE = 3*(grid_size-1)*grid_size*grid_size;

rows = zeros(2*numE, 1);
cols = zeros(2*numE, 1);
vals = zeros(2*numE, 1);

ind = 1;
for indE = 1:numE
    
    [i,j,k,d] = indE_to_ijkd(indE, grid_size);
    
    rows(ind) = indE;
    cols(ind) = ijk_to_indV(i,j,k,grid_size);
    vals(ind) = -1;
    ind = ind + 1;

    if d==0
        indV = ijk_to_indV(i+1, j, k, grid_size);
    elseif d==1
        indV = ijk_to_indV(i, j+1, k, grid_size);
    else
        indV = ijk_to_indV(i, j, k+1, grid_size);
    end
    rows(ind) = indE;
    cols(ind) = indV;
    vals(ind) = 1;
    ind = ind + 1;

end

D0 = sparse(rows, cols, vals, numE, numV);

end

