
function D1 = get_D1(grid_size)

numE = 3*(grid_size-1)*grid_size*grid_size;
numF = 3*(grid_size-1)*(grid_size-1)*grid_size;

% D1 = sparse(numF, numE);
rows = zeros(4*numF, 1);
cols = zeros(4*numF, 1);
vals = zeros(4*numF, 1);

ind = 1;
for indF = 1:numF
    
    [i,j,k,d] = indF_to_ijkd(indF, grid_size);
    
    rows(ind:ind+3) = indF;

    if d==0
        
        cols(ind) = ijkd_to_indE(i, j, k, 0, grid_size);
        vals(ind) = 1;

        cols(ind+1) = ijkd_to_indE(i, j, k, 1, grid_size);
        vals(ind+1) = -1;

        cols(ind+2) = ijkd_to_indE(i+1, j, k, 1, grid_size);
        vals(ind+2) = 1;

        cols(ind+3) = ijkd_to_indE(i, j+1, k, 0, grid_size);
        vals(ind+3) = -1;
        ind = ind + 4;

    elseif d==1

        cols(ind) = ijkd_to_indE(i, j, k, 0, grid_size);
        vals(ind) = -1;

        cols(ind+1) = ijkd_to_indE(i, j, k, 2, grid_size);
        vals(ind+1) = 1;

        cols(ind+2) = ijkd_to_indE(i, j, k+1, 0, grid_size);
        vals(ind+2) = 1;

        cols(ind+3) = ijkd_to_indE(i+1, j, k, 2, grid_size);
        vals(ind+3) = -1;
        ind = ind + 4;

    else

        cols(ind) = ijkd_to_indE(i, j, k, 1, grid_size);
        vals(ind) = 1;

        cols(ind+1) = ijkd_to_indE(i, j, k, 2, grid_size);
        vals(ind+1) = -1;

        cols(ind+2) = ijkd_to_indE(i, j+1, k, 2, grid_size);
        vals(ind+2) = 1;

        cols(ind+3) = ijkd_to_indE(i, j, k+1, 1, grid_size);
        vals(ind+3) = -1;
        ind = ind + 4;

    end
end

D1 = sparse(rows, cols, vals, numF, numE);

end
