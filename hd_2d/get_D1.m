
% v2
function D1 = get_D1(grid_size)

numE = 2*(grid_size-1)*grid_size;
numF = (grid_size-1)*(grid_size-1);

rows = zeros(4*numF, 1);
cols = zeros(4*numF, 1);
vals = zeros(4*numF, 1);

ind = 1;
for indF = 1:numF
    [i, j] = indF_to_ij(indF, grid_size);

    rows(ind:ind+3) = indF;
    cols(ind) = ijd_to_indE(i, j, 0, grid_size);
    vals(ind) = 1;

    cols(ind+1) = ijd_to_indE(i, j, 1, grid_size);
    vals(ind+1) = -1;

    cols(ind+2) = ijd_to_indE(i+1, j, 1, grid_size);
    vals(ind+2) = 1;

    cols(ind+3) = ijd_to_indE(i, j+1, 0, grid_size);
    vals(ind+3) = -1;
    ind = ind + 4;
end

D1 = sparse(rows, cols, vals, numF, numE);

end

