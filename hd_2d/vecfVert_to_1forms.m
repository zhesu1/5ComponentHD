
function omega = vecfVert_to_1forms(U)

    % U: grid_size x grid_size x 2
    grid_size = size(U, 1);
    numV = grid_size*grid_size;

    for indV = 1:numV
        [i, j] = indV_to_ij(indV, grid_size);
        u(indV,:) = U(j, i, :);
    end

    % u: #V x 2
    % omega: #E
    omega = vflds_to_1forms(u, grid_size);

end