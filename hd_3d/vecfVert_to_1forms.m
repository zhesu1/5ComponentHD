
function omega = vecfVert_to_1forms(U)
% U: N x N x N x 3

grid_size = size(U, 1);
numV = grid_size^3;

for indV = 1:numV
    [i, j, k] = indV_to_ijk(indV, grid_size);
    u(indV,:) = U(j, i, k, :);
end

% u: #V x 3
% omega: # E

omega = vflds_to_1forms(u, grid_size);

end