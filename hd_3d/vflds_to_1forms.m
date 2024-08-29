
function omega = vflds_to_1forms(u, grid_size)

numE = 3*(grid_size-1)*grid_size*grid_size;

omega = zeros(numE,1);

for indE = 1:numE

        [i, j, k, d] = indE_to_ijkd(indE, grid_size);

        indV = ijk_to_indV(i, j, k, grid_size);

        if d == 0
            indV1 = ijk_to_indV(i+1, j, k, grid_size);
            omega(indE) = (u(indV, 1) + u(indV1, 1))/2;
        elseif d == 1
            indV1 = ijk_to_indV(i, j+1, k, grid_size);
            omega(indE) = (u(indV, 2) + u(indV1, 2))/2;
        else
            indV1 = ijk_to_indV(i, j, k+1, grid_size);
            omega(indE) = (u(indV, 3) + u(indV1, 3))/2;
        end
end


end