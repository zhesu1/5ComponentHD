

% this function computes the projection to edges that interect the boundary
% and in addition the weights with repect to the end points for each edge.

function [P1E, P0E] = get_support_E(X, Y, bf_grid)

grid_size = size(X, 1);

% epsilon = 1e-5;
numE = 2*(grid_size-1)*grid_size;
numV = grid_size*grid_size;

ind_Vs = [];
ind_Es = [];
for indE = 1: numE

    [i, j, d] = indE_to_ijd(indE, grid_size);

    indV = ij_to_indV(i, j, grid_size);

    if d == 0 
        f0 = boundary_function(X(j, i), Y(j,i), bf_grid);
        f1 = boundary_function(X(j,i+1), Y(j,i+1), bf_grid);

        indV1 = ij_to_indV(i+1, j, grid_size);
    else
        f0 = boundary_function(X(j,i), Y(j, i), bf_grid);
        f1 = boundary_function(X(j+1, i), Y(j+1, i), bf_grid);

        indV1 = ij_to_indV(i, j+1, grid_size);
    end

    if f0<=0 || f1<=0
        ind_Es = [ind_Es, indE];
        ind_Vs = [ind_Vs, indV, indV1];
    end

end

ind_Vs = sort(unique(ind_Vs));
ind_Es = sort(unique(ind_Es));

P0E_rows = 1:length(ind_Vs);
P0E_cols = ind_Vs;
P0E_vals = ones(length(ind_Vs), 1);

P1E_rows = 1:length(ind_Es);
P1E_cols = ind_Es;
P1E_vals = ones(length(ind_Es), 1);

P0E = sparse(P0E_rows, P0E_cols, P0E_vals, length(ind_Vs), numV);
P1E = sparse(P1E_rows, P1E_cols, P1E_vals, length(ind_Es), numE);

end
