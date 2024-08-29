
function [P0, S0] = get_P0S0(X, Y, bf_grid, key)

grid_size = size(X, 1);
numV = grid_size*grid_size;

if nargin < 4
    key = "normal";
end

P0_rows = zeros(numV, 1);
P0_cols = zeros(numV, 1);
P0_vals = ones(numV, 1);

S0_rows = zeros(numV, 1);
S0_cols = zeros(numV, 1);
S0_vals = zeros(numV, 1);

ind = 1;
for indV = 1:numV
    [i, j] = indV_to_ij(indV, grid_size);

    if key == "normal"
        f0 = boundary_function(X(j, i), Y(j, i), bf_grid);

        if f0 <= 0
            S0_rows(ind) = indV;
            S0_cols(ind) = indV;
            S0_vals(ind) = 1;

            P0_rows(ind) = ind;
            P0_cols(ind) = indV;
            ind = ind + 1;
        end

    elseif key == "tangential"
        f0 = boundary_function(X(j,i) - 1/2, Y(j,i) - 1/2, bf_grid);
        f1 = boundary_function(X(j,i) + 1/2, Y(j,i) - 1/2, bf_grid);
        f2 = boundary_function(X(j,i) - 1/2, Y(j,i) + 1/2, bf_grid);
        f3 = boundary_function(X(j,i) + 1/2, Y(j,i) + 1/2, bf_grid);
        
        if (f0 <= 0 || f1 <= 0 || f2 <= 0 || f3 <= 0)
            S0_rows(ind) = indV;
            S0_cols(ind) = indV;
            S0_vals(ind) = get_dual_Farea(f0, f1, f2, f3);

            P0_rows(ind) = ind;
            P0_cols(ind) = indV;
            ind = ind + 1;
        end
    end
end

S0 = sparse(S0_rows(1:ind-1), S0_cols(1:ind-1), S0_vals(1:ind-1), numV, numV);
P0 = sparse(P0_rows(1:ind-1), P0_cols(1:ind-1), P0_vals(1:ind-1), ind-1, numV);

end

