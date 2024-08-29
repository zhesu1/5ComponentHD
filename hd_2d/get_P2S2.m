

function [P2, S2] = get_P2S2(X, Y, bf_grid, key)

grid_size = size(X, 1);

if nargin < 4
    key = "normal";
end

epsilon = 1e-5;
numF = (grid_size-1)*(grid_size-1);

P2_rows = zeros(numF, 1);
P2_cols = zeros(numF, 1);
P2_vals = ones(numF, 1);

S2_rows = zeros(numF, 1);
S2_cols = zeros(numF, 1);
S2_vals = zeros(numF, 1);

ind = 1;
for indF = 1: numF

    [i, j] = indF_to_ij(indF, grid_size);

    if key == "normal"

        f0 = boundary_function(X(j,i), Y(j,i), bf_grid);
        f1 = boundary_function(X(j,i+1), Y(j,i+1), bf_grid);
        f2 = boundary_function(X(j+1,i), Y(j+1,i), bf_grid);
        f3 = boundary_function(X(j+1, i+1), Y(j+1, i+1), bf_grid);

        if(f0<=0 || f1<=0 || f2<=0 || f3<=0) % normal % the dual cell volume is unchanged

            S2_rows(ind) = indF;
            S2_cols(ind) = indF;
            S2_vals(ind) = 1/max(get_primal_Farea(f0, f1, f2, f3), epsilon^2);

            P2_rows(ind) = ind;
            P2_cols(ind) = indF;
            ind = ind + 1;
        end
        
    elseif key == "tangential"

        f0 = boundary_function(X(j,i)+1/2, Y(j,i)+1/2, bf_grid);

        if f0 <= 0 % tangential % the primal cell volume is unchanged

            S2_rows(ind) = indF;
            S2_cols(ind) = indF;
            S2_vals(ind) = 1;

            P2_rows(ind) = ind;
            P2_cols(ind) = indF;
            ind = ind + 1;
        end
    end

end

S2 = sparse(S2_rows(1:ind-1), S2_cols(1:ind-1), S2_vals(1:ind-1), numF, numF);
P2 = sparse(P2_rows(1:ind-1), P2_cols(1:ind-1), P2_vals(1:ind-1), ind-1, numF);

end

