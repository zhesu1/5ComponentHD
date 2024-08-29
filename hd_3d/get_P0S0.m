
function [P0, S0] = get_P0S0(X, Y, Z, bf_grid, key)

grid_size = size(X, 1);
numV = grid_size^3;

if nargin < 5
    key = "normal";
end

P0_rows = zeros(numV, 1);
P0_cols = zeros(numV, 1);
P0_vals = ones(numV, 1);

S0_rows = zeros(numV, 1);
S0_cols = zeros(numV, 1);
S0_vals = zeros(numV, 1);

ind = 1;
for indV = 1: numV
    
    % disp(indV)

    [i, j, k] = indV_to_ijk(indV, grid_size);
    
    if key == "normal"
        
        v0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];

        f0 = boundary_function(v0, bf_grid);

        if f0<=0 % normal
            % the dual cell volume is unchanged
            S0_rows(ind) = indV;
            S0_cols(ind) = indV;
            S0_vals(ind) = 1;

            P0_rows(ind) = ind;
            P0_cols(ind) = indV;
            ind = ind + 1;
        end

    elseif key == "tangential"
    
        VC = [[X(j,i,k)-1/2, Y(j,i,k)-1/2, Z(j,i,k)-1/2];
              [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)-1/2];
              [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
              [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
              [X(j,i,k)-1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
              [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
              [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
              [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2]];

        bf_VC = zeros(8,1);
        for i=1:8
            bf_VC(i) = boundary_function(VC(i,:), bf_grid);
        end

        if any(sign(bf_VC) == -1)

            S0_rows(ind) = indV;
            S0_cols(ind) = indV;
            S0_vals(ind) = get_Cvol(VC, bf_VC, bf_grid);

            P0_rows(ind) = ind;
            P0_cols(ind) = indV;
            ind = ind + 1;

        end

    end

end

S0 = sparse(S0_rows(1:ind-1), S0_cols(1:ind-1), S0_vals(1:ind-1), numV, numV);
P0 = sparse(P0_rows(1:ind-1), P0_cols(1:ind-1), P0_vals(1:ind-1), ind-1, numV);

end
