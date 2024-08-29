
function [P3, S3] = get_P3S3(X, Y, Z, bf_grid, key)

grid_size = size(X, 1);

if nargin < 5
    key = "normal";
end

epsilon = 1e-5;
numC = (grid_size-1)^3;

P3_rows = zeros(numC, 1);
P3_cols = zeros(numC, 1);
P3_vals = ones(numC, 1);

S3_rows = zeros(numC, 1);
S3_cols = zeros(numC, 1);
S3_vals = zeros(numC, 1);

ind = 1;
for indC = 1: numC

    [i, j, k] = indC_to_ijk(indC, grid_size);

    if key == "normal"

        % get all vertices of this cell

        VC = [[X(j,i,k), Y(j,i,k), Z(j,i,k)];
              [X(j,i+1,k), Y(j,i+1,k), Z(j,i+1,k)];
              [X(j+1,i+1,k), Y(j+1,i+1,k), Z(j+1,i+1,k)];
              [X(j+1,i,k), Y(j+1,i,k), Z(j+1,i,k)];
              [X(j,i,k+1), Y(j,i,k+1), Z(j,i,k+1)];
              [X(j,i+1,k+1), Y(j,i+1,k+1), Z(j,i+1,k+1)];
              [X(j+1,i+1,k+1), Y(j+1,i+1,k+1), Z(j+1,i+1,k+1)];
              [X(j+1,i,k+1), Y(j+1,i,k+1), Z(j+1,i,k+1)]]; 

        bf_VC = zeros(8,1);
        for i=1:8
            bf_VC(i) = boundary_function(VC(i,:), bf_grid);
        end

        if any(sign(bf_VC) == -1) % normal 
            % the dual cell volume is unchanged
            S3_rows(ind) = indC;
            S3_cols(ind) = indC;
            S3_vals(ind) = 1/max(get_Cvol(VC, bf_VC, bf_grid),epsilon);

            P3_rows(ind) = ind;
            P3_cols(ind) = indC;
            ind = ind + 1;

        end
        
    elseif key == "tangential"

        v0 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];

        f0 = boundary_function(v0, bf_grid);

        if f0<=0 
            % the primal cell volume is unchanged
            S3_rows(ind) = indC;
            S3_cols(ind) = indC;
            S3_vals(ind) = 1;

            P3_rows(ind) = ind;
            P3_cols(ind) = indC;
            ind = ind + 1;

        end

    end

end

S3 = sparse(S3_rows(1:ind-1), S3_cols(1:ind-1), S3_vals(1:ind-1), numC, numC);
P3 = sparse(P3_rows(1:ind-1), P3_cols(1:ind-1), P3_vals(1:ind-1), ind-1, numC);

end

