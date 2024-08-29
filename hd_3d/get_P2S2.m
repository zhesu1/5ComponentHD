

function [P2, S2] = get_P2S2(X, Y, Z, bf_grid, key)

grid_size = size(X, 1);

if nargin < 5
    key = "normal";
end

epsilon = 1e-5;
numF = 3*(grid_size-1)*(grid_size-1)*grid_size;

P2_rows = zeros(numF, 1);
P2_cols = zeros(numF, 1);
P2_vals = ones(numF, 1);

S2_rows = zeros(numF, 1);
S2_cols = zeros(numF, 1);
S2_vals = zeros(numF, 1);

ind = 1;
for indF = 1: numF

    [i, j, k, d] = indF_to_ijkd(indF, grid_size);

    if key == "normal"

        if d == 0

            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j,i+1,k), Y(j,i+1,k), Z(j,i+1,k)];
            V2 = [X(j+1,i,k), Y(j+1,i,k), Z(j+1,i,k)];
            V3 = [X(j+1,i+1,k), Y(j+1,i+1,k), Z(j+1,i+1,k)];

        elseif d==1

            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j,i+1,k), Y(j,i+1,k), Z(j,i+1,k)];
            V2 = [X(j,i,k+1), Y(j,i,k+1), Z(j,i,k+1)];
            V3 = [X(j,i+1,k+1), Y(j,i+1,k+1), Z(j,i+1,k+1)];

        else

            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j+1,i,k), Y(j+1,i,k), Z(j+1,i,k)];
            V2 = [X(j,i,k+1), Y(j,i,k+1), Z(j,i,k+1)];
            V3 = [X(j+1,i,k+1), Y(j+1,i,k+1), Z(j+1,i,k+1)];

        end


        f0 = boundary_function(V0, bf_grid);
        f1 = boundary_function(V1, bf_grid);
        f2 = boundary_function(V2, bf_grid);
        f3 = boundary_function(V3, bf_grid);

        if (f0<=0 || f1<=0 || f2<=0 || f3<=0)
            % normal % the dual cell volume is unchanged
            S2_rows(ind) = indF;
            S2_cols(ind) = indF;
            S2_vals(ind) = 1/max(get_Farea(f0, f1, f2, f3), epsilon);

            P2_rows(ind) = ind;
            P2_cols(ind) = indF;
            ind = ind + 1;
        end

    elseif key == "tangential"

        if d == 0

            V0 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
            
        elseif d==1

            V0 = [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
           
        else

            V0 = [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];

        end

        f0 = boundary_function(V0, bf_grid);
        f1 = boundary_function(V1, bf_grid);
    
    
        if (f0<=0 || f1<=0) 
            % the primal cell volume is unchanged
            S2_rows(ind) = indF;
            S2_cols(ind) = indF;
            S2_vals(ind) = get_Elength(f0, f1);

            P2_rows(ind) = ind;
            P2_cols(ind) = indF;
            ind = ind + 1;
        end

    end

end

S2 = sparse(S2_rows(1:ind-1), S2_cols(1:ind-1), S2_vals(1:ind-1), numF, numF);
P2 = sparse(P2_rows(1:ind-1), P2_cols(1:ind-1), P2_vals(1:ind-1), ind-1, numF);

end
