
function [P1, S1] = get_P1S1(X, Y, Z, bf_grid, key)

grid_size = size(X, 1);

if nargin < 5
    key = "normal";
end

epsilon = 1e-5;
numE = 3*(grid_size-1)*grid_size*grid_size;

P1_rows = zeros(numE, 1);
P1_cols = zeros(numE, 1);
P1_vals = ones(numE, 1);

S1_rows = zeros(numE, 1);
S1_cols = zeros(numE, 1);
S1_vals = zeros(numE, 1);

ind = 1;
for indE = 1: numE

    [i, j, k, d] = indE_to_ijkd(indE, grid_size);

    if key == "normal" % consider primal edges

        if d==0
    
            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j,i+1,k), Y(j,i+1,k), Z(j,i+1,k)];
           
        elseif d==1
    
            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j+1,i,k), Y(j+1,i,k), Z(j+1,i,k)];
    
        else
    
            V0 = [X(j,i,k), Y(j,i,k), Z(j,i,k)];
            V1 = [X(j,i,k+1), Y(j,i,k+1), Z(j,i,k+1)];
    
        end
    
         f0 = boundary_function(V0, bf_grid);
         f1 = boundary_function(V1, bf_grid);
    
    
        if (f0<=0 || f1<=0) % normal
            % the dual cell volume is unchanged
            S1_rows(ind) = indE;
            S1_cols(ind) = indE;
            S1_vals(ind) = 1/max(get_Elength(f0, f1), epsilon);

            P1_rows(ind) = ind;
            P1_cols(ind) = indE;
            ind = ind + 1;

        end

    elseif key == "tangential" % consider dual faces

        if d==0
    
            V0 = [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)-1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
            V2 = [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
            V3 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
           
        elseif d==1
    
            V0 = [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)-1/2];
            V2 = [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
            V3 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
    
        else
    
            V0 = [X(j,i,k)-1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
            V1 = [X(j,i,k)+1/2, Y(j,i,k)-1/2, Z(j,i,k)+1/2];
            V2 = [X(j,i,k)-1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
            V3 = [X(j,i,k)+1/2, Y(j,i,k)+1/2, Z(j,i,k)+1/2];
    
        end
    
        f0 = boundary_function(V0, bf_grid);
        f1 = boundary_function(V1, bf_grid);
        f2 = boundary_function(V2, bf_grid);
        f3 = boundary_function(V3, bf_grid);

        if (f0<=0 || f1<=0 || f2<=0 || f3<=0) 
            % the primal cell volume is unchanged
            S1_rows(ind) = indE;
            S1_cols(ind) = indE;
            S1_vals(ind) = get_Farea(f0, f1, f2, f3);

            P1_rows(ind) = ind;
            P1_cols(ind) = indE;
            ind = ind + 1;

        end

    end
end

S1 = sparse(S1_rows(1:ind-1), S1_cols(1:ind-1), S1_vals(1:ind-1), numE, numE);
P1 = sparse(P1_rows(1:ind-1), P1_cols(1:ind-1), P1_vals(1:ind-1), ind-1, numE);

end
