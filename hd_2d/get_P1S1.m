
function [P1, S1] = get_P1S1(X, Y, bf_grid, key)

grid_size = size(X, 1);

if nargin < 4
    key = "normal";
end

epsilon = 1e-5;
numE = 2*(grid_size-1)*grid_size;

P1_rows = zeros(numE, 1);
P1_cols = zeros(numE, 1);
P1_vals = ones(numE, 1);

S1_rows = zeros(numE, 1);
S1_cols = zeros(numE, 1);
S1_vals = zeros(numE, 1);

ind = 1;
for indE = 1: numE

    [i, j, d] = indE_to_ijd(indE, grid_size);

    if key == "normal" 
        if d==0
            f0 = boundary_function(X(j, i), Y(j,i), bf_grid);
            f1 = boundary_function(X(j,i+1), Y(j,i+1), bf_grid);
        else
            f0 = boundary_function(X(j,i), Y(j, i), bf_grid);
            f1 = boundary_function(X(j+1, i), Y(j+1, i), bf_grid);
        end

        if (f0<=0 || f1<=0)

            % normal
            % the dual cell volume is unchanged
            S1_rows(ind) = indE;
            S1_cols(ind) = indE;
            S1_vals(ind) = 1/max(get_primal_Elength(f0, f1), epsilon);

            P1_rows(ind) = ind;
            P1_cols(ind) = indE;
            ind = ind + 1;

        end

    elseif key == "tangential" 

        if d==0
            f0 = boundary_function(X(j,i)+1/2, Y(j,i)-1/2, bf_grid);
            f1 = boundary_function(X(j,i)+1/2, Y(j,i)+1/2, bf_grid);
        else
            f0 = boundary_function(X(j,i)-1/2, Y(j,i)+1/2, bf_grid);
            f1 = boundary_function(X(j,i)+1/2, Y(j,i)+1/2, bf_grid);
        end
        
        if(f0<=0 || f1<=0)

            % tangential
            % the primal cell volume is unchanged
            S1_rows(ind) = indE;
            S1_cols(ind) = indE;
            S1_vals(ind) = get_dual_Elength(f0, f1);

            P1_rows(ind) = ind;
            P1_cols(ind) = indE;
            ind = ind + 1;

        end
    end

end

S1 = sparse(S1_rows(1:ind-1), S1_cols(1:ind-1), S1_vals(1:ind-1), numE, numE);
P1 = sparse(P1_rows(1:ind-1), P1_cols(1:ind-1), P1_vals(1:ind-1), ind-1, numE);

end

