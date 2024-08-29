
function omega = vflds_to_1forms(u, grid_size)
    % u: #V x 2
    
    numE = 2*(grid_size-1)*grid_size;

    omega = zeros(numE,1);

    for indE = 1:numE

        [i, j, direction] = indE_to_ijd(indE, grid_size);

        indV = ij_to_indV(i, j, grid_size);

        if direction == 0
            indV1 = ij_to_indV(i+1, j, grid_size);
            omega(indE) = (u(indV, 1) + u(indV1, 1))/2;
        else
            indV1 = ij_to_indV(i, j+1, grid_size);
            omega(indE) = (u(indV, 2) + u(indV1, 2))/2;
        end
    end
end



% function omega = vflds_to_1forms(u, grid_size)
%     % u: #V x 2
% 
%     numE = 2*(grid_size-1)*grid_size;
% 
%     omega = zeros(numE,1);
% 
%     for indE = 1:numE
% 
%         [i, j, direction] = indE_to_ijd(indE, grid_size);
% 
%         indV = ij_to_indV(i, j, grid_size);
% 
%         if direction == 0
%             indV1 = ij_to_indV(i+1, j, grid_size);
%             omega(indE) = (u(indV, 1) + u(indV1, 1))/2;
%         else
%             indV1 = ij_to_indV(i, j+1, grid_size);
%             omega(indE) = (u(indV, 2) + u(indV1, 2))/2;
%         end
%     end
% end