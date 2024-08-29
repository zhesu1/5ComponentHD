

function [i, j, direction] = indE_to_ijd(indE, grid_size)

   direction = floor((indE-1)/((grid_size-1)*grid_size));

   if direction == 0 % horizontal
       i = mod(indE-1, grid_size-1) + 1;
       j = floor((indE-1)/(grid_size-1)) + 1;
   else % vertical
       i = mod(indE-(grid_size-1)*grid_size-1, grid_size) + 1; 
       j = floor((indE-(grid_size-1)*grid_size-1)/grid_size) + 1;
   end


end

% % vectorized
% function [is, js, ds] = indE_to_ijd(indEs, grid_size)
% 
% ds = floor((indEs-1)/((grid_size-1)*grid_size));
% 
% ind_d0s = ds == 0;
% ind_d1s = ds == 1;
% 
% is = zeros(size(indEs));
% js = zeros(size(indEs));
% 
% is(ind_d0s) = mod(indEs(ind_d0s)-1, grid_size-1) + 1;
% js(ind_d0s) = floor((indEs(ind_d0s)-1)/(grid_size-1)) + 1;
% 
% is(ind_d1s) = mod(indEs(ind_d1s)-(grid_size-1)*grid_size-1, grid_size) + 1;
% js(ind_d1s) = floor((indEs(ind_d1s)-(grid_size-1)*grid_size-1)/grid_size) + 1;
% 
% end

