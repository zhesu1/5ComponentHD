
function [i, j, k, d] = indE_to_ijkd(indE, grid_size)

d = floor((indE-1)/((grid_size-1) * grid_size * grid_size));

if d==0 % x positive
    i = mod(indE-1, grid_size-1)+1;
    j = floor(mod(indE-1, grid_size*(grid_size-1)) / (grid_size-1))+1;
    k = floor((indE-1) / (grid_size*(grid_size-1)))+1;
    
elseif d==1 % y positive
    i = mod(indE-(grid_size-1)*grid_size*grid_size-1, grid_size)+1;
    j = floor(mod(indE-(grid_size-1)*grid_size*grid_size-1, (grid_size-1)*grid_size) / grid_size)+1;
    k = floor((indE-(grid_size-1)*grid_size*grid_size-1) / ((grid_size-1)*grid_size))+1;
    
else % z positive
    i = mod(indE-2*(grid_size-1)*grid_size*grid_size-1, grid_size)+1;
    j = floor(mod(indE-2*(grid_size-1)*grid_size*grid_size-1, grid_size*grid_size) / grid_size)+1;
    k = floor((indE-2*(grid_size-1)*grid_size*grid_size-1) / (grid_size*grid_size))+1;
end

end