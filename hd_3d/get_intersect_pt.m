
function V = get_intersect_pt(V0, V1, bf_grid)

% this function returns intersection points 

f0 = boundary_function(V0,bf_grid);
f1 = boundary_function(V1,bf_grid);

if f0 <= 0 && f1>0

    ratio = -f0/(f1-f0);
    V = (1-ratio)*V0 + ratio*V1;

elseif f0>0 && f1<=0

    ratio = -f1/(f0-f1);
    V = (1-ratio)*V1 + ratio*V0;
    
else % otherwise when both are the same sign return []
    V = [];
end
end