

function f = boundary_function(V, bf_grid)

x = V(1);
y = V(2);
z = V(3);

epsilon = 1e-5;
grid_size = size(bf_grid,1);


if abs(mod(x,1)) < epsilon && abs(mod(y, 1)) < epsilon && abs(mod(z, 1)) < epsilon
    f = bf_grid(round(y)+1, round(x)+1, round(z)+1);
else
    x0 = floor(x);
    y0 = floor(y);
    z0 = floor(z);

    if x0<0 || y0<0 || z0<0 || x0>=grid_size-1 || y0>=grid_size-1 || z0>=grid_size-1
        f = inf;
        return;
    end

    % disp([x0, y0])
    f = (bf_grid(y0+1, x0+1, z0+1) + bf_grid(y0+1, x0+2, z0+1) +...
         bf_grid(y0+2, x0+1, z0+1) + bf_grid(y0+2, x0+2, z0+1) +...
         bf_grid(y0+1, x0+1, z0+2) + bf_grid(y0+1, x0+2, z0+2) +...
         bf_grid(y0+2, x0+1, z0+2) + bf_grid(y0+2, x0+2, z0+2))/8;
end

f(f >= 0 & f<epsilon) = epsilon;
f(f < 0 & f>-epsilon) = -epsilon;
end
