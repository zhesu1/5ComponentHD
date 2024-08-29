
function b_mask = boundary_mask(x_grid, y_grid, bf_grid)

    f = zeros(size(x_grid));
    for i=1:size(x_grid,1)
        for j=1:size(x_grid,1)
            f(j,i) = boundary_function(x_grid(j,i), y_grid(j,i), bf_grid);
        end
    end
    b_mask = f <= 0;
end
