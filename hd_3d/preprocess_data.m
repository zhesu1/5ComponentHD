function [X, Y, Z, U] = preprocess_data(x_grid, y_grid, z_grid, u_grid, v_grid, w_grid)

X = x_grid;
Y = y_grid;
Z = z_grid;

% rescale the data to make grid length 1
h = X(1,2,1)-X(1,1,1);
X = (X-X(1,1,1))/h;
Y = (Y-Y(1,1,1))/h;
Z = (Z-Z(1,1,1))/h;

u_grid = u_grid/h;
v_grid = v_grid/h;
w_grid = w_grid/h;
U = cat(4, u_grid, v_grid, w_grid);

end