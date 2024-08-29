function [X, Y, U] = preprocess_data(x_grid, y_grid, u_grid, v_grid)

X = x_grid;
Y = y_grid;

% rescale the data to make grid length 1
h = X(1,2)-X(1,1);
X = (X-X(1,1))/h;
Y = (Y-Y(1,1))/h;

u_grid = u_grid/h;
v_grid = v_grid/h;
U = cat(3, u_grid, v_grid);

end
