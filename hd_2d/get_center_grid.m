
function [center_X, center_Y] = get_center_grid(X, Y)

center_X = (X(1:end-1, 1:end-1) + X(2:end, 2:end)) / 2;
center_Y = (Y(1:end-1, 1:end-1) + Y(2:end, 2:end)) / 2;
end