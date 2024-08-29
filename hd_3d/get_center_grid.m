
function [center_X, center_Y, center_Z] = get_center_grid(X, Y, Z)

center_X = (X(1:end-1, 1:end-1, 1:end-1) + X(2:end, 2:end, 2:end)) / 2;
center_Y = (Y(1:end-1, 1:end-1, 1:end-1) + Y(2:end, 2:end, 2:end)) / 2;
center_Z = (Z(1:end-1, 1:end-1, 1:end-1) + Z(2:end, 2:end, 2:end)) / 2;

end