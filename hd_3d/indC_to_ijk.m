
function [i, j, k] = indC_to_ijk(indC, N)

i = mod(indC-1, N-1) + 1;
j = floor(mod(indC-1, (N-1)*(N-1)) / (N-1)) + 1;
k = floor((indC-1) / ((N-1)*(N-1))) + 1;

end