
function indC = ijk_to_indC(i,j,k,N)

indC = i + (N-1)*(j-1) + (N-1)*(N-1)*(k-1);

end