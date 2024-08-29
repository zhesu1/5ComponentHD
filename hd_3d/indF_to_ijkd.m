
function [i, j, k, d] = indF_to_ijkd(indF, N)

d = floor((indF-1)/(N*(N-1)^2));

if d == 0 % z positive
    i = mod(indF-1, N-1) + 1;
    j = floor(mod(indF-1, (N-1)*(N-1)) / (N-1))+ 1;
    k = floor((indF-1)/((N-1)*(N-1))) + 1;
    
elseif d == 1 % y positive
    
    i = mod(indF-N*(N-1)^2-1, N-1) + 1;
    j = floor(mod(indF-N*(N-1)^2-1, N*(N-1))/(N-1)) + 1;
    k = floor((indF-N*(N-1)^2-1)/((N-1)*N)) + 1;

else % x positive
    
    i = mod(indF-2*N*(N-1)^2-1, N) + 1;
    j = floor(mod(indF-2*N*(N-1)^2-1, N*(N-1))/N) + 1;
    k = floor((indF-2*N*(N-1)^2-1)/((N-1)*N)) + 1;
    
end

end