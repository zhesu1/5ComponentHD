
function indF = ijkd_to_indF(i, j, k, d, N)

if d == 0 % z positive
    
    indF = i + (N-1)*(j-1) + (N-1)*(N-1)*(k-1);
     
elseif d==1 % y positive
   
    indF = i + (N-1)*(j-1) + (N-1)*N*(k-1) + N*(N-1)^2;

else % x positive
    
    indF = i + N*(j-1) + (N-1)*N*(k-1) + 2*N*(N-1)^2;
    
end
end
