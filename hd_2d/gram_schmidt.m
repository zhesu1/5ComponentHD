
function H1_orthonormal = gram_schmidt(H1, S1)

% 1 gram-schemidt process
H1_orthonormal = zeros(size(H1));
H1_orthonormal(:,1) = H1(:,1);
if size(H1,2)>1
    for i = 2:size(H1,2)
        H1_orthonormal(:,i) = H1(:,i);
        for j = 1:(i-1)
            H1_orthonormal(:,i) = H1_orthonormal(:,i) - ...
                H1_orthonormal(:,j)*(H1(:,i)'*S1*H1_orthonormal(:,j))/...
                (H1_orthonormal(:,j)'*S1*H1_orthonormal(:,j));
        end
    end
end
% normalize
for i=1:size(H1,2)
    norm_H1n_t_orthonormal_i = sqrt(H1_orthonormal(:,i)'*S1*H1_orthonormal(:,i));
    if norm_H1n_t_orthonormal_i ~= 0
        H1_orthonormal(:,i) = H1_orthonormal(:,i)/norm_H1n_t_orthonormal_i;
    end
end
end
