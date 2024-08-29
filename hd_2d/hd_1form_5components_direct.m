
function [omega_t_I, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components_direct(omega_ori, X, Y, bf_grid)

grid_size = size(X, 1);

% get all differential operators
D0 = get_D0(grid_size);
D1 = get_D1(grid_size);

[P0n, S0n] = get_P0S0(X, Y, bf_grid, "normal");%
[P1n, S1n] = get_P1S1(X, Y, bf_grid, "normal");%
[P2n, S2n] = get_P2S2(X, Y, bf_grid, "normal");%

D0n_int = P1n*D0*P0n'; %
D1n_int = P2n*D1*P1n';

S0n_int = P0n*S0n*P0n';
S1n_int = P1n*S1n*P1n'; %
S2n_int = P2n*S2n*P2n';

[P0t, S0t] = get_P0S0(X, Y, bf_grid, "tangential");
[P1t, S1t] = get_P1S1(X, Y, bf_grid, "tangential");%
[P2t, S2t] = get_P2S2(X, Y, bf_grid, "tangential");%

D1t_int = P2t*D1*P1t'; %
D0t_int = P1t*D0*P0t';

S0t_int = P0t*S0t*P0t';
S1t_int = P1t*S1t*P1t'; %
S2t_int = P2t*S2t*P2t'; %

S0t_int_inv = spdiags(1 ./ spdiags(S0t_int), 0, size(S0t_int, 1), size(S0t_int, 2));
S0n_int_inv = spdiags(1 ./ spdiags(S0n_int), 0, size(S0n_int, 1), size(S0n_int, 2));
S1t_int_inv = spdiags(1 ./ spdiags(S1t_int), 0, size(S1t_int, 1), size(S1t_int, 2));
S1n_int_inv = spdiags(1 ./ spdiags(S1n_int), 0, size(S1n_int, 1), size(S1n_int, 2));

L0n = D0n_int'*S1n_int*D0n_int;
L2t = S2t_int*D1t_int*S1t_int_inv*D1t_int'*S2t_int;

L1n = D1n_int'*S2n_int*D1n_int + S1n_int*D0n_int*S0n_int_inv*D0n_int'*S1n_int;
L1t = D1t_int'*S2t_int*D1t_int + S1t_int*D0t_int*S0t_int_inv*D0t_int'*S1t_int;

num_eigs = 5;

[H1t,eig_1t]=eigs(L1t, S1t_int, num_eigs, 'smallestabs'); 
[H1n,eig_1n]=eigs(L1n, S1n_int, num_eigs, 'smallestabs');

beta_1_t = sum(diag(eig_1t) <= 1e-10);
beta_1_n = sum(diag(eig_1n) <= 1e-10);
% disp(['beta 1 t: ' num2str(beta_1_t)])
% disp(['beta 1 n: ' num2str(beta_1_n)])

if beta_1_t ~= 0 
    H1t = H1t(:,1:beta_1_t);
else
    H1t = zeros(size(H1t, 1), 1);
end

if beta_1_n ~=0
    H1n = H1n(:,1:beta_1_n);
else 
    H1n = zeros(size(H1n, 1), 1);   
end

H1n = gram_schmidt(H1n, S1n_int);
H1t = gram_schmidt(H1t, S1t_int);
%% Hodge decomposition

omega_n = S1n_int_inv*P1n*omega_ori;
omega_t = P1t*omega_ori;

% 1st term
alpha_n = L0n \ (D0n_int'*S1n_int*omega_n);
dalpha_n = D0n_int*alpha_n;

% 2nd term
beta_t = L2t \ (S2t_int*D1t_int*omega_t);
deltabeta_t = S1t_int_inv*D1t_int'*S2t_int*beta_t;

% 3rd term
hn = H1n*(H1n'*S1n_int*omega_n);

% 4th term
ht = H1t*(H1t'*S1t_int*omega_t);

% projected back to the entire grid
omega_t_I = P1t'*omega_t;
dalpha_n_I = P1n'*S1n_int*dalpha_n; 
deltabeta_t_I = P1t'*deltabeta_t;
hn_I = P1n'*S1n_int*hn;
ht_I = P1t'*ht;

% 5th term
eta_I = P1t'*(omega_t - P1t*dalpha_n_I - deltabeta_t - P1t*hn_I - ht);

end
