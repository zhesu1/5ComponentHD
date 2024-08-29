

function [omega_t_I, dalpha_n_t_I, deltabeta_t_I,...
    hn_t_I, ht_I, eta_I] = hd_1form_5components_stepwise(omega_ori, X, Y, Z, bf_grid)

grid_size = size(X, 1);
% 
% fprintf('Get DEC operators \n')
D0 = get_D0(grid_size);
D1 = get_D1(grid_size);
D2 = get_D2(grid_size);

[P0n, S0n] = get_P0S0(X, Y, Z, bf_grid, "normal");%
[P1n, S1n] = get_P1S1(X, Y, Z, bf_grid, "normal");%
[P2n, S2n] = get_P2S2(X, Y, Z, bf_grid, "normal");%

D0n_int = P1n*D0*P0n'; %
D1n_int = P2n*D1*P1n';
S0n_int = P0n*S0n*P0n';
S1n_int = P1n*S1n*P1n'; %
S2n_int = P2n*S2n*P2n';

[P0t, S0t] = get_P0S0(X, Y, Z, bf_grid, "tangential");
[P1t, S1t] = get_P1S1(X, Y, Z, bf_grid, "tangential");%
[P2t, S2t] = get_P2S2(X, Y, Z, bf_grid, "tangential");%
[P3t, S3t] = get_P3S3(X, Y, Z, bf_grid, "tangential");

D0t_int = P1t*D0*P0t';
D1t_int = P2t*D1*P1t';
D2t_int = P3t*D2*P2t';

S0t_int = P0t*S0t*P0t';
S1t_int = P1t*S1t*P1t';%
S2t_int = P2t*S2t*P2t';%
S3t_int = P3t*S3t*P3t';

S0t_int_inv = spdiags(1./spdiags(S0t_int), 0, size(S0t_int, 1), size(S0t_int, 2));
S0n_int_inv = spdiags(1./ spdiags(S0n_int), 0, size(S0n_int, 1), size(S0n_int, 2));
S1t_int_inv = spdiags(1./spdiags(S1t_int), 0, size(S1t_int, 1), size(S1t_int, 2));
% S1n_int_inv = spdiags(1./spdiags(S1n_int), 0, size(S1n_int, 1), size(S1n_int, 2));

% L0n = D0n_int'*S1n_int*D0n_int;
L0t = D0t_int'*S1t_int*D0t_int;
L1t = D1t_int'*S2t_int*D1t_int + S1t_int*D0t_int*S0t_int_inv*D0t_int'*S1t_int;
L2t = D2t_int'*S3t_int*D2t_int + S2t_int*D1t_int*S1t_int_inv*D1t_int'*S2t_int;

L1n = D1n_int'*S2n_int*D1n_int + S1n_int*D0n_int*S0n_int_inv*D0n_int'*S1n_int;

num_eigs = 5;

[H1t,eig_1t]=eigs(L1t, S1t_int, num_eigs, 'smallestabs'); 
[H1n,eig_1n]=eigs(L1n, S1n_int, num_eigs, 'smallestabs');

beta_1_t = sum(diag(eig_1t) <= 1e-10);
beta_2_n = sum(diag(eig_1n) <= 1e-10);
% disp(['beta 1: ' num2str(beta_1_t)])
% disp(['beta 2: ' num2str(beta_2_n)])

% if beta_1_t ~= 0 
%     H1t = H1t(:,1:beta_1_t);
% else
%     H1t = zeros(size(H1t, 1), 1);
% end

if beta_2_n ~=0
    H1n = H1n(:,1:beta_2_n);
else 
    H1n = zeros(size(H1n, 1), 1);   
end

Proj_L0t = [sparse(size(S0t_int,1)-1, 1) speye(size(S0t_int,1)-1)];
bar_L0t = Proj_L0t*L0t*Proj_L0t';

% add more cases
if beta_2_n == 1
    Proj_L2t = [sparse(size(S2t_int,1)-1, 1) speye(size(S2t_int,1)-1)];
else
    Proj_L2t = speye(size(S2t_int,1));
end
bar_L2t = Proj_L2t*L2t*Proj_L2t';

%% 1 decompose omega_t = d\alpha_t + \delta\beta_t + ht

omega_t = P1t*omega_ori;

% 1 dalpha_t
alpha_t_projected = bar_L0t \ (Proj_L0t*D0t_int'*S1t_int*omega_t);
alpha_t = Proj_L0t'*alpha_t_projected;
dalpha_t = D0t_int*alpha_t;

% 2 delta_beta_t
beta_t_projected = bar_L2t \ (Proj_L2t*S2t_int*D1t_int*omega_t);
beta_t = Proj_L2t'*beta_t_projected;
deltabeta_t = S1t_int_inv*D1t_int'*S2t_int*beta_t;

% h_t
ht = (omega_t - dalpha_t - deltabeta_t);

% projected back to the entire grid
omega_t_I = P1t'*omega_t;
dalpha_t_I = P1t'*dalpha_t;
deltabeta_t_I = P1t'*deltabeta_t;
ht_I = P1t'*ht;

%% 2 decompose d\alpha_t = d\alpha_n + eta + hn

% fprintf('Build linear system \n')

[P1E, P0E] = get_support_E(X, Y, Z, bf_grid);

P_E2T = P0t*P0E';
D0E = P1E*D0*P0E';

L0t_E = P_E2T'*L0t*P_E2T;
bar_L = P0n*P0E'*(D0E'*D0E);

Proj_E = [sparse(size(L0t_E,1)-1, 1) speye(size(L0t_E,1)-1)];

L0t_E_proj = Proj_E*L0t_E*Proj_E';
bar_L_proj = bar_L*Proj_E';


A = [L0t_E_proj bar_L_proj'; 
    bar_L_proj sparse(size(bar_L_proj, 1), size(bar_L_proj, 1))];

b = [Proj_E*P_E2T'*L0t*alpha_t; 
    sparse(size(bar_L_proj, 1), 1)];

A = A + 1e-15*speye(size(A));
x = A\b;

hn_plus_eta = D0t_int*P_E2T*Proj_E'*x(1:size(L0t_E_proj,1), 1);

% 1 dalpha_n
dalpha_n_t = dalpha_t - hn_plus_eta;

%% get hn 

% get a basis for h_n on tangential support
T_N2T = P1t*P1n'*S1n_int;
bar_L_H1n = D0t_int'*S1t_int*T_N2T*H1n;
H1n_T_converted = D0t_int*Proj_L0t'*(bar_L0t\(Proj_L0t*bar_L_H1n));
H1n_T = sparse(size(H1n_T_converted, 1), size(H1n_T_converted, 2));

for ind_B_Hn = 1:size(H1n, 2)
    b = [Proj_E*P_E2T'*D0t_int'*S1t_int*H1n_T_converted(:, ind_B_Hn); 
        sparse(size(bar_L_proj, 1), 1)];
    x = A\b;
    H1n_T(:,ind_B_Hn) = D0t_int*P_E2T*Proj_E'*x(1:size(L0t_E_proj,1), 1);
end

H1n_orthonormal_T = gram_schmidt(H1n_T, S1t_int);

% 2 compute hn
hn_t = H1n_orthonormal_T*(H1n_orthonormal_T'*S1t_int*hn_plus_eta);

% 3 compute eta
eta_t = hn_plus_eta - hn_t;

% projected back to the entire grid
dalpha_n_t_I = P1t'*dalpha_n_t;
hn_t_I = P1t'*hn_t;
eta_I = P1t'*eta_t;

end
