

close all
clear, clc 

addpath('../hd_3d')
addpath('plt/')

%% load data
load('data/g1_1c_remeshed2.mat')

%% pertube the lvf for numerical stability such that values on grid points and centers are nonzeros

epsilon = 1e-5;
bf_grid(0 <= bf_grid & bf_grid<epsilon) = epsilon;
bf_grid(0 > bf_grid & bf_grid>-epsilon) = -epsilon;

%% Hodge decomposition

key = "direct";
% key = "stepwise";

[omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components(omega_ori, X, Y, Z, bf_grid, key);

%% plot
set_view = [-180,-50];

get_hd_streamlines_sdf(omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I, X, Y, Z, vert, ...
    tri, bf_grid, set_view)

