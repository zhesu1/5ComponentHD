
close all
clear
clc

addpath('../hd_2d')
addpath('plt/')

%% load data
load('data/bunny_2d.mat')

%% pertube the lvf for numerical stability such that values on grid points and centers are nonzeros

epsilon = 1e-5;
bf_grid(0 <= bf_grid & bf_grid<epsilon) = epsilon;
bf_grid(0 > bf_grid & bf_grid>-epsilon) = -epsilon;

%% Hodge decomposition

% key = "direct";
key = "stepwise";

[omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components(omega_ori, X, Y, bf_grid, key);

%% plot

get_hd_buildin_streamlines(omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I, X, Y, bf_grid)
