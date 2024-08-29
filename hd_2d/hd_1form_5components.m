
function [omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components(omega_ori, X, Y, bf_grid, key)

if key == "direct"
    [omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components_direct(omega_ori, X, Y, bf_grid);
elseif key == "stepwise"
    [omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I] = hd_1form_5components_stepwise(omega_ori, X, Y, bf_grid);
end

end