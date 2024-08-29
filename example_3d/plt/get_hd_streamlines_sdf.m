
function get_hd_streamlines_sdf(omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I, X, Y, Z, ...
    vert, tri, bf_grid, view0)

grid_size = size(X,1);

vecfCell_grid_omega = one_forms_to_vecfCell(omega,  grid_size);

vecfCell_grid_dalpha_n = one_forms_to_vecfCell(dalpha_n_I, grid_size);
vecfCell_grid_deltabeta_t = one_forms_to_vecfCell(deltabeta_t_I, grid_size);
vecfCell_grid_hn = one_forms_to_vecfCell(hn_I, grid_size);
vecfCell_grid_ht = one_forms_to_vecfCell(ht_I, grid_size);
vecfCell_grid_eta = one_forms_to_vecfCell(eta_I, grid_size);

%%

[center_X, center_Y, center_Z] = get_center_grid(X, Y, Z);
bf_grid_c = interp3(X, Y, Z, bf_grid, center_X, center_Y, center_Z);

epsilon1 = -0.1;

interior_mask_c = bf_grid_c < -epsilon1; %  set -0.1 to remove the artifacts of the vectors near the boundary

% centralize the mesh 
X = X - grid_size/2;
Y = Y - grid_size/2;
Z = Z - grid_size/2;
vert = vert - grid_size/2;

%% get masked vfs

vecfCell_grid_omega = get_masked_vfld(vecfCell_grid_omega, interior_mask_c);
vecfCell_grid_dalpha_n = get_masked_vfld(vecfCell_grid_dalpha_n, interior_mask_c);
vecfCell_grid_deltabeta_t = get_masked_vfld(vecfCell_grid_deltabeta_t, interior_mask_c);
vecfCell_grid_hn = get_masked_vfld(vecfCell_grid_hn, interior_mask_c);
vecfCell_grid_ht = get_masked_vfld(vecfCell_grid_ht, interior_mask_c);
vecfCell_grid_eta = get_masked_vfld(vecfCell_grid_eta, interior_mask_c);

%% plot

k = 3;

% the original vector field
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_omega, ...
    bf_grid, vert, k, "ori", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
title("Original")
cameratoolbar('SetCoordSys','x')
camlight('headlight')
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


% dalpha_n
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_dalpha_n, ...
    bf_grid, vert, k, "dalpha n", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
cameratoolbar('SetCoordSys','x')
camlight('headlight')
title("d\alpha_n")
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


% delta_beta_t
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_deltabeta_t, ...
    bf_grid, vert, 5, "deltabeta t", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
cameratoolbar('SetCoordSys','x')
camlight('headlight')
title("\delta\beta_t")
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


% hn
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_hn, ...
    bf_grid, vert, k, "hn", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
title("h_n")
cameratoolbar('SetCoordSys','x')
camlight('headlight')
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


% ht
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_ht, ...
    bf_grid, vert, 4, "ht", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
title("h_t")
cameratoolbar('SetCoordSys','x')
camlight('headlight')
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


% eta
figure;set(gcf,'color','w');
get_streamlines_sdf(X, Y, Z, vecfCell_grid_eta, ...
    bf_grid, vert, k, "eta", epsilon1);
plot_mesh(vert, tri)
colormap('sky');
shading interp;
view(view0)
title("\eta")
cameratoolbar('SetCoordSys','x')
camlight('headlight')
xlim([-grid_size/2, grid_size/2])
ylim([-grid_size/2, grid_size/2])
zlim([-grid_size/2, grid_size/2])


end

