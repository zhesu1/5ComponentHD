

function get_hd_buildin_streamlines(omega, dalpha_n_I, deltabeta_t_I, hn_I, ht_I, eta_I, X, Y, ...
    bf_grid)

grid_size = size(X,1);

vecfFace_grid_omega = one_forms_to_vecfFace(omega, grid_size);

vecfFace_grid_dalpha_n = one_forms_to_vecfFace(dalpha_n_I, grid_size);
vecfFace_grid_deltabeta_t = one_forms_to_vecfFace(deltabeta_t_I, grid_size);
vecfFace_grid_hn = one_forms_to_vecfFace(hn_I, grid_size);
vecfFace_grid_ht = one_forms_to_vecfFace(ht_I, grid_size);
vecfFace_grid_eta = one_forms_to_vecfFace(eta_I, grid_size);

%%
[center_X, center_Y] = get_center_grid(X, Y);
bf_grid_c = interp2(X, Y, bf_grid, center_X, center_Y);
interior_mask_c = bf_grid_c < 0;

%%  get masked vfs
vecfFace_grid_omega = get_masked_vfld(vecfFace_grid_omega, interior_mask_c);
vecfFace_grid_dalpha_n = get_masked_vfld(vecfFace_grid_dalpha_n, interior_mask_c);
vecfFace_grid_deltabeta_t = get_masked_vfld(vecfFace_grid_deltabeta_t, interior_mask_c);
vecfFace_grid_hn = get_masked_vfld(vecfFace_grid_hn, interior_mask_c);
vecfFace_grid_ht = get_masked_vfld(vecfFace_grid_ht, interior_mask_c);
vecfFace_grid_eta = get_masked_vfld(vecfFace_grid_eta, interior_mask_c);

%% plot

X = X + 1;
Y = Y + 1;

density = 3;

% the original vector field
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0)
set(gca,'YDir','reverse')
title("Original")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_omega, bf_grid, density, "Original");
axis off
axis equal


% dalpha_n
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0)
set(gca,'YDir','reverse')
title("d\alpha_n")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_dalpha_n, bf_grid, density, "curl free");
axis off
axis equal


% delta_beta_t
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0)
set(gca,'YDir','reverse')
title("\delta\beta_t")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_deltabeta_t, bf_grid, density, "div free");
axis off
axis equal


% hn
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0)
set(gca,'YDir','reverse')
title("h_n")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_hn, bf_grid, density, "hn");
axis off
axis equal


% ht
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0)
set(gca,'YDir','reverse')
title("h_t")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_ht, bf_grid, density, "ht");
axis off
axis equal


% eta
figure;set(gcf,'color','w');
hold on;
contourf(bf_grid, [0 0],"FaceAlpha",0);
set(gca,'YDir','reverse')
title("\eta")
get_interior_buildin_streamlines(X, Y, vecfFace_grid_eta, bf_grid, density, "eta");
axis off
axis equal

end










