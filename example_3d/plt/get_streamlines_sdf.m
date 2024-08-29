
function verts_all = get_streamlines_sdf(X, Y, Z, vf, bf_grid, vert, k, vf_name, epsilon)

epsilon_streamlines = epsilon; 
epsilon_seedPts = epsilon;


[center_X, center_Y, center_Z] = get_center_grid(X, Y, Z);


[startX,startY,startZ] = meshgrid(min(vert(:,1)):k:max(vert(:,1)),...
    min(vert(:,2)):k:max(vert(:,2)),...
    min(vert(:,3)):k:max(vert(:,3)));


bf_grid_start = interp3(X, Y, Z, bf_grid, startX, startY, startZ);
interior_mask_start = bf_grid_start < epsilon_seedPts;


% positive direction
verts0 = stream3(center_X, center_Y, center_Z, ...
    vf(:,:,:,1), vf(:,:,:,2), vf(:,:,:,3),...
    startX(interior_mask_start), startY(interior_mask_start), startZ(interior_mask_start),...
    [0.1 500]);

% negative direction
verts1 = stream3(center_X, center_Y, center_Z, ...
    -vf(:,:,:,1), -vf(:,:,:,2), -vf(:,:,:,3),...
    startX(interior_mask_start), startY(interior_mask_start), startZ(interior_mask_start),...
    [0.1 500]);


verts_all = cell(1,0);
for j = 1:length(verts0)
    % disp(j)
    pts_line = [flip(verts1{j}); verts0{j}(2:end,:)]; % remove duplicate seed points
    f_pts_line = interp3(X, Y, Z, bf_grid, pts_line(:,1), pts_line(:,2), pts_line(:,3));

    neg_indices = find(f_pts_line <= epsilon_streamlines);
    pts_inside_consecutive = all(diff(neg_indices) == 1); % check if the interior pts on streamlines are consecutive, otherwise remove the streamlines
    if pts_inside_consecutive
        pts_line = pts_line(f_pts_line < epsilon_streamlines,:); % keep only the part of streamlines inside the boundary
        l_pts_line = sum(sqrt(sum(diff(pts_line).^2, 2))); % length of the streamline
        if l_pts_line > 1
            verts_all = [verts_all, pts_line];
        end
    else
        breaks0 = find(diff(neg_indices) ~= 1) + 1; % find the breaks
        breaks = [1, breaks0', numel(neg_indices) + 1]; 

        for i = 1:numel(breaks)-1 % split into interior segments
            verts_i = pts_line(neg_indices(breaks(i):breaks(i + 1) - 1), :);
            l_verts_i = sum(sqrt(sum(diff(verts_i).^2, 2))); % length of the streamline

            if size(verts_i, 1)>10 && l_verts_i > 1 % add only long segments
                verts_all = [verts_all, verts_i];
            end
        end
    end

end

lineobj = streamline(verts_all);
set(lineobj, 'LineWidth', .5);

title(vf_name)
hold on;

daspect([1 1 1])

end


