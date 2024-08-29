
function get_interior_buildin_streamlines(X, Y, vf, bf_grid, k, vf_name)

[center_X, center_Y] = get_center_grid(X, Y);

[verts,averts] = streamslice(center_X, center_Y, ...
    vf(:,:,1), vf(:,:,2), k);

epsilon = -3; 

% remove exterior points on streamlines
verts_all = cell(1,0);
for i=1:length(verts)
    pts_line = verts{i};
    if size(pts_line,1)~=0
        f_pts_line = interp2(X, Y, bf_grid, pts_line(:,1),  pts_line(:,2));
        neg_indices = find(f_pts_line <= epsilon); 
        pts_inside_consecutive = all(diff(neg_indices) == 1); % check if the interior pts on streamlines are consecutive, otherwise remove the streamlines
        if pts_inside_consecutive
            pts_line = pts_line(f_pts_line < epsilon,:); % keep only the part of streamlines inside the boundary
            verts_all = [verts_all, pts_line];
        else
            breaks0 = find(diff(neg_indices) ~= 1) + 1; % find the breaks
            breaks = [1, breaks0', numel(neg_indices) + 1]; 
    
            for i = 1:numel(breaks)-1 % split into interior segments
                verts_i = pts_line(neg_indices(breaks(i):breaks(i + 1) - 1), :);
                if size(verts_i, 1)>10 % add only long segments
                    verts_all = [verts_all, verts_i];
                end
            end
        end
    end
end

% remove exterior arrows
averts_all = cell(1,0);
for i=1:length(averts)
    pts_arrow = averts{i};
    f_pts_arrow = interp2(X, Y, bf_grid, pts_arrow(:,1),  pts_arrow(:,2));
    neg_indices = find(f_pts_arrow <= 0);
    if length(neg_indices)==3
        averts_all = [averts_all, pts_arrow];
    end
end

% avoid the stream2 error that U, V must all be of size 2x2 and more
if length(verts_all)>1
    lineobj = streamline([verts_all, averts_all]);
    set(lineobj, 'LineWidth', .7);
end
% title(vf_name)
% plot(averts{1}(:,1), averts{1}(:,2), '+')
end


