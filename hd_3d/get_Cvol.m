

function Cvol = get_Cvol(VC, bf_VC, bf_grid)

% find the vertex indices where f are negative
% inds = find(bf_VC < 0);
inds = find(sign(bf_VC) == -1);
% disp(indC)

if length(find(bf_VC >= 0)) == 8 % function values on vertices all positive - outside the boundary
    % here we also consider the case that =0 to remove the degeneate case
    % for the convhull, which gives an error messge: Not enough unique points specified.

    Cvol = 0;

elseif length(find(bf_VC <= 0)) == 8

    Cvol = 1; % all are interior vertices

else
    
    ConvhVlist = get_v_convh(VC, inds, bf_grid);

    % rank(ConvhVlist)<=2 || size(unique(ConvhVlist,'rows'),1) <= 3
    if rank(ConvhVlist - ConvhVlist(1,:))<=2 % check if the interior points are in the same plane
        Cvol = 0;
    else
        [~,Cvol] = convhull(ConvhVlist);
    end

end


end


function ConvhVlist = get_v_convh(VC, inds, bf_grid)
ConvhVlist = [];
for i = 1: size(inds)
    [ind1, ind2, ind3] = return_adjindV_Cell(inds(i));
    V = VC(inds(i),:);
    V1 = get_intersect_pt(V, VC(ind1,:), bf_grid);
    V2 = get_intersect_pt(V, VC(ind2,:), bf_grid);
    V3 = get_intersect_pt(V, VC(ind3,:), bf_grid);
    ConvhVlist = [ConvhVlist; V; V1; V2; V3];
end
end


function [ind1, ind2, ind3] = return_adjindV_Cell(ind0)

% this function gives the indices of three points which connect by edges to the
% point with index ind0

if ind0 < 5
    ind1 = mod(ind0,4) + 1;
    ind2 = mod(ind0+2,4) + 1;
    ind3 = ind0+4;
else
    ind1 = mod(ind0,4) + 5;
    ind2 = mod(ind0+2,4) + 5;
    ind3 = mod(ind0+3,8) + 1;
end

end


