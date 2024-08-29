
function h = plot_mesh(vertex,face)

h = patch('vertices',vertex,'faces',face, ...
    'FaceVertexCData',zeros(size(vertex,1),1), ...
    'FaceColor','interp', 'facealpha',0.3);

lighting phong;
camproj('perspective');
axis square; 
axis off;
axis tight;
axis equal;