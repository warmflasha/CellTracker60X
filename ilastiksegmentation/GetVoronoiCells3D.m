%%
GetVoronoiCells3D(maskz,inuc)

stats3d = regionprops(maskz,inuc(:,:,1:4),'Centroid');
xyz_fin = cat(1,stats3d.Centroid);


%[V,C] = voronoin(X)
[V,C] = voronoin(xyz_fin);
% V rows represent the voronoi vertex
% C each cell has indices into V ( defines the verteces of the cell C{i}
% need to assign all verteces to  unique cells 
% the number of the unique cells will not be equal to the number of
% verteces, obv.

for i=1:length(C), disp(C{i}), end % disply the cell verteces

%vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

% inhull