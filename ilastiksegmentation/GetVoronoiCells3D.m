%%
function [maskzcyto] = GetVoronoiCells3D(maskz,inuc,icyto)
%%
stats3d = regionprops(newmask,inuc(:,:,1:4),'Centroid');
xyz_fin = cat(1,stats3d.Centroid);

%----------------
%[V,C] = voronoin(X)
[V,C] = voronoin(xyz_fin);
% V rows represent the voronoi vertex
% C each cell has indices into V ( defines the verteces of the cell C{i}
% need to assign all verteces to  unique cells 

% -----------------algorythm 1:
% 1. get the coordinates of the verteces of the first cell C{1}
% 2. see how many nuclei coordinates lie within  hull, defined by this cell
% 3. check all cells C{} and all nuc coordinates
% 4. find cells with more than one nucleus 
% 5. get their convex hulls
% 6. assign 


for i=1:length(C), disp(C{i}), end % disply the cell verteces
%for i=1:length(V), disp(V(i,:,:)), end
%vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

 
for i=1:length(C), disp(C{i}), end 
% 
% in = inhull(testpts,xyz,tess,tol)
V(C{j},:); % verteces that define the j Vor cell
% test point: corrdinate of the nucleus????:
xyz_fin(:,:);        % check if the first nucleus is in several cells
for j=1:size(C,1)
  in(j) = inhull(xyz_fin(j,:),V(C{j}>1,:));

end
%----------another algorythm:

%1. get the vertces coordinates (remove infinity points)
%2. get the nuc coordinate
%3. find the interpoint distance matrix between the nuc coordinate and the
%rest of the verteces OR bte nuc and the verteces that are within the cell C{j} ,
%returned by the voronoin
%4. get the set of N (?) points with the min distance
%5. create a convex hull around that point :[K,V] = convhull(...) returns the convex hull K and the corresponding area/volume V bounded by K.
%5. alternative: use inhull to check if the nuc coordinate is in hull within those points
%6. do for all nuclei
%7. identify nuclei which are within several hulls ( use inhull function
%here (after running all the data)
%8. assign nucleus to unique hull (set of verteces) based on hull volume???

if V(1,1) == Inf                                                                  %step 1, don't really need this since Inf will not contr, to min distance
    V(1,:)=[];
end
vertcell = cell(1,size(xyz_fin,1));
for k=1:size(xyz_fin,1)                                                           %step 2
nuc1 = xyz_fin(k,:); 
for j=1:size(xyz_fin,1)
d = ipdm(nuc1,V,'result','structure','subset','SmallestFew','limit',8);        %step 3,4   V(C{k}>1,:)= coordinated of verteces of the first voronoi cell, excl Inf
% in = inhull(testpts,xyz,tess,tol)

d.columnindex; %row index of the vertex coordinate
d.distance;
xyz = V(d.columnindex,:);   %actual vertex coordinates of the closest N points
[K,Vol] = convhulln(xyz);    % each row of K is a triangle defined in terms of the point indices.      % step 5,6

% in(k) = inhull(nuc1,xyz);                                                                             % alternative step 5,6
% if in(k) == 1                       % get the hull coordinates
% vert{k} = xyz;
% end
end
end
figure,trisurf(K,xyz(:,1),xyz(:,2),xyz(:,3));

% figure,plot(xyz(K,1),xyz(K,2),'r*',xyz(:,1),xyz(:,2),'b-');
% hold on
% plot(nuc1(1),nuc1(2),'m*');


