%%
function [maskzcyto] = GetVoronoiCells3D(maskz,inuc,icyto)
%%
 stats3d = regionprops(newmask,inuc(:,:,1:4),'Centroid');
% xyz_fin = cat(1,stats3d.Centroid);

clear nucmaskall
clear dist
clear dist_all
% mask1test = (newmask(:,:,1)==2);
% dist = bwdist(mask1test);
nelem = max(max(max(newmask))); %number of elements in all masks
nelemvect = 1:nelem;
min_ind = cell(1,size(newmask,3));
Inew = zeros(1024,1024,size(newmask,3));
%j=2;
range =cell(1,nelem);
ints =cell(1,nelem);

% 
for ii=1%:size(newmask,3)
    nucmaskall = newmask(:,:,ii);
    for j=1:nelem
        nuc_submask(:,:,j) = nucmaskall==nelemvect(j);
        if sum(sum(nuc_submask(:,:,j)))>0
            nuc_submask1(:,:,j)=nuc_submask(:,:,j);
        end                                             %at this point have the separate nuclei masks with neighboring cells, from which to take bwdist
    end
end
 
for k=1:size(nuc_submask1,3)
    dist(:,:,:,k) = bwdist(nuc_submask1(:,:,k));%
end

 dist_all = cat(3,dist);
 
 [~,min_ind] = min(dist_all,[],3);% dist transform number (~ row ) with the pixels having min distance to the nonzero object
 
 
 

%%
for k=1:size(nuc_submask1,3)
    figure, imshow(nuc_submask1(:,:,k),[]);
end
%%

for k=1:size(nuc_submask1,3)
    figure, imshow(dist_all(:,:,k),[]);
end

%%
for k=1:4
    figure, imshow(nucmaskall(:,:,k),[]);
end

%%
for k=1:4
    figure, imshow(Inew(:,:,k),[]);
end

%%

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
% 6. assign to different cells based on area?


for i=1:length(C), disp(C{i}), end % disply the cell verteces
%for i=1:length(V), disp(V(i,:,:)), end
%vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

 
% %----------another algorythm:
% 
% % 1. get the vertces coordinates (remove infinity points ?)
% % 2. get the nuc coordinate
% % 3. find the interpoint distance matrix between the nuc coordinate and the
% %   rest of the verteces (OR btw nuc and the verteces that are within the cell C{j} ,
% %   returned by the voronoin)
% % 4. get the set of N (?)  min distance points (vertices that are closest to
% %   nuc1)
% % 5. create a convex hull/boundary (?) around those points with the nuclear coordinate :[K,V] = convhull(...) returns the convex hull K and the corresponding area/volume V bounded by K.
% %5.  alternative: use inhull to check if the nuc coordinate is in hull within those points
% % 6. Do for all nuclei
% % 7. identify nuclei which are within several hulls ( for which the 3d hulls
% %   are intersecting
% % 8. assign nucleus to unique hull (set of verteces) based on hull volume???
% 
% if V(1,1) == Inf                                                                  %step 1, don't really need this since Inf will not contr, to min distance
%     V(1,:)=[];
% end
% vert = cell(1,size(xyz_fin,1));
% for k=1:size(xyz_fin,1)                                                           %step 2
% nuc1 = xyz_fin(k,:); 
%     
% d = ipdm(nuc1,V,'result','structure','subset','SmallestFew','limit',6);        %step 3,4   V(C{k}>1,:)= coordinated of verteces of the first voronoi cell, excl Inf
% 
% xyz = V(d.columnindex,:);    % actual vertex coordinates of the closest N points
% xyz1 = cat(1,nuc1,xyz);      
% %B = boundary(xyz1);
% [K,Vol] = convhulln(xyz1);    % each row of K is a triangle defined in terms of the point indices.      % step 5,6
% 
% % in(j) = inhull(nuc1,xyz1);                                                                             % alternative step 5,6
% %  if in(j) == 1                       % get the hull coordinates
% %vert{k} = K;
% %  end
% figure,plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.','MarkerSize',10), grid on
% hold on
% trisurf(K,xyz1(:,1),xyz1(:,2),xyz1(:,3),'Facecolor','y','FaceAlpha',0.3)
% zlim([1.5 4]);
% ylim([400 1100]);
% xlim([400 800]);
% end

% %%
% figure,trisurf(K,xyz1(:,1),xyz1(:,2),xyz1(:,3));
% hold on
% trisurf(K,xyz1(:,1),xyz1(:,2),xyz1(:,3),'Facecolor','red','FaceAlpha',0.1)
% %%
% figure,plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.','MarkerSize',10), grid on
% hold on
% trisurf(B,xyz1(:,1),xyz1(:,2),xyz1(:,3),'Facecolor','m','FaceAlpha',0.1)
% figure,plot3(xyz_fin(:,1),xyz_fin(:,2),xyz_fin(:,3),'r*','MarkerSize',10)
% 
% figure,plot3(V(2:end,1),V(2:end,2),V(2:end,3),'k*','MarkerSize',15), title('voronoi vertices coordinates')

