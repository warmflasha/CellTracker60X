function [MaskFin2] = Unmergetwonuclei(mask3)% input the nuclear mask
% the mask to input should be already thresholded from probability mask of ilastik
% function checks if the nuclear mask has two cells merged
% if yes, then splits them and returns a new mask where all the two merged
% nuclei are separated
% if there are no merged objects then returns the same mask as input

% 
% get the image with only the merged object
global userParam;
mask3new = bwareafilt(mask3,[userParam.areanuclow_unmerge 20000]);               % see if there are any merged nuclei
mask3old = bwareafilt(mask3,[userParam.areanuclow userParam.areanuclow_unmerge]);% create the mask without the merged object for later use
stats = bwconncomp(mask3new);
% check how many merged objects there are
if (sum(sum(mask3new)) == 0) || (stats.NumObjects >1); % if there are no merged objects, return the same nuc mask as was input
MaskFin2 = mask3;
return
end
% if there is only one cell in the image, don't split it
mask3test = bwareafilt(mask3,[userParam.areanuclow userParam.areanuchi]); 
stats2 = bwconncomp(mask3test);
if (stats2.NumObjects == 1)
    MaskFin2 = mask3;
return
end
% mask3new is the binart image with only merged object

%bw = bwlabel(mask3new);
bw2 = bwconncomp(mask3new);
stats3 = regionprops(bw2,'PixelList','ConvexHull','Area','ConvexArea','ConvexImage');
ch = stats3.ConvexHull;
chi = bwconvhull(mask3new);               % image that has the convex hull object filled
data_mc = stats3.PixelList;               % pixels of the merged cell, needed for later

% get the boundary of convex hull object
mask_ch = imerode(chi,strel('disk',1));
boundary_ch = chi&~mask_ch;
stboundary_ch = regionprops(boundary_ch,'PixelList');
%hold on
%plot(stboundary_ch.PixelList(:,1),stboundary_ch.PixelList(:,2),'c*','markersize',15);
data_ch = stboundary_ch.PixelList;       % pixel coordinates of the convex hull boundary

%hold on
%plot(stats3.Centroid(:,1),stats3.Centroid(:,2),'b*','markersize',15)
% get the boundary of the merged cell object   
mask4_2 = imdilate(mask3new,strel('disk',1));  
boundary2 = mask4_2&~mask3new;
statsboundary2 = regionprops(boundary2,'PixelList');
%hold on
%plot(statsboundary2.PixelList(:,1),statsboundary2.PixelList(:,2),'m*','markersize',15);

data_c2 = statsboundary2.PixelList;         % pixel coordinates of the boundary of the cell merged object
extra = chi&~mask3new;                      % extra stuff = difference between the filled convex hull area and the merged cell object

extrafilt = bwareafilt(extra,[10 userParam.areanuclow]);    % filter out very small intersections
stextra = bwconncomp(extrafilt);
stextra2 = regionprops(extrafilt,'Area','PixelIdxList'); 

% if after (object-filled convex hull) there are more than two ojects left,
% leave the largest two (to be used for the splitting)
if stextra.NumObjects > 2
a = [stextra2.Area];
a2 = sort(a,'descend');
edg1 = find(a == a2(1));
edg2 = find(a == a2(2));
extrafilt = zeros(1024,1024);
extrafilt(stextra2(edg1).PixelIdxList) = 1;
extrafilt(stextra2(edg2).PixelIdxList) = 1;
end
% if the shape was so weird, that only one object is left in extra, then return
if stextra.NumObjects == 1
MaskFin2 = mask3;
return
end
% get the boundaries of the two extra objects
%figure, imshow(extrafilt);hold on
extra2= imerode(extrafilt,strel('disk',1));                                   
extra2bndr = extrafilt&~extra2;
stats_extra = regionprops(extra2bndr,'PixelList','PixelIdxList');             
%  add the loop here over the found objects (current code will only split
%  one 2-nuc object,not more merges

% if the filtering did not work, return the same mask
if size(stats_extra,1)>2 || size(stats_extra,1)< 2
    MaskFin2 = mask3;
    return
else
% each object here will have its own boundary  
data_extra1 = stats_extra(1).PixelList;     
data_extra2 = stats_extra(2).PixelList;
end
% plot the extra stuff
% hold on
% plot(stats_extra(1).PixelList(:,1),stats_extra(1).PixelList(:,2),'r*','markersize',15);
% plot(stats_extra(2).PixelList(:,1),stats_extra(2).PixelList(:,2),'r*','markersize',15);
 
% separate the boundary of each object into the one that is the boundary
% with the cell or the one with the convex hull ( for two objects will get
% four separate boundaries
hull_bnd_only1 = intersect(data_ch,data_extra1,'rows');        
hull_bnd_only2 = intersect(data_ch,data_extra2,'rows');

obj_bnd_only1 = intersect(data_c2,data_extra1,'rows');
obj_bnd_only2 = intersect(data_c2,data_extra2,'rows');
% remove the points where those boundaries intersect 
[~,r1,c1] = intersect(hull_bnd_only1,obj_bnd_only1,'rows');
[~,r2,c2] = intersect(hull_bnd_only2,obj_bnd_only2,'rows');
hull_bnd_only1(r1,:)=[];
obj_bnd_only1(c1,:)=[];
hull_bnd_only2(r2,:)=[];
obj_bnd_only2(c2,:)=[];

% hold on; plot(obj_bnd_only1(:,1),obj_bnd_only1(:,2),'b*','markersize',15);
% plot(hull_bnd_only1(:,1),hull_bnd_only1(:,2),'y*','markersize',15);
% plot(hull_bnd_only2(:,1),hull_bnd_only2(:,2),'m*','markersize',15);
% plot(obj_bnd_only2(:,1),obj_bnd_only2(:,2),'g*','markersize',15);

% find the max disatce between the hull boundary and the cell object
% boundary, for each of the two objects (one unique max for each object)
% the two maxima define the points where the two cells are merged and where
% the cut needs to be made

h = size(hull_bnd_only1,1);
maxdist1 = zeros(h,1);
objrow = zeros(h,1);
% calculate the distance from the hull boundary to the cell object boundary along the x direction
for k=1:h;               
x1 = hull_bnd_only1(k,1);
y1 = hull_bnd_only1(k,2);
curr = find(obj_bnd_only1(:,2)== y1);%%%
if ~isempty(curr)
x2 = obj_bnd_only1(curr(end),1);
y2 = obj_bnd_only1(curr(end),2);
maxdist1(k,1) = sqrt(power((x1-x2),2)+power((y1-y2),2));
objrow(k,1) = curr(end);   % save the pixel coordinate of the point on the cell object boundary separately, since it is not the same as k-value for the hull boundary pixel coordinate
end
end
 m = max(maxdist1);
 % r has the row numbers of the pixel coordinate (on the hull boundary ) of the point with max distance to the cell object
 % objrow  has the row numbers of the pixel coordinate (on the cell object boundary) of the point with max distance to the hull boundary
 [r,~] = find(maxdist1 == m); 
 % calculate the distance from the hull boundary to the cell object boundary along the y direction 
                for k=1:h;               
                   x1 = hull_bnd_only1(k,1);
                   y1 = hull_bnd_only1(k,2);
                   curr = find(obj_bnd_only1(:,1)== x1);%%%%%
                   if ~isempty(curr)
                       x2 = obj_bnd_only1(curr(end),1);
                       y2 = obj_bnd_only1(curr(end),2);
                       maxdist1(k,1) = sqrt(power((x1-x2),2)+power((y1-y2),2));
                       objrow(k,1) = curr(end);   
                    end
                   
               end
               m12 = max(maxdist1);
               [r2,~] = find(maxdist1 == m12);
               
% calculate the distance from the hull boundary of the second object to the cell object boundary along the x direction        
h2 = size(hull_bnd_only2,1);
maxdist2 = zeros(h2,1);
objrow2 = zeros(h2,1);
for k=1:h2                    
x1 = hull_bnd_only2(k,1);
y1 = hull_bnd_only2(k,2);
curr2 = find(obj_bnd_only2(:,2)== y1);%%%%
if ~isempty(curr2)
x2 = obj_bnd_only2(curr2(end),1);
y2 = obj_bnd_only2(curr2(end),2);
%disp([num2str(y2) ]);
maxdist2(k) = sqrt(power((x1-x2),2)+power((y1-y2),2));
objrow2(k,1) = curr2(end);
end
end

 m2 = max(maxdist2);             
 [rr,~] = find(maxdist2 == m2);           % same as for the first object
 
% calculate the distance from the hull boundary of the second object to the cell object boundary along the y direction (y works too)
     for k=1:h2                    
         x1 = hull_bnd_only2(k,1);
         y1 = hull_bnd_only2(k,2);
         curr2 = find(obj_bnd_only2(:,1)== x1);%%%
         if ~isempty(curr2)
             x2 = obj_bnd_only2(curr2(end),1);
             y2 = obj_bnd_only2(curr2(end),2);
             %disp([num2str(y2) ]);
             maxdist2(k) = sqrt(power((x1-x2),2)+power((y1-y2),2));
             objrow2(k,1) = curr2(end);
         end
     end
  m22 = max(maxdist2);             
 [rr2,~] = find(maxdist2 == m22); 
 % COORDINATES OF CLOSEST POINTS,IF THE DISTANCES ARE CALCULATED ALONG THE X DIRECTION
% pt2, pt22 represent the pixel coordinates on the merged cell objects that need to be connected and where the cut
% should be made
pt1 = hull_bnd_only1(r(1),:);             % since some max distances are repeated , will take the first one
pt2 = obj_bnd_only1(objrow(r(1)),:); 
pt12 = hull_bnd_only2(rr(1),:);           % since some max distances are repeated , will take the first one
pt22 = obj_bnd_only2(objrow2(rr(1)),:);
%  hold on                                         
% plot(pt1(:,1),pt1(:,2),'m*','markersize',30);
% plot(pt2(:,1),pt2(:,2),'m*','markersize',30);
% plot(pt12(:,1),pt12(:,2),'m*','markersize',30);
% plot(pt22(:,1),pt22(:,2),'m*','markersize',30);

x1 = pt2(1);   
y1 = pt2(2);   
x2 = pt22(1);
y2 = pt22(2);
a = (y1-y2)/(x1-x2);                  % draw a line through those points
b = 0.5*((y1+y2) -a*(x1+x2));

if x1<x2
vect1 = x1:0.001:x2;
else
    vect1 = x2:0.001:x1;
end
    
ynew=round(a*vect1+b);
toelim = cat(2,vect1',ynew');                       % pixel coordinated of the line (line goes through the whole image)

%plot(toelim(:,1),toelim(:,2),'c*');
toelim2 = intersect(data_mc,toelim,'rows');         % find where the line intersects with the cell object
%plot(toelim2(:,1),toelim2(:,2),'c*');
%  COORDINATES OF CLOSEST POINTS,IF THE DISTANCES ARE CALCULATED ALONG THE Y DIRECTION
pt1 = hull_bnd_only1(r2(1),:);             
pt2 = obj_bnd_only1(objrow(r2(1)),:); 
pt12 = hull_bnd_only2(rr2(1),:);           
pt22 = obj_bnd_only2(objrow2(rr2(1)),:);

x1 = pt2(1);   
y1 = pt2(2);   
x2 = pt22(1);
y2 = pt22(2);
a = (y1-y2)/(x1-x2);                  % draw a line through those points
b = 0.5*((y1+y2) -a*(x1+x2));
if x1<x2
vect1 = x1:0.001:x2;
else
    vect1 = x2:0.001:x1;
end
ynew=round(a*vect1+b);
toelim_y = cat(2,vect1',ynew');        
toelim2_y = intersect(data_mc,toelim_y,'rows');         % find where the line intersects with the cello object
% Choose which line has the least intersecting points with the merged
% object = this is the line where the cut wll be made
if size(toelim2_y,2)>size(toelim,2)
toelimfin = toelim;
else
toelimfin = toelim2_y;
end

I = zeros(1024,1024);                                    % create an image with only that element                 
linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
I(linearInd)=1;
II = imdilate(I,strel('disk',1)); % 'disk',4             % dilate a little in order to create a merged line 
%vect1 = toelimfin(1,1):0.1:toelimfin(end,1);
%toelimfin2=a*vect1+b; 
%toelimfin2 = cat(2,vect1',toelimfin2');
%linearInd = sub2ind(size(I), toelimfin2(:,2), toelimfin2(:,1));
%I(linearInd)=[];


MaskFin = mask3new&~II;                                       % remove those pixels from the merged object
MaskFin2 = MaskFin + mask3old;                                % return to the original nuclear mask but with the cells unmerged
%figure, imshow(MaskFin2);
end