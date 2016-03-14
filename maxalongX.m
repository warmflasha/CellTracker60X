function [toelim2] = maxalongX(extrafilt,data_ch,data_mc,data_c2)

extra2= imerode(extrafilt,strel('disk',1));                                   
extra2bndr = extrafilt&~extra2;
stats_extra = regionprops(extra2bndr,'PixelList','PixelIdxList');             
%  add the loop here over the found objects (current code will only split
%  one 2-nuc object,not more merges

% if the filtering did not work, return the same mask
% each object here will have its own boundary  
data_extra1 = stats_extra(1).PixelList;     
data_extra2 = stats_extra(2).PixelList;

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

A = hull_bnd_only1;
B = obj_bnd_only1;
C = hull_bnd_only2;
D = obj_bnd_only2;


h = size(hull_bnd_only1,1);
maxdistX = zeros(h,1);
objrow = zeros(h,1);
% calculate the distance from the hull boundary to the cell object boundary along the x direction
for k=1:h;               
x1 = hull_bnd_only1(k,1);
y1 = hull_bnd_only1(k,2);
curr = find(obj_bnd_only1(:,2)== y1);%%%
if ~isempty(curr)
x2 = obj_bnd_only1(curr(end),1);
y2 = obj_bnd_only1(curr(end),2);
maxdistX(k,1) = sqrt(power((x1-x2),2)+power((y1-y2),2));
objrow(k,1) = curr(end);   % save the pixel coordinate of the point on the cell object boundary separately, since it is not the same as k-value for the hull boundary pixel coordinate
end
end
 mX = max(maxdistX);
 % r has the row numbers of the pixel coordinate (on the hull boundary ) of the point with max distance to the cell object
 % objrow  has the row numbers of the pixel coordinate (on the cell object boundary) of the point with max distance to the hull boundary
 [rX,~] = find(maxdistX == mX); 
                
% calculate the distance from the hull boundary of the second object to the cell object boundary along the x direction        
h2 = size(hull_bnd_only2,1);
maxdistX2 = zeros(h2,1);
objrow2 = zeros(h2,1);
for k=1:h2                    
x1 = hull_bnd_only2(k,1);
y1 = hull_bnd_only2(k,2);
curr2 = find(obj_bnd_only2(:,2)== y1);%%%%
if ~isempty(curr2)
x2 = obj_bnd_only2(curr2(end),1);
y2 = obj_bnd_only2(curr2(end),2);
%disp([num2str(y2) ]);
maxdistX2(k) = sqrt(power((x1-x2),2)+power((y1-y2),2));
objrow2(k,1) = curr2(end);
end
end
 mX2 = max(maxdistX2);             
 [rX2,~] = find(maxdistX2 == mX2);           % same as for the first object
 
  % COORDINATES OF CLOSEST POINTS,IF THE DISTANCES ARE CALCULATED ALONG THE X DIRECTION
% pt2, pt22 represent the pixel coordinates on the merged cell objects that need to be connected and where the cut
% should be made
% if isempty(objrow2(rX2(1))) || isempty(objrow(rX(1)))
%    toelim2 = [];
%    %return
% end 
if (objrow2(rX2(1)))~= 0 && objrow(rX(1))~=0
    ptXb = hull_bnd_only1(rX(1),:);             % since some max distances are repeated , will take the first one
    ptX = obj_bnd_only1(objrow(rX(1)),:);
    ptX2b = hull_bnd_only2(rX2(1),:);           % since some max distances are repeated , will take the first one
    ptX2 = obj_bnd_only2(objrow2(rX2(1)),:);

    
%  hold on                                         
% plot(ptXb(:,1),ptXb(:,2),'m*','markersize',30);
% plot(ptX(:,1),ptX(:,2),'b*','markersize',30);
% plot(ptX2b(:,1),ptX2b(:,2),'m*','markersize',30);
% plot(ptX2(:,1),ptX2(:,2),'b*','markersize',30);

x1 = ptX(1);   
y1 = ptX(2);   
x2 = ptX2(1);
y2 = ptX2(2);
a = (y1-y2)/(x1-x2);                  % draw a line through those points
b = 0.5*((y1+y2) -a*(x1+x2));

if x1<x2
vect1 = x1:0.0001:x2;
else
    vect1 = x2:0.0001:x1;
end
    
y_X=round(a*vect1+b);
toelim = cat(2,vect1',y_X');                
%plot(toelim(:,1),toelim(:,2),'c*');
toelim2 = intersect(data_mc,toelim,'rows');         % find where the line intersects with the cell object
% else
%     MaskFin2 = mask3;
%     return

else
   toelim2 = [];
   %return
end 
    
end