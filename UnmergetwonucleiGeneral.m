function [MaskFin2,maskfin] = UnmergetwonucleiGeneral(mask3)% input the nuclear mask
% the mask to input should be already thresholded binary from probability mask of ilastik
% function checks if the nuclear mask has two cells merged
% if yes, then splits them and returns a new mask where all the two merged
% nuclei are separated
% checks every nucleis in the mask, not limited to one object;
% if there are no merged objects then returns the same mask as input
clear MaskFin2
clear maskfin
% 
% get the image with only the merged object
global userParam;
mask3new = bwareafilt(mask3,[600 20000]);  %  userParam.areanuclow          
stats = bwconncomp(mask3new);
% mask3new is the binary image with only the merged object

nn = (stats.NumObjects);
masktmp = cell(1, nn);

 for ii=1:nn
 masktmp{ii} = zeros(1024,1024);
 masktmp{ii}(stats.PixelIdxList{ii}) = 1;
 MaskFin2{ii} = masktmp{ii};
 end
 for ii=1:nn
 [a,b,extrafilt,data_ch,data_mc,data_c2] = testmergenuclei(masktmp{ii});
if a<30 && b<30
    MaskFin2{ii} = masktmp{ii};
else
    
% get the boundaries of the two extra objects
%figure, imshow(extrafilt);hold on
extra2= imerode(extrafilt,strel('disk',1));                                   
extra2bndr = extrafilt&~extra2;
stats_extra = regionprops(extra2bndr,'PixelList','PixelIdxList');             
%  add the loop here over the found objects (current code will only split
%  one 2-nuc object,not more merges

% if the filtering did not work, return the same mask
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
if isempty(hull_bnd_only1)||isempty(hull_bnd_only2)||isempty(obj_bnd_only1)||isempty(obj_bnd_only2)
  MaskFin2{ii} = masktmp{ii};
  continue
end


% hold on; plot(obj_bnd_only1(:,1),obj_bnd_only1(:,2),'b*','markersize',15);
% plot(hull_bnd_only1(:,1),hull_bnd_only1(:,2),'y*','markersize',15);
% plot(hull_bnd_only2(:,1),hull_bnd_only2(:,2),'m*','markersize',15);
% plot(obj_bnd_only2(:,1),obj_bnd_only2(:,2),'g*','markersize',15);

% find the max disatce between the hull boundary and the cell object
% boundary, for each of the two objects (one unique max for each object)
% the two maxima define the points where the two cells are merged and where
% the cut needs to be made

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
end
%----------------------------------------------------
%  COORDINATES OF CLOSEST POINTS,IF THE DISTANCES ARE CALCULATED ALONG THE Y DIRECTION
% calculate the distance from the hull boundary of the second object to the cell object boundary along the y direction 
h2 = size(hull_bnd_only2,1);  
maxdistY2 = zeros(h2,1);
objrowY2 = zeros(h2,1);
for k=1:h2                    
         x1 = hull_bnd_only2(k,1);
         y1 = hull_bnd_only2(k,2);
         curr2 = find(obj_bnd_only2(:,1)== x1);%%%
         if ~isempty(curr2)
             x2 = obj_bnd_only2(curr2(end),1);
             y2 = obj_bnd_only2(curr2(end),2);
             %disp([num2str(y2) ]);
             maxdistY2(k) = sqrt(power((x1-x2),2)+power((y1-y2),2));
             objrowY2(k,1) = curr2(end);
         end
     end
  mY2 = max(maxdistY2);             
 [rY2,~] = find(maxdistY2 == mY2); 
       % pixel coordinated of the line (line goes through the whole image)


%plot(toelim2(:,1),toelim2(:,2),'c*');

% calculate the distance from the first hull boundary to the cell object boundary along the y direction 
h = size(hull_bnd_only1,1);  
maxdistY = zeros(h,1);
objrowY = zeros(h,1);
                for k=1:h;               
                   x1 = hull_bnd_only1(k,1);
                   y1 = hull_bnd_only1(k,2);
                   curr = find(obj_bnd_only1(:,1)== x1);%%%%%
                   if ~isempty(curr)
                       x2 = obj_bnd_only1(curr(end),1);
                       y2 = obj_bnd_only1(curr(end),2);
                       maxdistY(k,1) = sqrt(power((x1-x2),2)+power((y1-y2),2));
                       objrowY(k,1) = curr(end);   
                    end
                   
               end
               mY = max(maxdistY);
               [rY,~] = find(maxdistY == mY);
               if objrowY(rY(1))~= 0 && objrowY2(rY2(1))~=0
                   ptYb = hull_bnd_only1(rY(1),:);
                   ptY = obj_bnd_only1(objrowY(rY(1)),:);
                   ptY2b = hull_bnd_only2(rY2(1),:);
                   ptY2 = obj_bnd_only2(objrowY2(rY2(1)),:);
               
xx1 = ptY(1);   
yy1 = ptY(2);   
xx2 = ptY2(1);
yy2 = ptY2(2);
a = (yy1-yy2)/(xx1-xx2);                  % draw a line through those points
b = 0.5*((yy1+yy2) -a*(xx1+xx2));
if xx1<xx2
vect2 = xx1:0.0001:xx2;
else
    vect2 = xx2:0.0001:xx1;
end
y_Y=round(a*vect2+b);
toelim_y = cat(2,vect2',y_Y');        
toelim2_y = intersect(data_mc,toelim_y,'rows');         % find where the line intersects with the cello object
else 
MaskFin2{ii} = masktmp{ii};
return
end
% Choose which line has the least intersecting points with the merged
% object = this is the line where the cut wll be made
if ((objrow2(rX2(1)))~= 0 && objrow(rX(1))~=0)==0
    toelimfin = toelim2_y;
    I = zeros(1024,1024);                                                    
    linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
    I(linearInd)=1;
    II = imdilate(I,strel('disk',4)); % 'disk',4            

    MaskFin2{ii} = masktmp{ii};  %&~II                                     
    %MaskFin2 = MaskFin + mask3old;                               
else
if size(toelim2_y,1)>size(toelim2,1)
toelimfin = toelim2;
else
toelimfin = toelim2_y;
end

I = zeros(1024,1024);                                    % create an image with only that element                 
linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
I(linearInd)=1;
II = imdilate(I,strel('disk',4)); % 'disk',4             % dilate a little in order to create a merged line 

MaskFin2{ii} = masktmp{ii}; %&~II                                      % remove those pixels from the merged object
%MaskFin2 = MaskFin + mask3old;                                % return to the original nuclear mask but with the cells unmerged
%figure, imshow(MaskFin2);
end
 end
 maskfin = zeros(1024,1024);
 for ii=1:size(MaskFin2,2)
     stats = regionprops(MaskFin2{ii},'PixelIdxList');
      maskfin(stats.PixelIdxList) = 1;
 end
end