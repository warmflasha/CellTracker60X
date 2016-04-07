function [toelim2_y] = maxalongY(extrafilt,data_ch,data_mc,data_c2)



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
                   toelim2_y = [];
                   
               end

end