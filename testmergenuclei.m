function [area1,area2,rx,ry,ch_x,ch_y,extrafilt,data_ch,data_mc,data_c2] = testmergenuclei(mask)

global userParam 
mask3new = mask;
bw2 = bwconncomp(mask3new);
stats3 = regionprops(bw2,'PixelList','ConvexHull','Area','ConvexArea','ConvexImage');
ch = stats3.ConvexHull;
chi = bwconvhull(mask3new);               % image that has the convex hull object filled
data_mc = stats3.PixelList;               % pixels of the merged cell, needed for later
% if the cell isnvery small don't even check if it is merged
if  stats3.Area < userParam.areanuclow_unmerge
extrafilt = mask;%%%%%%%%%%%%%%%%%%%%%
area1 = 0;
area2 = 0;
rx=[];
ry=[];
ch_x=[];
ch_y=[];
data_ch = [];
data_mc = [];
data_c2 = [];

return
end
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

extrafilt = extra;   
stextra = bwconncomp(extrafilt);
stextra2 = regionprops(extrafilt,'Area','PixelIdxList'); 


% if the shape was so weird, that only one object is left in extra, then return
if stextra.NumObjects == 1
extrafilt = mask;%%%%%%%%%%%%%%%%%%%%%
area1 = 0;
area2 = 0;
rx=[];
ry=[];
ch_x=[];
ch_y=[];
data_ch = [];
data_mc = [];
data_c2 = [];

return
end

% if after (object-filled convex hull) there are more than two ojects left,
% leave the largest two (to be used for the splitting)

if stextra.NumObjects > 2
    a = [stextra2.Area];
    a1 = unique(a);
    a2 = sort(a1,'descend');
end
if size(a2,2) == 1
extrafilt = mask;%%%%%%%%%%%%%%%%%%%%%
area1 = 0;
area2 = 0;
rx=[];
ry=[];
ch_x=[];
ch_y=[];
data_ch = [];
data_mc = [];
data_c2 = [];
    return
end
if size(a2,2) > 1
    edg1 = find(a == a2(1));
    edg2 = find(a == a2(2));
    r = regionprops(im2bw(mask),'Centroid');
    rx = round(r.Centroid(1));
    ry = round(r.Centroid(2));
    rch = regionprops(chi,'Centroid');
    ch_x = round(r.Centroid(1));
    ch_y = round(r.Centroid(2));
    
    
    inew = zeros(1024,1024);
    inew(stextra2(edg1(1)).PixelIdxList) = 1;
    inew(stextra2(edg2(1)).PixelIdxList) = 1;
end

extrafilt = inew;
area1 = stextra2(edg1(1)).Area;
area2 = stextra2(edg2(1)).Area;
end