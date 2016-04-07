function [datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks_3Dsegm(mask1,mask2,img_nuc,img_cyto,paramfile,zrange)
% [datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks(mask1,mask2,img_nuc,img_cyto,paramfile)
% --------------------------------------------------------------
% takes masks produced by ilastik for two markers, and quantifies how much
% of marker 2 is overlapping with marker 1 vs not. Use case is if marker 1
% is nuclear marker while marker 2 is smad image. each case of marker 2
% must have marker 1 inside. 
%
% Inputs:
%   -mask1 - mask for marker 1
%   -mask2 - mask for marker 2
%   -img_nuc - actual image for marker 1
%   -img_cyto - actual image for marker 2
%   -paramfile - paramter file
% Outputs:
%   -datacell - output segmentation data in the usual format, one row per
%   cell
%   -Lnuc (final nuclear marker mask)
%   -Lcyto (final cytoplasmic mask (exclusive with nuclear marker)

global userParam;

try 
    eval(paramfile);
catch
    error('Error evaluating paramfile');
end

%process nuclear mask

% the nuc masks come already filtered and preprocessed
Lnuc = mask1; % this is a 3D mask now

%cytoplasmic mask

LcytoIl = (mask2);
Lcytonondil = LcytoIl;
LcytoIl = imdilate(LcytoIl,strel('disk',userParam.dilate_cyto)); 

%if no nuclei, exit
if sum(sum(sum(Lnuc))) == 0
    datacell = [];
    Lcytofin =zeros(size(Lnuc));
    return;
end

%raw images
I2 = img_cyto;
Inuc = img_nuc;

%this removes cytoplasms that don't have any nucleus inside.
% need this for the 3D case as well

cc = bwconncomp(LcytoIl);
cnuc = bwconncomp(Lnuc);                        
st = regionprops(cc,'PixelIdxList');         
stnuc = regionprops(cnuc,'PixelIdxList','Area');  %
goodinds = zeros(length(st),1);

for i = 1:length(stnuc);
    x =stnuc(i).PixelIdxList;
    for k=1:length(st);
        y = st(k).PixelIdxList;
        in = intersect(x,y);
        if ~isempty(in) && size(in,1)>userParam.areacytolow % the cytoplasms are at least 200 pixels ( bc if some junk was founf in nuclear channel it will likely have very small sytos
            goodinds(k,1) = k;
            
        end
    end
    
end
goodindsfin = nonzeros(goodinds);%
goodstats = st(goodindsfin);
% here need to leave the PixelIds of the goodinds and then convert back to
% the labeled image to get the final mask of the cytoplasms (still nuclei
% are in)
%
onebiglist = cat(1,goodstats.PixelIdxList);
Inew = zeros(1024,1024,size(mask1,3));
Inew(onebiglist) = true;
LcytoIl(Inew==0) =0;   %  this is the good cyto mask with only the cytos WITH nucleus
%Lnuc(Inew==0) =0; %   zero anything that may be in the nuc mask , buy not in cyto
%
%Next is similar processing to the above: now to remove the nuclei tha
%don't have the cytoplasms in them
cc = bwconncomp(LcytoIl); % use good cytoplasms now to screen bad nuclei
cnuc = bwconncomp(Lnuc);                        
st = regionprops(cc,'PixelIdxList');         
stnuc = regionprops(cnuc,'PixelIdxList','Area');  %
goodindsnuc = zeros(length(stnuc),1);

for i = 1:length(st);
    x =st(i).PixelIdxList;
        for k=1:length(stnuc);
        y = stnuc(k).PixelIdxList;
                in = intersect(x,y);
        if ~isempty(in)  
           goodindsnuc(k,1) = k;
           
        end
    end
    
end
goodindsfin = nonzeros(goodindsnuc); %
goodstats = stnuc(goodindsfin);

% here need to leave the PixelIds of the goodindsNUC 
%
onebiglist = cat(1,goodstats.PixelIdxList);
Inewnuc = zeros(1024,1024,size(mask1,3));
Inewnuc(onebiglist) = true;
Lnuc(Inewnuc==0) =0;   % now Lnuc has to have the same number of elements as the cytoplasm channel; Lnuc is still labeled
% now subtract those nuclei from filled cytoplasms to get final good
% labeled masks
% erode nuclei a little since sometimes causes problems
t = imerode(Lnuc,strel('disk',userParam.erode_nuc));      

LcytoIl(t>0)=0;                                           % remove nuclei from the cytoplasms
% return back to the non-dilated cyto masks
% and non-eroded nuc mask

LcytoIl(Lcytonondil ==0)=0;
LcytoIl(Lnuc > 0)=0; 
Lcytofin = LcytoIl;

% at this point have the set of 3D masks (nuc and cyto); I2 is the GFP
% channel 3D stack icyto(1024,1024,zrange)

I2proc = imopen(I2,strel('disk',userParam.small_rad));         % remove small bright stuff
I2proc = smoothImage(I2proc,userParam.gaussRadius,userParam.gaussSigma); %smooth
I2proc = presubBackground_self(I2proc);
% NOTE, the nuclear channel is not pre processed...

%get the NUCLEAR mean intensity for each labeled object

statsnuc = regionprops(Lnuc,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
statsnucw0 = regionprops(Lnuc,Inuc,'Area','Centroid','PixelIdxList','MeanIntensity');% these are the stats for the actual nuclear image(rfp)

badinds = [statsnuc.Area] < 0;   % don't need to filter anything since the number of elements in nuc and cyto is already matched in the code above      
badinds2 = [statsnucw0.Area] < 0;
statsnuc(badinds) = [];
statsnucw0(badinds2) = [];


%get the cytoplasmic mean intensity for each labeled object

statscyto = regionprops(Lcytofin,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
badinds = [statscyto.Area] < 0;  % don't need to filter anything since the number of elements in nuc and cyto is already matched in the code above 
statscyto(badinds) = [];


% ncells = length(statsN);
if size(Lcytofin,3) ==1
    xyz = round([statsnucw0.Centroid]);
xx =  xyz(1:2:end)';
yy =  xyz(2:2:end)';
zz =  (zrange*ones(1,size(statsnucw0,1)))';
xyzall = cat(2,xx,yy,zz);
else 
xyz = round([statsnucw0.Centroid]);
xx =  xyz(1:3:end)';
yy =  xyz(2:3:end)';
zz =  xyz(3:3:end)';
xyzall = cat(2,xx,yy,zz);
end

nuc_avrw0  = [statsnucw0.MeanIntensity]';%
nuc_areaw0  = [statsnucw0.Area]';%
nuc_avrw1 = round([statsnuc.MeanIntensity]');
cyto_area  = [statscyto.Area]';
cyto_avrw1  = round([statscyto.MeanIntensity]');
placeholder = -round(ones(length(xyzall(:,1)),1));

if isempty(statscyto)
    cyto_area = zeros(length(nuc_avrw1),1);
    cyto_avrw1 = cyto_area;
end

%this is done for whne all the previous clean up failed and still there is
%a mismatch between the nuc and cyto number of elements , only then remove
%that datapoint
if size(cyto_area,1) < size(nuc_areaw0,1) ||  size(cyto_area,1) > size(nuc_areaw0,1)
    datacell = [];
    return;
end

datacell=[xyzall(:,1) xyzall(:,2) xyzall(:,3) placeholder nuc_avrw0 nuc_avrw1 cyto_avrw1];%cyto_area


end