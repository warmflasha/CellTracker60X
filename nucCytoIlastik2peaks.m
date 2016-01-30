function [datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks(mask1,mask2,img_nuc,img_cyto,paramfile)
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
%Lnuc = mask1 <1;
Lnuc = imfill(mask1 > userParam.probthresh_nuc,'holes');% for probabilities exported

Lnuc =  bwareafilt(Lnuc',[userParam.areanuclow userParam.areanuchi]);
%Lnuc = imerode(Lnuc,strel('disk',1));% erode by 1 pixel if some masks had touching cells still

%cytoplasmic mask
%LcytoIl = mask2 < 1;  
LcytoIl = im2bw(mask2,userParam.probthresh_cyto);% for probabilities exported
LcytoIl = (LcytoIl');
Lcytonondil = LcytoIl;
LcytoIl = imdilate(LcytoIl,strel('disk',userParam.dilate_cyto)); %this should be made into a parameter

%if no nuclei, exit
if sum(sum(Lnuc)) == 0
    datacell = [];
    Lcytofin =zeros(size(Lnuc));
    return;
end

stats = regionprops(Lnuc,'Centroid','PixelIdxList');
xy = [stats.Centroid];
xx = xy(1:2:end);
yy=xy(2:2:end);

vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

%raw images
I2 = img_cyto;
Inuc = img_nuc;

%this removes cytoplasms that don't have any nucleus inside.
cc = bwconncomp(LcytoIl+Lnuc & ~vImg);%& ~vImg
cnuc = bwconncomp(Lnuc);
st = regionprops(cc,'PixelIdxList');
stnuc = regionprops(cnuc,'PixelIdxList');
goodinds = zeros(length(st),1);

for i = 1:length(stnuc)
    x =stnuc(i).PixelIdxList;
    for k=1:length(st)
        y = st(k).PixelIdxList;
        in = intersect(x,y);
        if ~isempty(in)
            goodinds(k,1) = k;
        end
    end
end
goodindsfin = nonzeros(goodinds);
goodstats = st(goodindsfin);

% here need to leave the PixelIds of the goodinds and then convert back to the binary image by using ind2sub subtract the
% Lnuc to get the final mask of the cytoplasms
%
onebiglist = cat(1,goodstats.PixelIdxList);
Inew = zeros(1024,1024);
Inew(onebiglist) = true;

% erode nuclei a little since sometimes causes problems
t = imerode(Lnuc,strel('disk',userParam.erode_nuc)); %this should also be a parameter

LcytoIl = (Inew & ~ t & ~vImg);    % cyto masks initially include both nuclei+cyto, so need to eliminate nuc, use voronoi to divide;
% return back to the non-dilated cyto masks
% and non-eroded nuc mask
LcytoIl = bwlabel(LcytoIl,8);
LcytoIl(Lcytonondil ==0)=0;
LcytoIl(Lnuc ==1)=0;
Lcytofin = LcytoIl;

% at this point should have an array of nuc and cyto masks(from ilastik
% and watershed respectively)

I2proc = imopen(I2,strel('disk',userParam.small_rad)); % remove small bright stuff
I2proc = smoothImage(I2proc,userParam.gaussRadius,userParam.gaussSigma); %smooth
I2proc = presubBackground_self(I2proc);
% NOTE, the nuclear channel is not pre processed...


%get the NUCLEAR mean intensity for each labeled object
cc_nuc = bwconncomp(Lnuc,8);
statsnuc = regionprops(cc_nuc,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
statsnucw0 = regionprops(cc_nuc,Inuc,'Area','Centroid','PixelIdxList','MeanIntensity');% these are the stats for the actual nuclear image(rfp)

badinds = [statsnuc.Area] < userParam.areanuclow; %this area should also be a parameter
badinds2 = [statsnucw0.Area] < userParam.areanuclow;
statsnuc(badinds) = [];
statsnucw0(badinds2) = [];


%get the cytoplasmic mean intensity for each labeled object
%cc_cyto = bwconncomp(Lcytofin);
statscyto = regionprops(Lcytofin,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
badinds = [statscyto.Area] < userParam.areacytolow; 
statscyto(badinds) = [];


% ncells = length(statsN);
xy = stats2xy(statsnucw0);
nuc_avrw0  = [statsnucw0.MeanIntensity]';%[statsN.NuclearAvr];
nuc_areaw0  = [statsnucw0.Area]';%[statsN.NuclearArea];
nuc_avrw1 = round([statsnuc.MeanIntensity]');
cyto_area  = [statscyto.Area]';
cyto_avrw1  = round([statscyto.MeanIntensity]');
placeholder = -round(ones(length(xy(:,1)),1));
%

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

datacell=[xy(:,1) xy(:,2) nuc_areaw0 placeholder nuc_avrw0 nuc_avrw1 cyto_avrw1];%cyto_area


end



