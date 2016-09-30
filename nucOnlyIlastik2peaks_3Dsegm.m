function [datacell,Lnuc] = nucOnlyIlastik2peaks_3Dsegm(mask1,img_nuc,paramfile,imagedir,selectZ,pos,timegroup)
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

%if no nuclei, exit
if sum(sum(sum(Lnuc))) == 0
    datacell = [];
       return;
end

%raw images
%I2 = img_cyto;
Inuc = img_nuc;


%I2proc = imopen(I2,strel('disk',userParam.small_rad));         % remove small bright stuff
%I2proc = smoothImage(I2proc,userParam.gaussRadius,userParam.gaussSigma); %smooth 
%I2proc = presubBackground_self(I2proc);%I2proc           old 

                                                          % new method to subtract background
% now need to get the raw images in the non-nuclearmarker channel (actual
% staining)
ff = readAndorDirectoryANmod(imagedir);
allchanels=cell(1,(size(ff.w,2))) ;

for jj=1:(size(ff.w,2)) % start from 1 so to select the nuc image plane correctly (selectZ)
[imgsnuc_reader]   =  getrawimgfilesselectZ(imagedir,selectZ,(pos-1),timegroup,jj);  
k = 1;% first and only time point
 
    for m = 1:size(imgsnuc_reader,2) %
        planenuc = imgsnuc_reader{m}.getIndex(0,0, k - 1) + 1;
        inuc(:,:,m) = bfGetPlane(imgsnuc_reader{m},planenuc);
         
    end
    % here need the adjusment if these are more than one slices
      allchanels{jj} = inuc;
   
end

I2proc = zeros(1024,1024,size(ff.w,2));
for ll=1:size(ff.w,2)
I2proc(:,:,ll) = simplebg([],Lnuc,allchanels{ll}); % here I2proc contains the bg-subtracted images for the good zplane for all existing chanels
end

% get the stats form each channel using the Lnuc mask and put into the
% datcell matrix

statsnuc = regionprops(Lnuc,I2proc(:,:,1),'Area','Centroid','PixelIdxList','MeanIntensity');% stas for the nuclear marker

if size(Lnuc,3) ==1
    xyz = round([statsnuc.Centroid]);
xx =  xyz(1:2:end)';
xyzall = zeros(size(xx,1),3); % initialize the atrix for the data once the number of rows is known (ftom size of xx)
yy =  xyz(2:2:end)';
zz =  selectZ*ones(size(xx,1),1);
xyzall = cat(2,xx,yy,zz);
else 
xyz = round([statsnuc.Centroid]);
xx =  xyz(1:3:end)';
xyzall = zeros(size(xx,1),3); % initialize the matrix for the data once the number of rows is known (ftom size of xx)
yy =  xyz(2:3:end)';
zz =  xyz(3:3:end)';
xyzall = cat(2,xx,yy,zz);
end

nuc_avr  = [statsnuc.MeanIntensity]';%
nuc_area  = [statsnuc.Area]';%

datacell=[xyzall(:,1) xyzall(:,2) xyzall(:,3) nuc_area nuc_avr ];%

for k=2:size(ff.w,2) % (populate the data for the rest of the chanels (so start the loop from 2)
statsother= regionprops(Lnuc,I2proc(:,:,k),'Area','Centroid','PixelIdxList','MeanIntensity');
datacell(:,5+k-1) = [statsother.MeanIntensity]';
end



end