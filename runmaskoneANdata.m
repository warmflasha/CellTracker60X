function peaks = runmaskoneANdata(ilastikdir,ilastikdircyto, imagedir,pos,tpt, timegroup,chan,paramfile)
%%
% ilastikdir: directory path of ilastik 2d segmentation probability density maps
% imagedir: directory path of the nuclear channel raw images 

ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/masks_zcyto');
ilastikdir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/Masks_z2');
imagedir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/rawimages_');
pos = 23;
tpt =40;
timegroup = 1;
chan = [1 2];
%%

[pnuc, inuc] = readmaskfiles1(ilastikdir,imagedir, pos,tpt, timegroup,chan(1));%

%%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
%%
% 
setUserParam3DsegmentationAN;
global userParam;

pmasks = primaryfilter(pnuc,userParam.logfilter, userParam.bthreshfilter, userParam.diskfilter, userParam.area1filter);

%%
[zrange, smasks] = secondaryfilter(pmasks, userParam.minstartobj, userParam.minsolidity, userParam.diskfilter, userParam.area2filter);
%%

[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange,CC] = traceobjectsz(smasks, userParam.matchdistance, zrange, userParam.zmatch);

%%
% [maskz] = 3Dlblmask(CC,nucleilist);

%%

[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userParam.overlapthresh);

%%

% plot the distinct nuclei
lb = labelmatrix(masterCCn);
a = label2rgb(lb);
figure, imshow(a);hold on
bw = regionprops(lb,'Centroid','PixelIdxList','Area');
bw1 = regionprops(masterCC,'Centroid','PixelIdxList','Area');
badinds = [bw.Area]< userParam.area2filter ;
bw(badinds) = [];

PILsSourcen;   % final planes where all the nuclei are 

Inew = zeros(1024,1024);
for j=1:size(bw,1)
Inew(bw(j).PixelIdxList) = j;  % make labeled objects from the final PixelIdxList
smasknew(:,:,j) = Inew;
end
figure, imshow(Inew,[]);

bw1 = regionprops(Inew,'Area','Centroid','PixelIdxList');
xy = [bw1.Centroid];
x = xy(1:2:end);
y =  xy(2:2:end);
hold on; plot(x,y,'*');%color(j,:)




