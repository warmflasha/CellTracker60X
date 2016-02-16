function peaks = runmaskoneANdata(ilastikdir,ilastikdircyto, imagedir,pos,tpt, timegroup,chan,paramfile)
%%
% ilastikdir: directory path of ilastik 2d segmentation probability density maps
% imagedir: directory path of the nuclear channel raw images 


[pnuc, inuc] = readmaskfiles1(ilastikdir,imagedir, pos,tpt, timegroup,chan(1));%

for ii=1:size(inuc,3)
    img_now = inuc(:,:,ii);
    if ii==1
        max_img=img_now;
    else
        max_img=max(img_now,max_img);
    end
end
max_inuc = max_img;

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

[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userParam.overlapthresh);


%%
[a,b]=3Dlblmask(CC,nucleilist);



%%



% plot the distinct nuclei
lb = labelmatrix(masterCC);
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

%%
%ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/masks_zcyto');
[pcyto, icyto] = readmaskfiles1(ilastikdircyto,imagedir, pos,tpt, timegroup,chan(2));%
for ii=1:size(icyto,3)
    img_now = icyto(:,:,ii);
    if ii==1
        max_img=img_now;
    else
        max_img=max(img_now,max_img);
    end
end
max_icyto = max_img;

%%
% run the peaks separately on those slices (nuc and cyto) and combine into
% one peaks
paramfile = 'setUserParamLiveImagingAN';
zdata = cell(1,size(zrange,1));
n = cell(1,size(zrange,1));
c = cell(1,size(zrange,1));
for s=1:size(zrange,2);
%mask1 = pnuc(:,:,goodslices(s))';
mask1 = Inew'; % this is the labeled object already thresholded and processed
mask2 = pcyto(:,:,zrange(s))';  %sum(pcyto,3)' these need to be processed as usual
img_nuc = max_inuc;%inuc(:,:,goodslices(s));%;
img_cyto = max_icyto;%icyto(:,:,goodslices(s));%



[datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks_3Dsegm(mask1,mask2,img_nuc,img_cyto,paramfile);
zdata{s} = datacell;
n{s} = Lnuc;
c{s} = Lcytofin;
end
for s=1:(size(goodslices,1)-1)
fulldata = cat(1,zdata{s},zdata{s+1});

end



