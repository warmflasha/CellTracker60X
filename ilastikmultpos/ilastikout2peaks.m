function ilastikout2peaks(segfilepath,samplepath)
% making ilastik output compatible with rest of the output
%samplepath = '/Users/sapnac18/Desktop/train_ilastik/multtif/';
dirinfo = dir(segfilepath);
startpos = postart(segfilepath);

sampledir = dir(samplepath);
samplepos = postart(samplepath);

mkdir(segfilepath, 'npeaks');

sampleno = 1;

%%
for file = startpos:size(dirinfo,1)
%% preprocessing
immask = h5read(strcat(segfilepath,'/',dirinfo(file).name), '/exported_data');
nchannel = size(immask,1);
ntime = size(immask,4);
%%
%Case 1: nchannel = 1

for time = 1:ntime
    mask = immask(1,:,:,time);
    mask_s = squeeze(mask);

    mask_b = im2bw(mask_s, 1); % cells:1, background:2
    mask_bc = imcomplement(mask_b);% we need background as '0'
    mask_bct = mask_bc'; % ilastik default settings have axes exchanged somehow.

%%
% peaks array
props = regionprops(mask_bct, 'Centroid', 'Area', 'PixelList');

nobj = size(props,1);

for i = 1:nobj
    peaks1(i,1:2) = props(i).Centroid;
    peaks1(i,3) = props(i).Area;
    objpxl{i} = props(i).PixelList;
    peaks1(i,4) = -1; 
    peaks1(:,5) = 1;
end
%%
% adding nuc intensity values to column 6

samp_im = strcat(samplepath, '/', sampledir(samplepos).name);

imch1 = imread(samp_im, time);

for i = 1:nobj
    pxl1 = objpxl{i};
    clear nucpxl;
    nucpxl = 0;
    for j = 1:size(pxl1,1)
     nucpxl = imch1(pxl1(j))+ nucpxl;
    end
    peaks1(i,6) = nucpxl;
end


%%
%selecting cells. 
% removing small spots- based on area? Lets trydetermining area threshold
thresh = 37;
peaks2 = peaks1(peaks1(:,3)>thresh,:);
peaks{time} = peaks2;


end    

fname = sprintf('/npeaks/fileseg%01d.mat', sampleno);
fname = strcat(segfilepath, fname);
save(fname, 'peaks');

sampleno = sampleno + 1;
samplepos = samplepos + 1;

end
