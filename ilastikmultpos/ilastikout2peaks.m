function ilastikout2peaks(segfilepath,samplepath,segchannel)
% making ilastik output compatible with rest of the output
%samplepath = '/Users/sapnac18/Desktop/train_ilastik/multtif/';
dirinfo = dir(segfilepath);
startpos = postart(segfilepath);

sampledir = dir(samplepath);
samplepos = postart(samplepath);

sampleno = 1;

%% deciphering file details(no. of channels, no. of time points etc. from metadata)
samplefl = strcat(samplepath, '/', sampledir(samplepos).name);
reader = bfGetReader(samplefl);
nch = getSizeC(reader);
ntime = getSizeT(reader);
nz = getSizeZ(reader);
%%

mkdir(segfilepath, 'npeaks');

%%
for file = startpos:size(dirinfo,1)
%% preprocessing
immask = h5read(strcat(segfilepath ,dirinfo(file).name), '/exported_data');
%nchannel = size(immask,1);
%ntime = size(immask,4);
%%
%Case 1: nchannel = 1

ntime = 2;
for time = 1:ntime
    mask = immask(segchannel,:,:,time);
    mask_s = squeeze(mask);
    thresh = graythresh(mask_s); 
    mask_b = im2bw(mask_s, thresh); % cells:1, background:2
    %mask_bc = imcomplement(mask_b);% we need background as '0'
    mask_bct = mask_b'; % ilastik default settings have axes exchanged somehow.

%%
% peaks array
props = regionprops(mask_bct, 'Centroid', 'Area', 'PixelList');

nobj = size(props,1);

for i = 1:nobj
    peaks1(i,1:2) = props(i).Centroid;
    peaks1(i,3) = props(i).Area;
    objpxl{i} = props(i).PixelList;
    peaks1(i,4) = -1; 
    %peaks1(:,5) = 1;
end
%%
% adding nuc intensity values to column 5

samp_im = strcat(samplepath, sampledir(samplepos).name);

for ch_no = 1:nch
 imch{ch_no} = imread(samp_im, nch*(time-1)+ch_no); 
 chint{ch_no} = 0;
end

%%
for i = 1:nobj
    pxl1 = objpxl{i};
   
    chint(1,:) = {0};
    for j = 1:size(pxl1,1)
        for ch_no = 1:nch
         chint{ch_no} = imch{ch_no}(pxl1(j))+ chint{ch_no};
        end
    end
    for ch_no = 1:nch
      peaks1(i,4+ch_no) = chint{ch_no};
    end
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
