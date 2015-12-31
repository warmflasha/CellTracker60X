function ilastikout2peaks(segfilepath,samplepath,segchannel)

% making ilastik output compatible with rest of the output
% segfilepath: path corresponding to ilastik segmentation files
% samplepath: path corresponding to sample raw images
% segchannel: for a multi channel image, specify which channel was used for segmentation  


segdirinfo = dir(segfilepath);
startpos = postart(segfilepath);% getting the start position of the files.

sampledir = dir(samplepath);
samplepos = postart(samplepath); 

sampleno = 1;

%% deciphering file details :- no. of channels, no. of time points etc. from metadata)

samplefile = strcat(samplepath, '/', sampledir(samplepos).name);
reader = bfGetReader(samplefile);
nch = getSizeC(reader);
ntime = getSizeT(reader);
nz = getSizeZ(reader);

%%
mkdir(segfilepath, 'npeaks');

%%
for file = startpos:size(segdirinfo,1)
 
    immask = h5read(strcat(segfilepath ,segdirinfo(file).name), '/exported_data');
    
    %%
    clear peaks;
    %Case 1: nchannel = 1
   
    for time = 1:ntime
        mask = immask(segchannel,:,:,time);
        mask_s = squeeze(mask);
        mask_b = im2bw(mask_s,1); % object1:1, object2:2
        mask_b = imcomplement(mask_b);% if object 1 refers to background, comment this statement.
        mask_t = mask_b'; % ilastik default settings have axes transposed
        
        mask_final = bwareaopen(mask_t,100);
        
        %%
        % peaks array
        props = regionprops(mask_final, 'Centroid', 'Area', 'PixelList');
        nobj = size(props,1);
        
        for i = 1:nobj
            peaks1(i,1:2) = props(i).Centroid;
            peaks1(i,3) = props(i).Area;
            %objpxl{i} = props(i).PixelList;
            peaks1(i,4) = -1;
        end
        %%
        % adding nuc intensity values to column 5, channel1 intensity to
        % column6 and so on.
        
        samp_im = strcat(samplepath, sampledir(samplepos).name);
        
        for channel = 1:nch
            sampleimage = imread(samp_im, nch*(time-1)+channel);
            for i = 1:nobj
                pixellist = props(i).PixelList;
                
                ch_intensity = 0;
                
                for j = 1:size(pixellist,1)
                    ch_intensity(j) = sampleimage(pixellist(j,1), pixellist(j,2));
                end
                
                peaks1(i,4+channel) = sum(ch_intensity);
            end
                
        end
         
        peaks{time} = peaks1;
    end
    
    fname = sprintf('/npeaks/fileseg%01d.mat', sampleno);
    fname = strcat(segfilepath, fname);
    save(fname, 'peaks');
    
    sampleno = sampleno + 1;
    samplepos = samplepos + 1;
    
end
