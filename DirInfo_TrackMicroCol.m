
% info on where the images are and whidh position / timepoint to run

ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/masks_zcyto');
ilastikdir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/Masks_z2');
imagedir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/');%rawimages_
pos = 23;
tpt =57;
timegroup = 1;
chan = [1 2];

paramfile = 'setUserParamLiveImagingAN';
paramfile2 = 'setUserParam3DsegmentationAN';
%%
%run all timepoints and save into one peaks
% main fnction:

%[peaks,Lnuc,Lcytofin] = runmaskoneANdata(ilastikdir,ilastikdircyto, imagedir,pos,tpt, timegroup,chan,paramfile)

[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec2,['newCytoMasks' num2str(pos) '_']);
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikDirec1,['newNucMasks' num2str(pos) '_']);%

timegroups = 1;

ilastikCytoAll = ilastikCytoAll((pos*timegroups+1):(pos+1)*timegroups); % 1-4; 5-8;...
ilastikNucAll = ilastikNucAll((pos*timegroups+1):(pos+1)*timegroups);


nTprev = 0;
for j = 1:length(ilastikCytoAll)
    
    %get the ilastik masks
    ilastikNuc = fullfile(ilastikDirec1,ilastikNucAll(j).name);
    ilastikCyto= fullfile(ilastikDirec2,ilastikCytoAll(j).name);
    
    nuc_mask_all = h5read(ilastikNuc,'/exported_data');
    %nuc_mask_all = squeeze(nuc_mask_all);
    
    nuc_mask_all = nuc_mask_all(2,:,:,:);% for probabilities exported
    nuc_mask_all = squeeze(nuc_mask_all);
    
    
    cyto_mask_all = h5read(ilastikCyto,'/exported_data');
    %cyto_mask_all = squeeze(cyto_mask_all);% same comment
    cyto_mask_all = cyto_mask_all(2,:,:,:);% for probabilities exported
    cyto_mask_all = squeeze(cyto_mask_all);
    
    
    %get the image readers
    ff1 = readAndorDirectory(imageDirec1);
    ff2 = readAndorDirectory(imageDirec2);% AN
    filename1 = getAndorFileName(ff1,pos,ff1.t(j),zplane,chan(1));
    filename2 = getAndorFileName(ff2,pos,ff2.t(j),zplane,chan(2));
    img_nuc_reader = bfGetReader(filename1);
    img_cyto_reader = bfGetReader(filename2);
    
    nT = img_nuc_reader.getSizeT;
    
    for k = 1:nT
        
        plane1 = img_nuc_reader.getIndex(0,0, k - 1) + 1;
        nuc_img = bfGetPlane(img_nuc_reader,plane1);
        
        plane1 = img_cyto_reader.getIndex(0,0, k - 1) + 1;
        nuc_cyto = bfGetPlane(img_cyto_reader,plane1);
        
               
       % [outdat, Lnuc,Lcytofin] = nucCytoIlastik2peaks(nuc_mask_all(:,:,k),cyto_mask_all(:,:,k),nuc_img,nuc_cyto,paramfile);%
        peaks{nTprev+k} = outdat;
        
        if sum(sum(Lnuc)) == 0
            imgfiles(nTprev+k).compressNucMask = [];
        else
            imgfiles(nTprev+k).compressNucMask = compressBinaryImg(Lnuc);
        end
        if sum(sum(Lcytofin)) == 0
            imgfilescyto(nTprev+k).compressNucMask = [];
            
        else
            imgfilescyto(nTprev+k).compressNucMask = compressBinaryImg(Lcytofin);
        end
        
    end
    nTprev = nTprev + nT;
end

% save(outfile,'peaks','imgfiles');
save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');

end
