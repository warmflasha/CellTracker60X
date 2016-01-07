function peaks = nucCytoIlastik2peaksLoop(ilastikDirec,imageDirec1,imageDirec2,zplane,pos,chan,paramfile,outfile)
% peaks = nucCytoIlastik2peaksLoop(ilastikDirec,imageDirec,zplane,pos,chan,paramfile,outfile)
% ----------------------------------------------------------------------------
% runs nucCytoIlastik2peaks in a loop on all images at a particular
% position in a directory. Each image can have multiple time points but
% must be z slices and channels must be saved in separate images. 
% Inputs:
%   -ilastikDirec - ilastik directory with masks, one mask file per image
%   file
%   -imageDirec - directory with images (assumes Andor format)
%   -zplane - zplane to use (in Andor numbers, 0 possible)
%   -pos - position to use (in Andor numbers, 0 possible)
%   -chan - array of channels (1st is nuclear, 2nd to quantify, 0 possible)
%   -paramfile - parameter file to use
%   -outfile - output will be saved here
% Outputs:
%   -peaks - output segmentation in the usual format. also saved in outfile


% {CytoMaskt2}_{0000}.h5
% need to change the find string in the name 
if pos <10
[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec,'Cyto',['Mask' num2str(pos)]);% all cyto masks for the frame 0 ( four time groups)'Cyto','{0002}'['Outfile_000' num2str(pos) '_t']
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikDirec,'Nuc',['Mask' num2str(pos)]);% all nuc masks for the frame 0 ( four time groups)
end
if pos >=10
[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec,'Cyto',['Mask' num2str(pos)]);% all cyto masks for the frame 0 ( four time groups)'Cyto','{0002}'['Outfile_000' num2str(pos) '_t']
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikDirec,'Nuc',['Mask' num2str(pos)]);% all nuc masks for the frame 0 ( four time groups)
end    
nTprev = 0;
for j = 1:length(ilastikCytoAll)
    
    %get the ilastik masks
    ilastikNuc = fullfile(ilastikDirec,ilastikNucAll(j).name);
    ilastikCyto= fullfile(ilastikDirec,ilastikCytoAll(j).name);
    
    nuc_mask_all = h5read(ilastikNuc,'/exported_data');
    %nuc_mask_all = squeeze(nuc_mask_all);% this line is for the data that
    %was exported from the non headles mode run
    
    nuc_mask_all = nuc_mask_all(2,:,:,:);
    nuc_mask_all = squeeze(nuc_mask_all);
    
    
    cyto_mask_all = h5read(ilastikCyto,'/exported_data');
    %cyto_mask_all = squeeze(cyto_mask_all);% same comment
    cyto_mask_all = cyto_mask_all(2,:,:,:);
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
        
%  nuc_mask_all = nuc_mask_all(1,:,:,k); % 
%  nuc_mask_all = squeeze(nuc_mask_all);
% the following two lines are done within the nucCytoIlastik2peaks, but
% need to adjust to the appropriate h5 structure (different in ilastik is
% ran from feadless

%  nuc_mask_all = nuc_mask_all<1;
%  Lnuc =  bwareafilt(nuc_mask_all',[userParam.areanuclow userParam.areanuchi]);
   % same processing for the cyto channel     
        
        [outdat, Lnuc] = nucCytoIlastik2peaks(nuc_mask_all(:,:,k),cyto_mask_all(:,:,k),nuc_img,nuc_cyto,paramfile);%
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).compressNucMask = compressBinaryImg(Lnuc);
    end
    nTprev = nTprev + nT;
    
   % save(outfile,'peaks','imgfiles');
    save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles');
        
end
