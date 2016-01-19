function peaks = nucCytoIlastik2peaksLoop(ilastikDirec1,ilastikDirec2,imageDirec1,imageDirec2,zplane,pos,chan,paramfile,outfile)
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


[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec2,'CytoMask');
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikDirec1,'NucMask');% 

timegroups = 3;

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
        
%  nuc_mask_all = nuc_mask_all(1,:,:,k); % 
%  nuc_mask_all = squeeze(nuc_mask_all);
% the following two lines are done within the nucCytoIlastik2peaks, but
% need to adjust to the appropriate h5 structure (different in ilastik is
% ran from feadless

%  nuc_mask_all = nuc_mask_all<1;
%  Lnuc =  bwareafilt(nuc_mask_all',[userParam.areanuclow userParam.areanuchi]);
   % same processing for the cyto channel     
        
        [outdat, Lnuc,Lcytofin] = nucCytoIlastik2peaks(nuc_mask_all(:,:,k),cyto_mask_all(:,:,k),nuc_img,nuc_cyto,paramfile);%
        peaks{nTprev+k} = outdat;
        if sum(sum(Lnuc)) == 0 || sum(sum(Lcytofin)==0)
         imgfiles(nTprev+k).compressNucMask = [];
         imgfilescyto(nTprev+k).compressNucMask = [];
        else
          imgfiles(nTprev+k).compressNucMask = compressBinaryImg(Lnuc);
          imgfilescyto(nTprev+k).compressNucMask = compressBinaryImg(Lcytofin);
        end
    end
    nTprev = nTprev + nT;
    
   % save(outfile,'peaks','imgfiles');
    save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');
        
end
