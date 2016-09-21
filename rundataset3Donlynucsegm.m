function rundataset3Donlynucsegm(ilastikdirnuc,imagedir,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,chanal)

% run all timepoints and save into one peaks
% main fnction
% pl = number of different ilastik files corresponding to the same position
% (either z planes or time groups, etc.)
% strnuc = keyword within the folder of nuc masks to separate the masks
% corresponding only to single position 'pos' ( e.g. give position may have masks
% separately for each x z plane)
%

    % get the specific position ilastik masks (all z projections)
[ilastikNucAll] = FindPositionMasks(ilastikdirnuc,pl,pos,strnuc);
imgfiles = struct;


nTprev = 0;
   
% read raw images
[imgsnuc_reader]   =  getrawimgfiles(imagedir,pl, (pos-1),timegroup,chanal(1));   % get the raw images for that position and merge them into a 3d format
 
nT = imgsnuc_reader{1}.getSizeT;                                                    % how many time point are within given time group

for k =1  % loop over time points the fixed images have only one datapoint
        
    % read pnuc and pcyto separately from images
    [pnuc]=readmaskfiles1(ilastikNucAll,k,2); % readmaskfiles1(ilastikNucAll,k,lblN)    lblN = which ilastik label to use as cell nuc. and which as cyto
        
    for m = 1:size(imgsnuc_reader,2) %
        planenuc = imgsnuc_reader{m}.getIndex(0,0, k - 1) + 1;
        inuc(:,:,m) = bfGetPlane(imgsnuc_reader{m},planenuc);   % make sure take the correct plane m (not just the first one)
               
    end
         
        [outdat,Lnuc] = runmaskoneANdataonlynucsegm(pnuc,inuc,paramfile,paramfile3D,imagedir,pl,pos,timegroup);
        
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).NucMask = (Lnuc(:,:,round(size(Lnuc,3)/2))); %(Lnuc(:,:,1:size(Lnuc,3)))%round(size(Lnuc,3)/2) compressBinaryImg(Lnuc(:,:,3)) round(size(Lnuc,3)/2) compressBinaryImg
        disp([k (pos-1)]);
    
end

% nTprev = nTprev + nT;% this is only if there are several time groups
% to put into the same peaks


save([ num2str(pos-1) '_' num2str(outfile)],'peaks','imgfiles');



