
function rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto,chanal)

% run all timepoints and save into one peaks
% main fnction
% pl = number of different ilastik files corresponding to the same position
% (either z planes or time groups, etc.)
% strnuc = keyword within the folder of nuc masks to separate the masks
% corresponding only to single position 'pos' ( e.g. give position may have masks
% separately for each x z plane)
%

[ilastikCytoAll] = FindPositionMasks(ilastikdircyto,pl,pos,strcyto);    % get the specific position ilastik masks (all z projections)
[ilastikNucAll] = FindPositionMasks(ilastikdirnuc,pl,pos,strnuc);

imgfilescyto = struct;
imgfiles = struct;

% need this chunck if working with multiple positions and timegroups
% ilastikCytoAll = ilastikCytoAll((pos*timegroups+1):(pos+1)*timegroups); % 1-4; 5-8;...
% ilastikNucAll = ilastikNucAll((pos*timegroups+1):(pos+1)*timegroups);
nTprev = 0;
%j = 1;%:length(ilastikNucAll)      % loop over the same position's time groups; not length(ilastikNucAll) is not the number of time groups
   
% read raw images

[imgsnuc_reader]   =  getrawimgfiles(imagedir1,pl, (pos-1),timegroup,chanal(1));        % get the raw images for that position and merge them into a 3d format
[imgscyto_reader] =   getrawimgfiles(imagedir2,pl, (pos-1),timegroup,chanal(1));    
nT = imgsnuc_reader{1}.getSizeT;                                                    % how many time point are within given time group

%nT = 81;% only for the february dataset (usable 82 timepoints)
for k =1:nT                                        % loop over time points within a given time group
        
    % read pnuc and pcyto separately from images
    [pnuc]=readmaskfiles1(ilastikNucAll,k,2); % readmaskfiles1(ilastikNucAll,k,lblN)    lblN = which ilastik label to use as cell nuc. and which as cyto
    [pcyto]=readmaskfiles1(ilastikCytoAll,k,2);
    
    for m = 1:size(imgscyto_reader,2) %
        planenuc = imgsnuc_reader{m}.getIndex(0,0, k - 1) + 1;
        inuc(:,:,m) = bfGetPlane(imgsnuc_reader{m},planenuc);
        planecyto = imgscyto_reader{m}.getIndex(0,0, k - 1) + 1;
        icyto(:,:,m) = bfGetPlane(imgscyto_reader{m},planecyto);
        
    end
         if pos == 39 || pos == 40 
             [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc(:,:,1),pcyto(:,:,1),inuc(:,:,1),icyto(:,:,1),timegroup,paramfile,paramfile3D);
         else
        [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc,pcyto,inuc,icyto,timegroup,paramfile,paramfile3D);
         end
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).NucMask = (Lnuc(:,:,round(size(Lnuc,3)/2))); %(Lnuc(:,:,1:size(Lnuc,3)))%round(size(Lnuc,3)/2) compressBinaryImg(Lnuc(:,:,3)) round(size(Lnuc,3)/2) compressBinaryImg
        imgfilescyto(nTprev+k).Cyto =(Lcytofin(:,:,round(size(Lcytofin,3)/2)));% (Lcytofin(:,:,1:size(Lcytofin,3)));%round(size(Lcytofin,3)/2)
        disp([k (pos-1)]);
    
end

% nTprev = nTprev + nT;% this is only if there are several time groups
% to put into the same peaks


save([ num2str(pos-1) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');



