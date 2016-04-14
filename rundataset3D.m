
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

[imgsnuc_reader]   =  getrawimgfiles(imagedir1,pl, pos,timegroup,chanal(1));        % get the raw images for that position and merge them into a 3d format
[imgscyto_reader] =   getrawimgfiles(imagedir2,pl, pos,timegroup,chanal(1));
nT = imgsnuc_reader{1}.getSizeT;                                                    % how many time point are within given time group

nT = 81;% only for the february dataset (usable 82 timepoints)
for k = 1:10;%1:nT                                                                   % loop over time points within a given time group
        
    % read pnuc and pcyto separately from images
    [pnuc]=readmaskfiles1(ilastikNucAll,k);
    [pcyto]=readmaskfiles1(ilastikCytoAll,k);
    
    for m = 1:size(imgscyto_reader,2)
        planenuc = imgsnuc_reader{m}.getIndex(0,0, k - 1) + 1;
        inuc(:,:,m) = bfGetPlane(imgsnuc_reader{m},planenuc);
        planecyto = imgscyto_reader{m}.getIndex(0,0, k - 1) + 1;
        icyto(:,:,m) = bfGetPlane(imgscyto_reader{m},planecyto);
        
    end
    
        [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc,pcyto, inuc,icyto,timegroup,paramfile,paramfile3D);
        
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).NucMask = (Lnuc(:,:,round(size(Lnuc,3)/2))); %compressBinaryImg(Lnuc(:,:,3)) round(size(Lnuc,3)/2) compressBinaryImg
        imgfilescyto(nTprev+k).Cyto = (Lcytofin(:,:,round(size(Lnuc,3)/2)));%round(size(Lcytofin,3)/2)
        disp([k pos]);
    
end

% nTprev = nTprev + nT;% this is only if there are several time groups
% to put into the same peaks


save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');



