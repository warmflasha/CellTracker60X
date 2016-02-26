
function rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir,pos,paramfile,timegroup,outfile,paramfile3D)

% run all timepoints and save into one peaks
% main fnction

[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikdircyto,['cytoframe' num2str(pos) '_']);   % get the specific position ilastik masks (all z projections)
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikdirnuc,['frame' num2str(pos) '_']);%

imgfilescyto = struct;
imgfiles = struct;

% need this chunck if working with multiple positions and timegroups
% ilastikCytoAll = ilastikCytoAll((pos*timegroups+1):(pos+1)*timegroups); % 1-4; 5-8;...
% ilastikNucAll = ilastikNucAll((pos*timegroups+1):(pos+1)*timegroups);
nTprev = 0;
for j = 1:length(ilastikNucAll)      % loop over the same position's time groups
    
    nT = GetNumberTimePointsAN(imagedir,pos,timegroup);    % how many time point are within given time group
    
    for k = 1:nT                                                                       % loop over time points within a given time group
        [pnuc, inuc]   =  readmaskfiles1(ilastikNucAll,imagedir, pos,k, j,chan(1));        % get the raw images for that position and merge them into a 3d format
        [pcyto, icyto] =  readmaskfiles1(ilastikCytoAll,imagedir, pos,k, j,chan(2));
        
        [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc,pcyto, inuc,icyto,timegroup,paramfile,paramfile3D);
        
        % [outdat, Lnuc,Lcytofin] = nucCytoIlastik2peaks(nuc_mask_all(:,:,k),cyto_mask_all(:,:,k),nuc_img,nuc_cyto,paramfile);%
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).NucMask = Lnuc; %compressBinaryImg
        imgfilescyto(nTprev+k).NucMask = Lcytofin;
        
    end
    nTprev = nTprev + nT;
    save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');

end



