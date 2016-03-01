
function rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto)

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
for j = 1%:length(ilastikNucAll)      % loop over the same position's time groups; not length(ilastikNucAll) is not the number of time groups
    
    nT = GetNumberTimePointsAN(imagedir,pos,timegroup);    % how many time point are within given time group
        
    
    for k = 1:nT                                                                       % loop over time points within a given time group
        [pnuc, imgsnuc]   =  readmaskfiles1(ilastikNucAll,imagedir, pos,k, j,chan(1));        % get the raw images for that position and merge them into a 3d format
        [pcyto, imgscyto] =  readmaskfiles1(ilastikCytoAll,imagedir2, pos,k, j,chan(2));
        for m = 1:pl
            img_now = imgsnuc{m}{1}{k,1};
            img_now_cyto = imgscyto{m}{1}{k,1};
            inuc(:,:,m) = img_now;
            icyto(:,:,m) = img_now_cyto;
        end
        
        [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc,pcyto, inuc,icyto,timegroup,paramfile,paramfile3D);
        
       
        peaks{nTprev+k} = outdat;
        
        imgfiles(nTprev+k).NucMask = Lnuc; %compressBinaryImg
        imgfilescyto(nTprev+k).NucMask = Lcytofin;
        
    end
    nTprev = nTprev + nT;
    save([ num2str(pos) '_' num2str(outfile)],'peaks','imgfiles','imgfilescyto');

end



