function [mrnapercell] = mrnapercell(nucleilist, stats, mrnamatfile, colonyno)

% getting the center of nuclei in 3D

nucleiinfo = cell(1,size(nucleilist,1));

for i = 1:size(nucleiinfo,2)
    
    object = nuclein1(i,:);
    objectmatchcol = find(~isnan(object));
    
    objectmatchz = objectmatchcol + zstart -1;
    objectsmatched = nuclein1(i,objectmatchcol);
    
    for match = 1:numel(objectsmatched)
        centroidmatch = stats{objectmatchz(match)}(objectsmatched(match)).Centroid;
        nucleiinfo{i}(match,:) = [centroidmatch objectmatchz(match)];
    end
    
    xcenter = 0.5*(max(nucleiinfo{i}(:,1))+min(nucleiinfo{i}(:,1)));
    ycenter = 0.5*(max(nucleiinfo{i}(:,2))+min(nucleiinfo{i}(:,2)));
    zcenter = 0.5*(max(nucleiinfo{i}(:,3))+min(nucleiinfo{i}(:,3)));
    
    nucleicenter(i,:) = [xcenter, ycenter, zcenter];
end

%%
% assigning mRNA's based on closest center in 3D
tic;
cellmrna3d = cell(1,size(nucleilist,1));

load(mrnamatfile);
spots = spotinfomat(spotinfomat(:,1) == colonyno,:);

filterspots = spots;
filterspots(any(isnan(spots),2), :) = [];

spotspos = [filterspots(:,3:4) filterspots(:,2)];
tic;
%%
cellmrna = cellmrna2d;



for i = 1: size((filterspots),1)
    
    x0 = spotspos(i,1);
    y0 = spotspos(i,2);
    z0 = spotspos(i,3);
    
    %
    mydist = sqrt((nucleicenter(:,1) - x0).^2 + (nucleicenter(:,2) - y0).^2 + (nucleicenter(:,3) - z0).^2);
    
    [dist, celln] = min(mydist);
    
    if(~isempty(cellmrna{celln}))
        
        nmrnacell = size(cellmrna{celln}{1},1);
        newrow = nmrnacell + 1;
    else
        newrow = 1;
    end
    
    cellmrna{celln}{1}(newrow,:) = filterspots(i,:);
    nucleimatch3d{i} = celln;
    
end


toc;

%%
for i = 1:size(cellmrna,2)
    
    if(isempty(cellmrna{i}))
        mrnapercell(i) = 0;
    else
        
        cmrna = cellmrna{i}{1};
        cmrna(any(isnan(cmrna),2),:) = [];
        
        nspots = size(cmrna,1);
        
        cellmrna{i}{1}(any(isnan(cellmrna{i}{1}),2),:) = [];
        
        clear mrna;
        
        for j = 1:nspots
            
            mrna(j) = cmrna(j,5);
        end
        
        mrnapercell(i) = floor(sum(mrna));
    end
end
