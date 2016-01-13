function peaks = mrnapercells(nucleilist, stats, mrnamatfile, colonyno, zrange, channels)

% getting the center of nuclei in 3D
zstart = zrange(1);
zend = zrange(end);

nucleiinfo = cell(1,size(nucleilist,1));

for i = 1:size(nucleiinfo,2)
    
    object = nucleilist(i,:);
    objectmatchcol = find(~isnan(object));
    
    objectmatchz = objectmatchcol + zstart -1;
    objectsmatched = nucleilist(i,objectmatchcol);
    
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
peaks = nucleicenter;

for i = 1:numel(channels)
    
    mrnafile = strcat(mrnamatfile, '/', sprintf('ch%dallspots.mat', channels(i)));
    load(mrnafile);
    
    spots = spotinfomat(spotinfomat(:,1) == colonyno,:);
    
    spots(any(isnan(spots),2), :) = [];
    
    spotspos = [spots(:,3:4), spots(:,2)];
    
    mrnapercell = zeros(1, size(nucleicenter,1));
    
    %%
    
    for i = 1: size((spots),1)
        
        x0 = spotspos(i,1);
        y0 = spotspos(i,2);
        z0 = spotspos(i,3);
        
        %
        mydist = sqrt((nucleicenter(:,1) - x0).^2 + (nucleicenter(:,2) - y0).^2 + (nucleicenter(:,3) - z0).^2);
        
        [dist, celln] = min(mydist);
        
       
         mrnapercell(celln) = mrnapercell(celln)+floor(spots(i,5));
       
        
    end
    
    
    
    
    %%
     ncolpeaks = size(peaks,2) + 1;
     
   for i = 1:size(nucleicenter,1)
        peaks(i,ncolpeaks) = mrnapercell(i);
    end
end

%%