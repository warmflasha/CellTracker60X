function TrackGroupuCol60X(strdir,paramfiletrack)
global userParam;


ff = dir(strdir);
for k=1:size(ff,1)
    outfile = ff(k).name;
    load(outfile,'peaks','imgfiles','imgfilescyto');
    if ~exist('imgfiles','var') || ~exist('imgfilescyto','var')
        load(outfile,'peaks');
    end
    for k=1:length(peaks)
        if ~isempty(peaks{k})
            a = find(isnan(peaks{k}(:,1))) ;
            peaks{k}(a,:) = [];
        end
    end
    if ~exist('imgfiles','var') || ~exist('imgfilescyto','var')
        save(outfile,'peaks','-append');
        
    else
        save(outfile,'peaks','imgfiles','imgfilescyto','-append');
    end
    
    runTracker(outfile,paramfiletrack);
    userParam.colonygrouping = 130 ;   % 60X: 130
    cellsToDynColonies(outfile);
end
end