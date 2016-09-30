
function [goodfile] = FindPositionMasks(ilastikdir,pl,pos,str)
% small function that find the ilastik masks corresponding to the selected
% position (pos) ; serches for these files within the directory ilastikdir
% much simpler than FilesFoldersFromKeyWord function
% str is the strin that specifies what to search for ( the same as str in
% the renameIlastikOutputFiles
% pl = number of planes or possible files correcponding to the same
% position ( could be either separate z planes of the same position or
% separate time groups of the same position (depends on the raw data type)

% output returned as a  cell array of strings 

%RenameIlastikOutputFiles

ff = dir(ilastikdir);
goodfile = cell(1,pl);
for k=1:size(ff,1)
    if ~isdir(ff(k).name)
        K = strfind(ff(k).name,[ str '_'  num2str(pos) '.' ]);% str '_'  num2str(pos) '_'
        if ~isempty(K)
            goodfile{k} = ff(k).name;
        end
    end
end
badinds = cellfun(@isempty,goodfile);
goodfile(badinds)=[];
end