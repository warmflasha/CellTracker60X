

function RenameIlastikOutputFiles(ilastikdir,positions,pl,str)
%
% function to rename the ilastik masks based on the position and zplanes
% number
% pl = number of z planes
% positions = vector with the position numbers (not necessaryly
% consecutive
% ilastikdir = directlry with ilastik masks, generated from headless mode
% (numbered consecutively)
% str = string that specified the basic name for the masks IlastikNucMasks_f0000_

[~, ff]=folderFilesFromKeyword(ilastikdir, ['CytoMasks3Dtile2_40x']);
% ff = dir(ilastikdir);
oldnames = cell(1, size(ff,2));
oldnames = {ff.name};

q = 1;

for j=1:length(positions)  %loop over positions names 
pos = positions(j); 

for i=q:(q+pl-1) 
torename = oldnames{i}; % needs to be string
newname = [str num2str(pos) '_z' num2str(i-q) '.h5'];
%disp([torename '   ' newname]);
movefile(torename, newname)
end
q = q+pl;

end
 
end
