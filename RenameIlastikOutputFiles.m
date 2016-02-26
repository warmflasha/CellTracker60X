

function RenameIlastikOutputFiles(ilastikdir,positions)
%
% function to rename the ilastik masks based on the position and zplane
% number


ff = dir(ilastikdir);
names = {ff(~[ff.isdir]).ff};% get the names of files within directory

size(ff) % total number of masks in the dir

pl = 5;% number of z planes
q = 1;

for j=1:length(positions)  %loop over positions names 
pos = positions(j); 

for i=q:(q+pl-1)   

oldname = [ff names{i}];
if q ==1
newname = ([ ff 'position_' num2str(pos) '_z' num2str(i)]);
end
if q>1
newname = [ff 'position_' num2str(pos) '_z' num2str(i-5)];
end
movefile(oldname, newname)
end

q = q+5;
end
end

