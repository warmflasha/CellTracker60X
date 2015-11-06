% run all the files supplied by Ilastik 
% need to be in the directory with a number of cyto and nuc masks, obtined
% from ilastik
% the output is a number of outall files, containing the peaks and colonies
% 

function RunFullTimeSerias60X_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

% pos should correspond to the name of the ilatikfile
for j = 1:length(currdir)
[peaks,dims,NucMasks,CytoMasks,colonies] = RunTimeSeries60XuColoniesAN(ilastikfile(j),ilastikfilecyto(j),pos,zplane,direc,flag);

save('Outfile_j','peaks','NucMasks','CytoMasks','colonies');
end
end
% end up with a number of outall files for each position(all time points)