% run all the files supplied by Ilastik 
% need to be in the directory with a number of cyto and nuc masks, obtined
% from ilastik
% the output is a number of outall files, containing the peaks and colonies
% 

function RunFullTimeSerias60X_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

% pos should correspond to the name of the ilatikfile
% dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/AnalysisResults_Imaging1(earlyAugust2015)/IlastikDataFiles/');
% cd(dir);
% 
 dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/AnalysisResults_Imaging1(earlyAugust2015)/IlastikDataFIles/');
% 
 %dirinfo = dir('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/AnalysisResults_Imaging1(earlyAugust2015)/IlastikDataFIles');
% numfiles = size(dirinfo,1);
% filename = cell(1,numfiles);
% for i = 1:numfiles
% if dirinfo(i).isdir == 0
% filename{i} = dirinfo(i).name;
% end
% end
[nums, filescyto]=folderFilesFromKeyword(dir,'Cyto');
[~, filesnuc]=folderFilesFromKeyword(dir,'Nuc');

for j = 1:length(nums)
        
ilastikfile = filesnuc(j).name;
ilastikfilecyto= filescyto(j).name;
    
[peaks,dims,NucMasks,CytoMasks,colonies] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag);

% check how the colonies array is returned
save(Outfile,'peaks','NucMasks','CytoMasks','colonies');
end
end
% end up with a number of outall files for each position(all time points)