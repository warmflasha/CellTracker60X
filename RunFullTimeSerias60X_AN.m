% run all the files supplied by Ilastik 
% need to be in the directory with a number of cyto and nuc masks, obtined
% from ilastik
% the output is a number of outall files, containing the peaks and colonies
% 

function [peaks,dims,colonies,imgfiles] = RunFullTimeSerias60X_AN(dir,zplane,direc,flag)


 %dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/AnalysisResults_Imaging1(earlyAugust2015)/IlastikDataFIles/');
 cd(dir);


[nums, filescyto]=folderFilesFromKeyword(dir,'Cyto');
[~, filesnuc]=folderFilesFromKeyword(dir,'Nuc');

for j = 1:length(nums)
        
ilastikfile = filesnuc(j).name;
ilastikfilecyto= filescyto(j).name;
    
[peaks,dims,colonies,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,nums(j),zplane,direc,flag);

%save(['Outfile_' num2str(nums(j)) ],'peaks','colonies','imgfiles');

end

end
% end up with a number of outall files for each position(all time points)