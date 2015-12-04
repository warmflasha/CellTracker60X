% run all the files supplied by Ilastik 
% need to be in the directory with a number of cyto and nuc masks, obtined
% from ilastik
% the output is a number of outall files, containing the peaks and colonies
% 

function [peaks,dims,imgfilescyto,imgfiles] = RunFullTimeSerias60X_AN(dir,zplane,direc,pos,dt,tg)% tg = time group

 % dir - has the ilastim .h5 files
 %dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/Nov12ImaginfResults/');
 %cd(dir);


[nums, filescyto]=folderFilesFromKeyword(dir,'Cyto','{0002}');% all cyto masks for the frame 0 ( four time groups)
[~, filesnuc]=folderFilesFromKeyword(dir,'Nuc','{0002}');% all nuc masks for the frame 0 ( four time groups)

for j = 1:length(nums)
        
ilastikfile = filesnuc(j).name;
ilastikfilecyto= filescyto(j).name;
    
[peaks,dims,imgfilescyto,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,dt,tg(j));
disp(['running position' num2str(nums(j)) ]);

save(['Outfile_0002_t' num2str(nums(j)) ],'peaks','dims','imgfiles','imgfilescyto');

end

end
% end up with a number of outall files for each position(all time points)