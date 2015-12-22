% run all the files supplied by Ilastik 
% need to be in the directory with a number of cyto and nuc masks, obtined
% from ilastik
% the output is a number of outall files, containing the peaks and colonies
% 

function [peaks,dims,imgfilescyto,imgfiles] = RunFullTimeSerias60X_AN(dir,zplane,direc,pos,dt,tg)

 % dir - has the ilastik masks, the  .h5 files
 % dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/Nov12ImaginfResults/');
 % zplane =  need to specify which image in z plane from the raw data to take 
 % direc =  directory with actual images from live-cell imaging experiment
 % pos = the number of frame processed TO DO: need to loop over those
 % dt =  data type: if 1 - means that all the deta is separate and the code
 % will not use bioformats to extract the images
 % if td = 0, the n the code will decompose gtouped data frames and extract
 % the correct images
 % tg = time group, this should be a vector, its length should be the
 % number of time froups the data was devided into upon saving
 % TO DO: 

[nums, filescyto]=folderFilesFromKeyword(dir,'Cyto',['{000' num2str(pos) '}']);% all cyto masks for the frame 0 ( four time groups)'Cyto','{0002}'['Outfile_000' num2str(pos) '_t']
[~, filesnuc]=folderFilesFromKeyword(dir,'Nuc',['{000' num2str(pos) '}']);% all nuc masks for the frame 0 ( four time groups)

for j = 1:length(nums)
        
ilastikfile = filesnuc(j).name;
ilastikfilecyto= filescyto(j).name;
    
[peaks,dims,imgfilescyto,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,dt,tg(j));
disp(['running position' num2str(nums(j)) ]);

save(['Outfile_000' num2str(pos) '_t' num2str(nums(j))],'peaks','dims','imgfiles','imgfilescyto');

end

end
% end up with a number of outall files for each position(all time points)