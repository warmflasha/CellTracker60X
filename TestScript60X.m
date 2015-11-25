
% test the watershed + ilastik analysis of 60X data


direccurr = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/AnalysisResults_Imaging1(earlyAugust2015)/IlastikDataFiles');
cd(direccurr);
direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/SingleCellSignalingAN_20150805_123245 PM');
direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/ANmicrocoloniesNov12(4)_20151116_102334AM');
%direc = ('/Volumes/data/Anastasia/LiveCellImagingGFPs4RFPh2b/SingleCellSignalingAN_20150805_123245 PM');
%direc = ('Z:\Anastasia\LiveCellImagingGFPs4RFPh2b/SingleCellSignalingAN_20150805_123245 PM');
flag = 0;
ilastikfile = 'NucMask0000batch.h5';
ilastikfilecyto = 'CytoMask0000batch.h5';
pos = 0;
zplane = 4;
img = 1;

%to run img = 1,first image


 [datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,img,dt,tg,imgs,imgs_nuc);


% to run all time points
% 
% [peaks,dims,imgfilescyto,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,dt);
%%
% 
[nums, files]=folderFilesFromKeyword(dir,'Outfile_');


for j=1:length(nums)

matfile = files(j).name;
load(matfile,'peaks');
% runTrackerEDS(matfile,'newTrackParam');
% load(matfile);

colonies=peaksToMicroColoniesAN(peaks);% for each time frame % here the colonies is a cell array : each cell is a colony object

 %save(['Outfile_' num2str(nums(j)) ],'peaks','dims','imgfiles','imgfilescyto','colonies','cells');
 %save('Outfile_15','peaks','dims','imgfiles','imgfilescyto','colonies','cells');




end
