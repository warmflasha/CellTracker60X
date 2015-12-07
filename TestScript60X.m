
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

