
% info on where the images are and whidh position / timepoint to run

ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/masks_zcyto');
ilastikdirnuc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/Masks_z2');
imagedir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/rawimages_');%rawimages_
pos = 23;
tpt =1;
timegroup = 1;
chan = [1 2];

paramfile = 'setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN';

outfile = '3D_27tpts_xyz.mat';

%positions = [1 2 5 7 8 12 14 15 18 19 20 24 25 27 29 30]; %positions that were processed