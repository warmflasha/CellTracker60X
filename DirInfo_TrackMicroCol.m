
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

outfile = '3D_45tpts.mat';